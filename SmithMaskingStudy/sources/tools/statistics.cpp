#include "tools/statistics.h"

#include <iomanip>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <random>

#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/accumulators/statistics/weighted_covariance.hpp>
#include <boost/accumulators/statistics/weighted_kurtosis.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/weighted_skewness.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>

#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/distributions/normal.hpp>

#include "utils/math/math.h"
#include "utils/params.h"
#include "utils/console.h"
#include "tools/logger.h"
#include "NDFs/Microfacet.hpp"

using namespace boost::accumulators;
using namespace boost::math::statistics;
using namespace gdt;


int StatisticsTool::wHead = 23;
int StatisticsTool::wCell = 23;
int StatisticsTool::nColumns = 3;

StatisticsTool::StatisticsTool(const TriangleMesh* _mesh, ErrorStats error, int n_features) : mesh(_mesh), m_Error(error)
{
    switch (n_features)
    {
    case 5:
        compute5Statistics();
        break;
    case 25:
        compute25Statistics();
        break;
    default:
        computeStatistics();
        break;
    }
}

StatisticsTool::~StatisticsTool() { }

bool closeToEdge(const TriangleMesh* mesh, const gdt::vec3sc& point)
{
    return (mesh->bounds.closest_distance(point).x < Parameters::userParams.sideEffectParams.borderPercentage * mesh->bounds.span().x / 2.f
        || mesh->bounds.closest_distance(point).y < Parameters::userParams.sideEffectParams.borderPercentage * mesh->bounds.span().y / 2.f);
}

#define STATS_HEIGHT
#define STATS_AREA
#define STATS_NORMAL_CART
#define STATS_NORMAL_SPHE
#define STATS_SLOPE
#define STATS_CORR
#define STATS_CLUSTER

scal StatisticsTool::computeAnisotropy() const
{
    Discrete D(*mesh, nullptr, Parameters::userParams.sideEffectParams.borderPercentage);
    const scal phiStart = D.phiStart();
    const scal phiEnd = D.phiEnd();
    const scal phiSize = D.phiSize();
    const scal thetaStart = D.thetaStart();
    const scal thetaEnd = D.thetaEnd();
    const scal thetaSize = D.thetaSize();
    scal dPhi = (phiEnd - phiStart) / phiSize;
    scal dTheta = (thetaEnd - thetaStart) / thetaSize;
    scal integral = 0;
    for (int t = 0; t < thetaSize; ++t) {
        scal sin_theta = sin(thetaStart + t * dTheta);
        for (int p = 0; p < phiSize; ++p) {
            if (p == 0)                integral += (D[p + 1][t] - D[p    ][t]) * sin_theta;
            else if (p == phiSize - 1) integral += (D[p    ][t] - D[p - 1][t]) * sin_theta;
            else                       integral += (D[p + 1][t] - D[p - 1][t]) * sin_theta;
        }
    }
    integral *= dTheta / dPhi; // *= dTheta * dPhi / (2 * dPhi)
    return integral;
}

void StatisticsTool::compute5Statistics()
{
    const time_point time_start = std::chrono::high_resolution_clock::now();

    accumulator_set<ld, stats<tag::immediate_weighted_mean>, ld > accNormalTheta;
    std::vector<ld> thetas;
    std::vector<ld> heights;
    accumulator_set<ld, stats<tag::covariance<ld, tag::covariate1> > > accCovThetaHeight;

    // macrosurface area
    const scal border = 1. - Parameters::userParams.sideEffectParams.borderPercentage;
    const ld A = static_cast<ld>(mesh->bounds.span().x * border * mesh->bounds.span().y * border);

    // statistics for faces
    for (int faceID = 0; faceID < mesh->triangle_normal.size(); ++faceID)
    {
        const gdt::vec3i& idx = mesh->index[faceID];
        const gdt::vec3sc& faceN = mesh->triangle_normal[faceID];
        const gdt::vec3sc& v1 = mesh->vertex[idx[0]];
        const gdt::vec3sc& v2 = mesh->vertex[idx[1]];
        const gdt::vec3sc& v3 = mesh->vertex[idx[2]];
        const gdt::vec3sc barycenter = (v1 + v2 + v3) / (scal)3;

        if (closeToEdge(mesh, barycenter)) { continue; }

        const ld faceA = static_cast<ld>(mesh->area[faceID]) / A; // microfacet area relative to A

        const ld height = static_cast<ld>(barycenter.z);
        const ld theta = static_cast<ld>(faceN.theta());

        heights.push_back(height);
        accNormalTheta(theta, weight = faceA);
        thetas.push_back(theta);
        accCovThetaHeight(theta, covariate1 = height);
    }

    // statistics for vertices
    for (int vertexID = 0; vertexID < mesh->vertex.size(); ++vertexID)
    {
        const vec3sc& p = mesh->vertex[vertexID];
        const vec3sc& n = mesh->vertex_normal[vertexID];

        if (closeToEdge(mesh, p)) { continue; }

        const ld height = static_cast<ld>(p.z);

        heights.push_back(height);
        accCovThetaHeight(n.theta(), covariate1 = p.z);
    }

    // theta
    {
        m_nTPMean = vec2_ld(extract::weighted_mean(accNormalTheta), 0);

        auto const Q1 = thetas.size() / 4;
        auto const Q2 = thetas.size() / 2;
        auto const Q3 = Q1 + Q2;

        std::nth_element(thetas.begin(), thetas.begin() + Q1, thetas.end());
        std::nth_element(thetas.begin() + Q1 + 1, thetas.begin() + Q2, thetas.end());
        std::nth_element(thetas.begin() + Q2 + 1, thetas.begin() + Q3, thetas.end());

        m_TIQR = thetas[Q3] - thetas[Q1];
    }

    // heights
    {
        auto const Q1 = heights.size() / 4;
        auto const Q2 = heights.size() / 2;
        auto const Q3 = Q1 + Q2;

        std::nth_element(heights.begin(), heights.begin() + Q1, heights.end());
        std::nth_element(heights.begin() + Q1 + 1, heights.begin() + Q2, heights.end());
        std::nth_element(heights.begin() + Q2 + 1, heights.begin() + Q3, heights.end());

        m_heightIQR = heights[Q3] - heights[Q1];

        m_heightRange.lower = 0;
    }

    // anisotropy
    scal anisotropy = computeAnisotropy();

    // correlation
    {
        ld covThetaHeight = extract::covariance(accCovThetaHeight);
        m_THCorr = covThetaHeight / (m_nTPStd[0] * m_heightStd);
    }


    const time_point time_end = std::chrono::high_resolution_clock::now();
    const Duration time_duration = Duration(time_start, time_end);
    LOG_MESSAGE("5 features computed in " + time_duration.str());
    Console::out << "5 features computed in " << time_duration.str() << std::endl;
}

void StatisticsTool::compute25Statistics()
{
    const time_point time_start = std::chrono::high_resolution_clock::now();

    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_mean, tag::variance> > accHeight;
    std::vector<ld> heights;

    accumulator_set<ld, stats<tag::variance, tag::sum > > accArea;

    accumulator_set<ld, stats<tag::max, tag::immediate_weighted_mean, tag::weighted_variance, tag::weighted_skewness, tag::weighted_kurtosis >, ld > accNormalTheta;
    std::vector<ld> thetas;

    accumulator_set<ld, stats<tag::covariance<ld, tag::covariate1> > > accCovAreaTheta;
    accumulator_set<ld, stats<tag::covariance<ld, tag::covariate1> > > accCovThetaHeight;

    accumulator_set<ld, stats<tag::sum > > accLowHArea;
    std::vector<scal> centroids(3);
    std::vector<std::set<int>> sets = mesh->heightSeparation((int)centroids.size(), centroids); // each set contains the face's ids
    int min_idx = std::min_element(centroids.begin(), centroids.end()) - centroids.begin();
    std::set<int>& low_cluster = sets[min_idx];

    // macrosurface area
    const scal border = 1. - Parameters::userParams.sideEffectParams.borderPercentage;
    const ld A = static_cast<ld>(mesh->bounds.span().x * border * mesh->bounds.span().y * border);

    // statistics for faces
    for (int faceID = 0; faceID < mesh->triangle_normal.size(); ++faceID)
    {
        const gdt::vec3i& idx = mesh->index[faceID];
        const gdt::vec3sc& faceN = mesh->triangle_normal[faceID];
        const gdt::vec3sc& v1 = mesh->vertex[idx[0]];
        const gdt::vec3sc& v2 = mesh->vertex[idx[1]];
        const gdt::vec3sc& v3 = mesh->vertex[idx[2]];
        const gdt::vec3sc barycenter = (v1 + v2 + v3) / (scal)3;

        if (closeToEdge(mesh, barycenter)) { continue; }

        const ld faceA = static_cast<ld>(mesh->area[faceID]) / A; // microfacet area relative to A

        const ld height = static_cast<ld>(barycenter.z);
        const ld theta = static_cast<ld>(faceN.theta());

        const ld dh = static_cast<ld>(std::max({ abs(v1.z - v2.z), abs(v2.z - v3.z), abs(v3.z - v1.z) }));
        accHeight(height);
        heights.push_back(height);

        accArea(faceA);

        accNormalTheta(theta, weight = faceA);
        const ld phi = static_cast<ld>(faceN.phi());
        thetas.push_back(theta);

        accCovAreaTheta(faceA, covariate1 = theta);
        accCovThetaHeight(theta, covariate1 = height);

        if (low_cluster.contains(faceID)) {
            const ld cosTheta = static_cast<ld>(dot(faceN, mesh->meso_normal));
            accLowHArea(faceA * cosTheta);
        }
    }

    // statistics for vertices
    for (int vertexID = 0; vertexID < mesh->vertex.size(); ++vertexID)
    {
        const vec3sc& p = mesh->vertex[vertexID];
        const vec3sc& n = mesh->vertex_normal[vertexID];

        if (closeToEdge(mesh, p)) { continue; }

        const ld height = static_cast<ld>(p.z);

        accHeight(height);
        heights.push_back(height);
        accCovThetaHeight(n.theta(), covariate1 = p.z);
    }

    // theta
    {
        m_nTPMean = vec2_ld(extract::weighted_mean(accNormalTheta), 0);
        vec2_ld nTPVar(extract::weighted_variance(accNormalTheta), 0);
        m_nTPStd = { sqrt(nTPVar.x), sqrt(nTPVar.y) };
        m_nTPRanges.upper = vec2_ld(extract::max(accNormalTheta), 0);

        m_TCV = m_nTPStd[0] / m_nTPMean[0];
        m_TMAD = median_absolute_deviation(thetas.begin(), thetas.end());

        auto const Q1 = thetas.size() / 4;
        auto const Q2 = thetas.size() / 2;
        auto const Q3 = Q1 + Q2;

        std::nth_element(thetas.begin(), thetas.begin() + Q1, thetas.end());
        std::nth_element(thetas.begin() + Q1 + 1, thetas.begin() + Q2, thetas.end());
        std::nth_element(thetas.begin() + Q2 + 1, thetas.begin() + Q3, thetas.end());

        m_TQ1 = thetas[Q1];
        m_TMedian = thetas[Q2];
        m_TQ3 = thetas[Q3];
        m_TIQR = m_TQ3 - m_TQ1;
        m_TQCD = m_TIQR / (m_TQ3 + m_TQ1);
    }

    // shape parameters
    {
        m_skewness = extract::weighted_skewness(accNormalTheta);
        m_kurtosis = extract::weighted_kurtosis(accNormalTheta);
    }

    // areas
    {
        m_totalA = extract::sum(accArea);
        m_areaStd = sqrt(extract::variance(accArea));
    }

    // cluster
    {
        m_lowH = extract::sum(accLowHArea); // faceA are already relatives to A, just as if the patch was 1 unit
    }

    // heights
    {
        m_heightRange.lower = extract::min(accHeight);
        m_heightRange.upper = extract::max(accHeight) - m_heightRange.lower;
        m_heightMean = extract::mean(accHeight) - m_heightRange.lower;
        m_heightStd = sqrt(extract::variance(accHeight));

        m_heightCV = m_heightStd / m_heightMean;
        m_heightMAD = median_absolute_deviation(heights.begin(), heights.end()); // median(|h - median(heights)|) = median(|(h-min) - median(heights-min)|) = median(|h - min - median(heights) + min|)

        auto const Q1 = heights.size() / 4;
        auto const Q2 = heights.size() / 2;
        auto const Q3 = Q1 + Q2;

        std::nth_element(heights.begin(), heights.begin() + Q1, heights.end());
        std::nth_element(heights.begin() + Q1 + 1, heights.begin() + Q2, heights.end());
        std::nth_element(heights.begin() + Q2 + 1, heights.begin() + Q3, heights.end());

        m_heightQ1 = heights[Q1] - m_heightRange.lower;
        m_heightMedian = heights[Q2] - m_heightRange.lower;
        m_heightQ3 = heights[Q3] - m_heightRange.lower;
        m_heightIQR = m_heightQ3 - m_heightQ1;
        m_heightQCD = m_heightIQR / (m_heightQ3 + m_heightQ1);

        m_heightRange.lower = 0;
    }

    // anisotropy
    {
        scal anisotropy = computeAnisotropy();
    }

    // correlations
    {
        ld covThetaHeight = extract::covariance(accCovThetaHeight);
        m_THCorr = covThetaHeight / (m_nTPStd[0] * m_heightStd);

        ld covAreaTheta = extract::covariance(accCovAreaTheta);
        m_ATCorr = covAreaTheta / (m_areaStd * m_nTPStd[0]);
    }

    const time_point time_end = std::chrono::high_resolution_clock::now();
    const Duration time_duration = Duration(time_start, time_end);
    LOG_MESSAGE("25 features computed in " + time_duration.str());
    Console::out << "25 features computed in " << time_duration.str() << std::endl;
}


void StatisticsTool::computeStatistics()
{
#ifdef STATS_HEIGHT
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_mean, tag::variance> > accHeight;
    accumulator_set<ld, stats<tag::max> > accDeltaHeight;
    std::vector<ld> heights;
#endif // STATS_HEIGHT

#ifdef STATS_AREA
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_mean, tag::variance, tag::sum > > accArea;
#endif // STATS_AREA

#ifdef STATS_NORMAL_CART
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_weighted_mean, tag::weighted_variance >, ld > accNormalX;
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_weighted_mean, tag::weighted_variance >, ld > accNormalY;
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_weighted_mean, tag::weighted_variance >, ld > accNormalZ;
    accumulator_set<ld, stats<tag::weighted_covariance<ld, tag::covariate1> >, ld > accCovNormalXY;
#endif // STATS_NORMAL_CART

#ifdef STATS_NORMAL_SPHE
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_weighted_mean, tag::weighted_variance, tag::weighted_skewness, tag::weighted_kurtosis >, ld > accNormalTheta;
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_weighted_mean, tag::weighted_variance >, ld > accNormalPhi;
    accumulator_set<ld, stats<tag::weighted_covariance<ld, tag::covariate1> >, ld > accCovNormalTP;
    std::vector<ld> thetas;
#endif // STATS_NORMAL_SPHE

#ifdef STATS_SLOPE
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_weighted_mean, tag::weighted_variance >, ld > accSlopeX;
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_weighted_mean, tag::weighted_variance >, ld > accSlopeY;
    accumulator_set<ld, stats<tag::weighted_covariance<ld, tag::covariate1> >, ld > accCovSlope;
#endif // STATS_SLOPE

#ifdef STATS_CORR
    accumulator_set<ld, stats<tag::covariance<ld, tag::covariate1> > > accCovAreaTheta;
    accumulator_set<ld, stats<tag::covariance<ld, tag::covariate1> > > accCovThetaHeight;
#endif // STATS_CORR

#ifdef STATS_CLUSTER
    accumulator_set<ld, stats<tag::sum > > accLowHArea;
    std::vector<scal> centroids(3);
    std::vector<std::set<int>> sets = mesh->heightSeparation((int)centroids.size(), centroids); // each set contains the face's ids
    int min_idx = std::min_element(centroids.begin(), centroids.end()) - centroids.begin();
    Console::info << min_idx << std::endl;
    std::set<int>& low_cluster = sets[min_idx];
#endif // STATS_CLUSTER

    // macrosurface area
    const scal border = 1. - Parameters::userParams.sideEffectParams.borderPercentage;
    const ld A = static_cast<ld>(mesh->bounds.span().x * border * mesh->bounds.span().y * border);

    // statistics for faces
    for (int faceID = 0; faceID < mesh->triangle_normal.size(); ++faceID)
    {
        const gdt::vec3i& idx = mesh->index[faceID];
        const gdt::vec3sc& faceN = mesh->triangle_normal[faceID];
        const gdt::vec3sc& v1 = mesh->vertex[idx[0]];
        const gdt::vec3sc& v2 = mesh->vertex[idx[1]];
        const gdt::vec3sc& v3 = mesh->vertex[idx[2]];
        const gdt::vec3sc barycenter = (v1 + v2 + v3) / (scal)3;

        if (closeToEdge(mesh, barycenter)) { continue; }

        const ld faceA    = static_cast<ld>(mesh->area[faceID]) / A; // microfacet area relative to A

#if defined(STATS_HEIGHT) || defined(STATS_CORR)
        const ld height   = static_cast<ld>(barycenter.z);
#endif // STATS_HEIGHT || STATS_CORR
#if defined(STATS_NORMAL_SPHE) || defined(STATS_CORR)
        const ld theta    = static_cast<ld>(faceN.theta());
#endif // STATS_NORMAL_SPHE || STATS_CORR
#if defined(STATS_NORMAL_CART) || defined(STATS_SLOPE)
        const ld Nx       = static_cast<ld>(faceN.x);
        const ld Ny       = static_cast<ld>(faceN.y);
        const ld Nz       = static_cast<ld>(faceN.z);
#endif // STATS_NORMAL_CART || STATS_SLOPE

#ifdef STATS_HEIGHT
        const ld dh = static_cast<ld>(std::max({ abs(v1.z - v2.z), abs(v2.z - v3.z), abs(v3.z - v1.z) }));
        accHeight(height);
        accDeltaHeight(dh);
        heights.push_back(height);
#endif // STATS_HEIGHT

#ifdef STATS_AREA
        accArea(faceA);
#endif // STATS_AREA

#ifdef STATS_NORMAL_CART
        accNormalX(Nx, weight = faceA);
        accNormalY(Ny, weight = faceA);
        accNormalZ(Nz, weight = faceA);
        accCovNormalXY(Nx, weight = faceA, covariate1 = Ny);
#endif // STATS_NORMAL_CART

#ifdef STATS_NORMAL_SPHE
        accNormalTheta(theta, weight = faceA);
        const ld phi = static_cast<ld>(faceN.phi());
        accNormalPhi(phi, weight = faceA);
        accCovNormalTP(theta, weight = faceA, covariate1 = phi);
        thetas.push_back(theta);
#endif // STATS_NORMAL_SPHE

#ifdef STATS_SLOPE
        accSlopeX(-Nx / Nz, weight = faceA);
        accSlopeY(-Ny / Nz, weight = faceA);
        accCovSlope(-Nx / Nz, weight = faceA, covariate1 = -Ny / Nz);
#endif // STATS_SLOPE

#ifdef STATS_CORR
        accCovAreaTheta(faceA, covariate1 = theta);
        accCovThetaHeight(theta, covariate1 = height);
#endif // STATS_CORR

#ifdef STATS_CLUSTER
        if (low_cluster.contains(faceID)) {
            const ld cosTheta = static_cast<ld>(dot(faceN, mesh->meso_normal));
            accLowHArea(faceA * cosTheta);
        }
#endif // STATS_CLUSTER
    }
   
    // statistics for vertices
#if defined(STATS_HEIGHT) || defined(STATS_CORR)
    for (int vertexID = 0; vertexID < mesh->vertex.size(); ++vertexID)
    {
        const vec3sc& p = mesh->vertex[vertexID];
        const vec3sc& n = mesh->vertex_normal[vertexID];

        if (closeToEdge(mesh, p)) { continue; }

        const ld height = static_cast<ld>(p.z);

#ifdef STATS_HEIGHT
        accHeight(height);
        heights.push_back(height);
#endif // STATS_HEIGHT
#ifdef STATS_CORR
        accCovThetaHeight(n.theta(), covariate1 = p.z);
#endif // STATS_CORR
    }
#endif // STATS_HEIGHT || STATS_CORR

#ifdef STATS_NORMAL_CART
    {
        m_nXYZMean = vec3_ld(
            extract::weighted_mean(accNormalX),
            extract::weighted_mean(accNormalY),
            extract::weighted_mean(accNormalZ));
        vec3_ld nXYZVar(
            extract::weighted_variance(accNormalX),
            extract::weighted_variance(accNormalY),
            extract::weighted_variance(accNormalZ));
        m_nXYZStd = { sqrt(nXYZVar.x), sqrt(nXYZVar.y), sqrt(nXYZVar.z) };
        ld nXYCov = extract::weighted_covariance(accCovNormalXY);
        m_nXYCorr = nXYCov / (m_nXYZStd.x * m_nXYZStd.y);
        m_nXYRanges.lower = vec2_ld(
            extract::min(accNormalX),
            extract::min(accNormalY));
        m_nXYRanges.upper = vec2_ld(
            extract::max(accNormalX),
            extract::max(accNormalY));
    }
#endif // STATS_NORMAL_CART

#ifdef STATS_NORMAL_SPHE
    {
        m_nTPMean = vec2_ld(
            extract::weighted_mean(accNormalTheta),
            extract::weighted_mean(accNormalPhi));
        vec2_ld nTPVar(
            extract::weighted_variance(accNormalTheta),
            extract::weighted_variance(accNormalPhi));
        m_nTPStd = { sqrt(nTPVar.x), sqrt(nTPVar.y) };
        ld nTPCov = extract::weighted_covariance(accCovNormalTP);
        m_nTPCorr = nTPCov / (m_nTPStd.x * m_nTPStd.y);
        m_nTPRanges.lower = vec2_ld(
            extract::min(accNormalTheta),
            extract::min(accNormalPhi));
        m_nTPRanges.upper = vec2_ld(
            extract::max(accNormalTheta),
            extract::max(accNormalPhi));
    }

    // theta
    {
        m_TCV = m_nTPStd[0] / m_nTPMean[0];
        m_TMAD = median_absolute_deviation(thetas.begin(), thetas.end());

        auto const Q1 = thetas.size() / 4;
        auto const Q2 = thetas.size() / 2;
        auto const Q3 = Q1 + Q2;

        std::nth_element(thetas.begin(), thetas.begin() + Q1, thetas.end());
        std::nth_element(thetas.begin() + Q1 + 1, thetas.begin() + Q2, thetas.end());
        std::nth_element(thetas.begin() + Q2 + 1, thetas.begin() + Q3, thetas.end());

        m_TQ1     = thetas[Q1];
        m_TMedian = thetas[Q2];
        m_TQ3     = thetas[Q3];
        m_TIQR    = m_TQ3 - m_TQ1;
        m_TQCD    = m_TIQR / (m_TQ3 + m_TQ1);
    }

    // shape parameters
    {
        m_skewness = extract::weighted_skewness(accNormalTheta);
        m_kurtosis = extract::weighted_kurtosis(accNormalTheta);
    }
#endif //  STATS_THETA

#ifdef STATS_SLOPE
    // slopes
    {
        m_slopeMean = vec2_ld(
            extract::weighted_mean(accSlopeX),
            extract::weighted_mean(accSlopeY));
        vec2_ld sVar(
            extract::weighted_variance(accSlopeX),
            extract::weighted_variance(accSlopeY));
        m_slopeStd = { sqrt(sVar.x), sqrt(sVar.y) };
        ld sCov = extract::weighted_covariance(accCovSlope);
        m_slopeCorr = sCov / (m_slopeStd.x * m_slopeStd.y);
        m_slopeRanges.lower = vec2_ld(
            extract::min(accSlopeX),
            extract::min(accSlopeY));
        m_slopeRanges.upper = vec2_ld(
            extract::max(accSlopeX),
            extract::max(accSlopeY));
    }
#endif // STATS_SLOPE

#ifdef STATS_AREA
    {
        m_totalA = extract::sum(accArea);
        m_areaMean = extract::mean(accArea);
        m_areaStd = sqrt(extract::variance(accArea));
        m_areaRange.lower = extract::min(accArea);
        m_areaRange.upper = extract::max(accArea);
    }
#endif // STATS_AREA

#ifdef STATS_CLUSTER
    {
        m_lowH = extract::sum(accLowHArea); // faceA are already relatives to A, just as if the patch was 1 unit
    }
#endif // STATS_CLUSTER

#ifdef STATS_HEIGHT
    {
        m_heightRange.lower = extract::min(accHeight);
        m_heightRange.upper = extract::max(accHeight) - m_heightRange.lower;
        m_heightMean = extract::mean(accHeight) - m_heightRange.lower;
        m_heightStd = sqrt(extract::variance(accHeight));
        m_maxDeltaHeight = extract::max(accDeltaHeight);

        m_heightCV  = m_heightStd / m_heightMean;
        m_heightMAD = median_absolute_deviation(heights.begin(), heights.end()); // median(|h - median(heights)|) = median(|(h-min) - median(heights-min)|) = median(|h - min - median(heights) + min|)

        auto const Q1 = heights.size() / 4;
        auto const Q2 = heights.size() / 2;
        auto const Q3 = Q1 + Q2;

        std::nth_element(heights.begin(),          heights.begin() + Q1, heights.end());
        std::nth_element(heights.begin() + Q1 + 1, heights.begin() + Q2, heights.end());
        std::nth_element(heights.begin() + Q2 + 1, heights.begin() + Q3, heights.end());

        m_heightQ1     = heights[Q1] - m_heightRange.lower;
        m_heightMedian = heights[Q2] - m_heightRange.lower;
        m_heightQ3     = heights[Q3] - m_heightRange.lower;
        m_heightIQR    = m_heightQ3 - m_heightQ1;
        m_heightQCD    = m_heightIQR / (m_heightQ3 + m_heightQ1);

        m_heightRange.lower = 0;
    }

    // Scores
    {
        // Height ratio
        m_heightSums = std::accumulate(
            mesh->vertex.begin(), mesh->vertex.end(), gdt::vec2i(0, 0),
            [this](gdt::vec2i s, gdt::vec3sc p) {
                const ld height = static_cast<ld>(p.z);
                return height < m_heightMean ? gdt::vec2i(s[0] + 1, s[1]) : gdt::vec2i(s[0], s[1] + 1);
            }
        );

        // Slope score


        // Anisotropy score
    }
#endif // STATS_HEIGHT

#ifdef STATS_CORR
    {
        ld covThetaHeight = extract::covariance(accCovThetaHeight);
        m_THCorr = covThetaHeight / (m_nTPStd[0] * m_heightStd);

        ld covAreaTheta = extract::covariance(accCovAreaTheta);
        m_ATCorr = covAreaTheta / (m_areaStd * m_nTPStd[0]);
    }
#endif // STATS_CORR

    LOG_STATS(*this);
}

// EXPORT
void StatisticsTool::CSVHeader(csv::CSVWriter* writer) {
    std::vector<csv::elem> headers;

    headers.push_back({ csv::elem::Tag::STRING, "index" });

    headers.push_back({ csv::elem::Tag::STRING, "error" });
    headers.push_back({ csv::elem::Tag::STRING, "meanError" });
    headers.push_back({ csv::elem::Tag::STRING, "stdError" });
    headers.push_back({ csv::elem::Tag::STRING, "minError" });
    headers.push_back({ csv::elem::Tag::STRING, "maxError" });

    headers.push_back({ csv::elem::Tag::STRING, "heightSumInf" });
    headers.push_back({ csv::elem::Tag::STRING, "heightSumSup" });
    headers.push_back({ csv::elem::Tag::STRING, "slopeScore" });
    headers.push_back({ csv::elem::Tag::STRING, "anisotropyScore" });

    headers.push_back({ csv::elem::Tag::STRING, "minH" });
    headers.push_back({ csv::elem::Tag::STRING, "maxH" });
    headers.push_back({ csv::elem::Tag::STRING, "meanH" });
    headers.push_back({ csv::elem::Tag::STRING, "stdH" });
    headers.push_back({ csv::elem::Tag::STRING, "maxDeltaH" });

    headers.push_back({ csv::elem::Tag::STRING, "cvH" });
    headers.push_back({ csv::elem::Tag::STRING, "madH" });
    headers.push_back({ csv::elem::Tag::STRING, "q1H" });
    headers.push_back({ csv::elem::Tag::STRING, "q2H" });
    headers.push_back({ csv::elem::Tag::STRING, "q3H" });
    headers.push_back({ csv::elem::Tag::STRING, "iqrH" });
    headers.push_back({ csv::elem::Tag::STRING, "qcdH" });

    headers.push_back({ csv::elem::Tag::STRING, "totalA" });
    headers.push_back({ csv::elem::Tag::STRING, "minA" });
    headers.push_back({ csv::elem::Tag::STRING, "maxA" });
    headers.push_back({ csv::elem::Tag::STRING, "meanA" });
    headers.push_back({ csv::elem::Tag::STRING, "stdA" });

    headers.push_back({ csv::elem::Tag::STRING, "lowH" });

    headers.push_back({ csv::elem::Tag::STRING, "minNX" });
    headers.push_back({ csv::elem::Tag::STRING, "maxNX" });
    headers.push_back({ csv::elem::Tag::STRING, "meanNX" });
    headers.push_back({ csv::elem::Tag::STRING, "stdNX" });

    headers.push_back({ csv::elem::Tag::STRING, "minNY" });
    headers.push_back({ csv::elem::Tag::STRING, "maxNY" });
    headers.push_back({ csv::elem::Tag::STRING, "meanNY" });
    headers.push_back({ csv::elem::Tag::STRING, "stdNY" });

    headers.push_back({ csv::elem::Tag::STRING, "minT" });
    headers.push_back({ csv::elem::Tag::STRING, "maxT" });
    headers.push_back({ csv::elem::Tag::STRING, "meanT" });
    headers.push_back({ csv::elem::Tag::STRING, "stdT" });

    headers.push_back({ csv::elem::Tag::STRING, "cvT" });
    headers.push_back({ csv::elem::Tag::STRING, "madT" });
    headers.push_back({ csv::elem::Tag::STRING, "q1T" });
    headers.push_back({ csv::elem::Tag::STRING, "q2T" });
    headers.push_back({ csv::elem::Tag::STRING, "q3T" });
    headers.push_back({ csv::elem::Tag::STRING, "iqrT" });
    headers.push_back({ csv::elem::Tag::STRING, "qcdT" });

    headers.push_back({ csv::elem::Tag::STRING, "skewness" });
    headers.push_back({ csv::elem::Tag::STRING, "kurtosis" });

    headers.push_back({ csv::elem::Tag::STRING, "corrXY" });
    headers.push_back({ csv::elem::Tag::STRING, "corrAT" });
    headers.push_back({ csv::elem::Tag::STRING, "corrTH" });

    writer->writeRow(headers);
}

void StatisticsTool::toCSV(csv::CSVWriter* writer) {
    std::vector<csv::elem> stats;

    stats.push_back({ csv::elem::Tag::STRING, mesh->name }); // index

    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_Error.integral }); // Error
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_Error.mean }); // mean(Error)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_Error.std }); // std(Error)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_Error.min }); // min(Error)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_Error.max }); // max(Error)

    stats.push_back({ csv::elem::Tag::INT, m_heightSums[0] }); // Height sum inf
    stats.push_back({ csv::elem::Tag::INT, m_heightSums[1] }); // Height sum sup
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_slopeScore }); // Slope score
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_anisotropyScore }); // Anisotropy score

    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_heightRange.lower }); // min(H)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_heightRange.upper }); // max(H)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_heightMean }); // mean(H)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_heightStd }); // std(H)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_maxDeltaHeight }); // maxDeltaH

    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_heightCV }); // cv(H)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_heightMAD }); // mad(H)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_heightQ1 }); // q1(H)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_heightMedian }); // q2(H)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_heightQ3 }); // q3(H)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_heightIQR }); // iqr(H)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_heightQCD }); // qcd(H)

    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_totalA }); // total(A)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_areaRange.lower }); // min(A)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_areaRange.upper }); // max(A)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_areaMean }); // mean(A)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_areaStd }); // std(A)

    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_lowH }); // percentage of projected lowest cluster

    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nXYRanges.lower.x }); // min(N.x)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nXYRanges.upper.x }); // max(N.x)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nXYZMean.x }); // mean(N.x)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nXYZStd.x }); // std(N.x)

    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nXYRanges.lower.y }); // min(N.y)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nXYRanges.upper.y }); // max(N.y)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nXYZMean.y }); // mean(N.y)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nXYZStd.y }); // std(N.y)

    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nTPRanges.lower.x }); // min(T)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nTPRanges.upper.x }); // max(T)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nTPMean.x }); // mean(T)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nTPStd.x }); // std(T)

    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_TCV }); // cv(T)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_TMAD }); // mad(T)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_TQ1 }); // q1(T)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_TMedian }); // q2(T)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_TQ3 }); // q3(T)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_TIQR }); // iqr(T)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_TQCD }); // qcd(T)

    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_skewness }); // skewness
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_kurtosis }); // kurtosis

    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_nXYCorr }); // Corr(x ; y)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_ATCorr }); // Corr(A ; T)
    stats.push_back({ csv::elem::Tag::LONG_DOUBLE, m_THCorr }); // Corr(T ; H)

    writer->writeRow(stats);
}


void StatisticsTool::print() {
    Console::out << *this << std::endl;
}


std::ostream& operator<<(std::ostream& output, const StatisticsTool& stats) {
    output << Console::timePad << "Statistics: "  << stats.mesh->name << std::endl;
    output << Console::timePad << "-- Heights " << std::endl;
    output << Console::timePad << "   - m_heightRange: " << stats.m_heightRange << std::endl;
    output << Console::timePad << "   - m_heightMean: " << stats.m_heightMean << std::endl; // mean(H)
    output << Console::timePad << "   - m_heightStd: " << stats.m_heightStd << std::endl; // std(H)
    output << Console::timePad << "   - m_maxDeltaHeight: " << stats.m_maxDeltaHeight << std::endl; // maxDeltaH
    output << Console::timePad << "-- Areas " << std::endl;
    output << Console::timePad << "  -- totalA: " << static_cast<ld>(stats.mesh->surfaceArea) / static_cast<ld>(stats.mesh->macroArea) << std::endl; // total(A)
    output << Console::timePad << "  -- m_areaRange: " << stats.m_areaRange << std::endl; // min(A)
    output << Console::timePad << "  -- m_areaMean: " << stats.m_areaMean << std::endl; // mean(A)
    output << Console::timePad << "  -- m_areaStd: " << stats.m_areaStd << std::endl; // std(A)
    output << Console::timePad << "-- Clusters " << std::endl;
    output << Console::timePad << "  -- m_lowH: " << stats.m_lowH << std::endl; // lowH
    output << Console::timePad << "-- Normals " << std::endl;
    output << Console::timePad << "  -- m_nXYRanges: " << stats.m_nXYRanges << std::endl; // max(N.x)
    output << Console::timePad << "  -- m_nXYZMean " << stats.m_nXYZMean << std::endl; // mean(N.x)
    output << Console::timePad << "  -- m_nXYZStd: " << stats.m_nXYZStd << std::endl; // std(N.x)
    output << Console::timePad << "  -- m_nTPRanges: " << stats.m_nTPRanges << std::endl; // max(T)
    output << Console::timePad << "  -- m_nTPMean " << stats.m_nTPMean << std::endl; // mean(T)
    output << Console::timePad << "  -- m_nTPStd: " << stats.m_nTPStd << std::endl; // std(T)
    output << Console::timePad << "-- Shape " << std::endl;
    output << Console::timePad << "  -- m_skewness: " << stats.m_skewness << std::endl; // skewness
    output << Console::timePad << "  -- m_kurtosis: " << stats.m_kurtosis << std::endl; // kurtosis
    output << Console::timePad << "-- Correlations " << std::endl;
    output << Console::timePad << "  -- m_nXYCorr: " << stats.m_nXYCorr << std::endl; // Corr(x ; y)
    output << Console::timePad << "  -- m_ATCorr: " << stats.m_ATCorr << std::endl; // Corr(A ; T)
    output << Console::timePad << "  -- m_THCorr: " << stats.m_THCorr << std::endl; // Corr(T ; H)

    return output;
}