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

StatisticsTool::StatisticsTool(const TriangleMesh* _mesh, scal error, int n_features) : mesh(_mesh), m_error(error)
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
    return (mesh->bounds.closest_distance(point).x < Parameters::get()->currentParams()->sideEffectParams.borderPercentage * mesh->bounds.span().x / 2.f
        || mesh->bounds.closest_distance(point).y < Parameters::get()->currentParams()->sideEffectParams.borderPercentage * mesh->bounds.span().y / 2.f);
}

#define STATS_PHI
#define STATS_THETA
#define STATS_AREA
#define STATS_HEIGHT
#define STATS_CORR

scal StatisticsTool::computeAnisotropy() const
{
    Discrete D(*mesh, nullptr, Parameters::get()->currentParams()->sideEffectParams.borderPercentage);
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

    // Compute for: m_phiAnisotropy, m_thetaIQR, m_areaTotal, m_heightStd, m_heightCb

    m_phiAnisotropy = computeAnisotropy();                          // for m_phiAnisotropy
    std::vector<ld> thetas;                                         // for m_thetaIQR
    accumulator_set<ld, stats<tag::sum > > accArea;                 // for m_areaTotal
    accumulator_set<ld, stats<tag::variance> > accHeight;           // for m_heightStd
    std::vector<ld> heights;                                        // 
    accumulator_set<ld, stats<tag::sum > > accHeightCb;             // for m_heightCb
    std::vector<scal> centroids(3);                                 //
    std::vector<std::set<int>> sets = mesh->heightSeparation(       //
        (int)centroids.size(),                                      //
        centroids                                                   //
    );                                                              //
    int min_idx = std::min_element(                                 //
        centroids.begin(),                                          //
        centroids.end()                                             //
    ) - centroids.begin();                                          //
    std::set<int>& low_cluster = sets[min_idx];                     //

    // macrosurface area
    const scal border = 1. - Parameters::get()->currentParams()->sideEffectParams.borderPercentage;
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

        thetas.push_back(theta);
        accArea(faceA);
        accHeight(height);
        heights.push_back(height);
        if (low_cluster.contains(faceID)) {
            const ld cosTheta = static_cast<ld>(dot(faceN, mesh->meso_normal));
            accHeightCb(faceA * cosTheta);
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
    }

    // m_thetaIQR
    {
        auto const Q1 = thetas.size() / 4;
        auto const Q2 = thetas.size() / 2;
        auto const Q3 = Q1 + Q2;

        std::nth_element(thetas.begin(), thetas.begin() + Q1, thetas.end());
        std::nth_element(thetas.begin() + Q1 + 1, thetas.begin() + Q2, thetas.end());
        std::nth_element(thetas.begin() + Q2 + 1, thetas.begin() + Q3, thetas.end());

        m_thetaIQR = thetas[Q3] - thetas[Q1];
    }

    // m_areaTotal
    {
        m_areaTotal = extract::sum(accArea);
    }

    // m_heightStd, m_heightCb
    {
        m_heightStd = sqrt(extract::variance(accHeight));
        m_heightCb = extract::sum(accHeightCb);
    }

    const time_point time_end = std::chrono::high_resolution_clock::now();
    const Duration time_duration = Duration(time_start, time_end);
    LOG_MESSAGE("5 features computed in " + time_duration.str());
    Console::out << "5 features computed in " << time_duration.str() << std::endl;
}

void StatisticsTool::compute25Statistics()
{
    const time_point time_start = std::chrono::high_resolution_clock::now();

    // Compute for: m_phiAnisotropy,
    // m_thetaMean, m_thetaStd, m_thetaMax, m_thetaCV, m_thetaMAD, m_thetaMedian, m_thetaQ3, m_thetaIQR, m_thetaQCD, m_thetaSkewness, m_thetaKurtosis,
    // m_areaTotal, m_areaStd,
    // m_heightMean, m_heightStd, m_heightMax, m_heightCV, m_heightMedian, m_heightQ3, m_heightIQR, m_heightMAD, m_heightCb,
    // m_corrAT, m_corrTH
    m_phiAnisotropy = computeAnisotropy();
    accumulator_set<ld, stats<tag::max, tag::immediate_weighted_mean, tag::weighted_variance, tag::weighted_skewness, tag::weighted_kurtosis>, ld > accTheta;
    std::vector<ld> thetas;
    accumulator_set<ld, stats<tag::variance, tag::sum > > accArea;
    std::vector<ld> areas;
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_mean, tag::variance> > accHeight;
    std::vector<ld> heights;
    // Set up the accumulator for C_b
    accumulator_set<ld, stats<tag::sum > > accHeightCb;
    std::vector<scal> centroids(3);
    std::vector<std::set<int>> sets = mesh->heightSeparation((int)centroids.size(), centroids); // each set contains the face's ids
    int min_idx = std::min_element(centroids.begin(), centroids.end()) - centroids.begin();
    std::set<int>& low_cluster = sets[min_idx];
    accumulator_set<ld, stats<tag::covariance<ld, tag::covariate1> > > accCovAT;
    accumulator_set<ld, stats<tag::covariance<ld, tag::covariate1> > > accCovTH;

    // macrosurface area
    const scal border = 1. - Parameters::get()->currentParams()->sideEffectParams.borderPercentage;
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

        accTheta(theta, weight = faceA);
        thetas.push_back(theta);
        accArea(faceA);
        areas.push_back(faceA);
        accHeight(height);
        heights.push_back(height);
        if (low_cluster.contains(faceID)) {
            const ld cosTheta = static_cast<ld>(dot(faceN, mesh->meso_normal));
            accHeightCb(faceA * cosTheta);
        }
        accCovAT(faceA, covariate1 = theta);
        accCovTH(theta, covariate1 = height);
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
        accCovTH(n.theta(), covariate1 = p.z);
    }


    // m_thetaMean, m_thetaStd, m_thetaMax, m_thetaCV, m_thetaMAD, m_thetaMedian, m_thetaQ3, m_thetaIQR, m_thetaQCD, m_thetaSkewness, m_thetaKurtosis
    {
        m_thetaMean = ld(extract::weighted_mean(accTheta));
        m_thetaStd = sqrt(extract::weighted_variance(accTheta));
        m_thetaRange.upper = ld(extract::max(accTheta));

        m_thetaCV = m_thetaStd / m_thetaMean;
        m_thetaMAD = median_absolute_deviation(thetas.begin(), thetas.end());

        auto const Q1 = thetas.size() / 4;
        auto const Q2 = thetas.size() / 2;
        auto const Q3 = Q1 + Q2;

        std::nth_element(thetas.begin(), thetas.begin() + Q1, thetas.end());
        std::nth_element(thetas.begin() + Q1 + 1, thetas.begin() + Q2, thetas.end());
        std::nth_element(thetas.begin() + Q2 + 1, thetas.begin() + Q3, thetas.end());

        m_thetaMedian = thetas[Q2];
        m_thetaQ3 = thetas[Q3];
        m_thetaIQR = m_thetaQ3 - thetas[Q1];
        m_thetaQCD = m_thetaIQR / (m_thetaQ3 + thetas[Q1]);

        // shape parameters
        m_thetaSkewness = extract::weighted_skewness(accTheta);
        m_thetaKurtosis = extract::weighted_kurtosis(accTheta);
    }


    // m_areaTotal, m_areaStd,
    {
        m_areaTotal = extract::sum(accArea);
        m_areaStd = sqrt(extract::variance(accArea));
    }

    // m_heightMean, m_heightStd, m_heightMax, m_heightCV, m_heightMedian, m_heightQ3, m_heightIQR, m_heightMAD, m_heightCb,
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

        m_heightMedian = heights[Q2] - m_heightRange.lower;
        m_heightQ3 = heights[Q3] - m_heightRange.lower;
        m_heightIQR = m_heightQ3 - heights[Q1] + m_heightRange.lower;

        m_heightRange.lower = 0;

        m_heightCb = extract::sum(accHeightCb); // faceA are already relatives to A, just as if the patch was 1 unit
    }

    // m_corrAT, m_corrTH
    {
        ld covTH = extract::covariance(accCovTH);
        m_THCorr = covTH / (m_thetaStd * m_heightStd);

        ld covAT = extract::covariance(accCovAT);
        m_ATCorr = covAT / (m_areaStd * m_thetaStd);
    }

    const time_point time_end = std::chrono::high_resolution_clock::now();
    const Duration time_duration = Duration(time_start, time_end);
    LOG_MESSAGE("25 features computed in " + time_duration.str());
    Console::out << "25 features computed in " << time_duration.str() << std::endl;
}


void StatisticsTool::computeStatistics()
{
#ifdef STATS_PHI
    m_phiAnisotropy = computeAnisotropy();
#endif // STATS_PHI
#ifdef STATS_THETA
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_weighted_mean, tag::weighted_variance, tag::weighted_skewness, tag::weighted_kurtosis>, ld > accTheta;
    std::vector<ld> thetas;
#endif // STATS_THETA
#ifdef STATS_AREA
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_mean, tag::variance, tag::sum > > accArea;
    std::vector<ld> areas;
#endif // STATS_AREA
#ifdef STATS_HEIGHT
    accumulator_set<ld, stats<tag::min, tag::max, tag::immediate_mean, tag::variance> > accHeight;
    std::vector<ld> heights;
    // Set up the accumulator for C_b
    accumulator_set<ld, stats<tag::sum > > accHeightCb;
    std::vector<scal> centroids(3);
    std::vector<std::set<int>> sets = mesh->heightSeparation((int)centroids.size(), centroids); // each set contains the face's ids
    int min_idx = std::min_element(centroids.begin(), centroids.end()) - centroids.begin();
    std::set<int>& low_cluster = sets[min_idx];
#endif // STATS_HEIGHT
#ifdef STATS_CORR
    accumulator_set<ld, stats<tag::covariance<ld, tag::covariate1> > > accCovAT;
    accumulator_set<ld, stats<tag::covariance<ld, tag::covariate1> > > accCovTH;
    accumulator_set<ld, stats<tag::covariance<ld, tag::covariate1> > > accCovAH;
#endif // STATS_CORR

    // macrosurface area
    const scal border = 1. - Parameters::get()->currentParams()->sideEffectParams.borderPercentage;
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

        const ld height   = static_cast<ld>(barycenter.z);
        const ld theta    = static_cast<ld>(faceN.theta());

#ifdef STATS_THETA
        accTheta(theta, weight = faceA);
        thetas.push_back(theta);
#endif // STATS_THETA
#ifdef STATS_AREA
        accArea(faceA);
        areas.push_back(faceA);
#endif // STATS_AREA
#ifdef STATS_HEIGHT
        accHeight(height);
        heights.push_back(height);
        if (low_cluster.contains(faceID)) {
            const ld cosTheta = static_cast<ld>(dot(faceN, mesh->meso_normal));
            accHeightCb(faceA * cosTheta);
        }
#endif // STATS_HEIGHT
#ifdef STATS_CORR
        accCovAT(faceA, covariate1 = theta);
        accCovTH(theta, covariate1 = height);
        accCovAH(faceA, covariate1 = height);
#endif // STATS_CORR
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
        accCovTH(n.theta(), covariate1 = p.z);
#endif // STATS_CORR
    }
#endif // STATS_HEIGHT || STATS_CORR


#ifdef STATS_THETA
    {
        m_thetaMean = ld(extract::weighted_mean(accTheta));
        m_thetaStd = sqrt(extract::weighted_variance(accTheta));
        m_thetaRange.lower = ld(extract::min(accTheta));
        m_thetaRange.upper = ld(extract::max(accTheta));

        m_thetaCV = m_thetaStd / m_thetaMean;
        m_thetaMAD = median_absolute_deviation(thetas.begin(), thetas.end());

        auto const Q1 = thetas.size() / 4;
        auto const Q2 = thetas.size() / 2;
        auto const Q3 = Q1 + Q2;

        std::nth_element(thetas.begin(), thetas.begin() + Q1, thetas.end());
        std::nth_element(thetas.begin() + Q1 + 1, thetas.begin() + Q2, thetas.end());
        std::nth_element(thetas.begin() + Q2 + 1, thetas.begin() + Q3, thetas.end());

        m_thetaQ1     = thetas[Q1];
        m_thetaMedian = thetas[Q2];
        m_thetaQ3     = thetas[Q3];
        m_thetaIQR = m_thetaQ3 - m_thetaQ1;
        m_thetaQCD = m_thetaIQR / (m_thetaQ3 + m_thetaQ1);

        // shape parameters
        m_thetaSkewness = extract::weighted_skewness(accTheta);
        m_thetaKurtosis = extract::weighted_kurtosis(accTheta);
    }
#endif //  STATS_THETA


#ifdef STATS_AREA
    {
        m_areaTotal = extract::sum(accArea);
        m_areaMean = extract::mean(accArea);
        m_areaStd = sqrt(extract::variance(accArea));
        m_areaRange.lower = extract::min(accArea);
        m_areaRange.upper = extract::max(accArea);

        m_areaCV = m_areaStd / m_areaMean;
        m_areaMAD = median_absolute_deviation(areas.begin(), areas.end());

        auto const Q1 = areas.size() / 4;
        auto const Q2 = areas.size() / 2;
        auto const Q3 = Q1 + Q2;

        std::nth_element(areas.begin(), areas.begin() + Q1, areas.end());
        std::nth_element(areas.begin() + Q1 + 1, areas.begin() + Q2, areas.end());
        std::nth_element(areas.begin() + Q2 + 1, areas.begin() + Q3, areas.end());

        m_areaQ1 = areas[Q1];
        m_areaMedian = areas[Q2];
        m_areaQ3 = areas[Q3];
        m_areaIQR = m_areaQ3 - m_areaQ1;
        m_areaQCD = m_areaIQR / (m_areaQ3 + m_areaQ1);
    }
#endif // STATS_AREA

#ifdef STATS_HEIGHT
    {
        m_heightRange.lower = extract::min(accHeight);
        m_heightRange.upper = extract::max(accHeight) - m_heightRange.lower;
        m_heightMean = extract::mean(accHeight) - m_heightRange.lower;
        m_heightStd = sqrt(extract::variance(accHeight));

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

        m_heightCb = extract::sum(accHeightCb); // faceA are already relatives to A, just as if the patch was 1 unit
    }
#endif // STATS_HEIGHT

#ifdef STATS_CORR
    {
        ld covTH = extract::covariance(accCovTH);
        m_THCorr = covTH / (m_thetaStd * m_heightStd);

        ld covAT = extract::covariance(accCovAT);
        m_ATCorr = covAT / (m_areaStd * m_thetaStd);

        ld covAH = extract::covariance(accCovAH);
        m_AHCorr = covAH / (m_areaStd * m_heightStd);
    }
#endif // STATS_CORR

    LOG_STATS(*this);
}

// EXPORT
void StatisticsTool::CSVHeader(csv::CSVWriter* writer) {
    std::vector<csv::elem> headers;

    headers.push_back({ csv::elem::Tag::STRING, "error" });

    headers.push_back({ csv::elem::Tag::STRING, "anisotropy" });

    headers.push_back({ csv::elem::Tag::STRING, "meanT" });
    headers.push_back({ csv::elem::Tag::STRING, "stdT" });
    headers.push_back({ csv::elem::Tag::STRING, "minT" });
    headers.push_back({ csv::elem::Tag::STRING, "maxT" });

    headers.push_back({ csv::elem::Tag::STRING, "cvT" });
    headers.push_back({ csv::elem::Tag::STRING, "madT" });
    headers.push_back({ csv::elem::Tag::STRING, "q1T" });
    headers.push_back({ csv::elem::Tag::STRING, "q2T" });
    headers.push_back({ csv::elem::Tag::STRING, "q3T" });
    headers.push_back({ csv::elem::Tag::STRING, "iqrT" });
    headers.push_back({ csv::elem::Tag::STRING, "qcdT" });

    headers.push_back({ csv::elem::Tag::STRING, "skewness" });
    headers.push_back({ csv::elem::Tag::STRING, "kurtosis" });

    headers.push_back({ csv::elem::Tag::STRING, "totalA" });
    headers.push_back({ csv::elem::Tag::STRING, "meanA" });
    headers.push_back({ csv::elem::Tag::STRING, "stdA" });
    headers.push_back({ csv::elem::Tag::STRING, "minA" });
    headers.push_back({ csv::elem::Tag::STRING, "maxA" });

    headers.push_back({ csv::elem::Tag::STRING, "cvA" });
    headers.push_back({ csv::elem::Tag::STRING, "madA" });
    headers.push_back({ csv::elem::Tag::STRING, "q1A" });
    headers.push_back({ csv::elem::Tag::STRING, "q2A" });
    headers.push_back({ csv::elem::Tag::STRING, "q3A" });
    headers.push_back({ csv::elem::Tag::STRING, "iqrA" });
    headers.push_back({ csv::elem::Tag::STRING, "qcdA" });

    headers.push_back({ csv::elem::Tag::STRING, "meanH" });
    headers.push_back({ csv::elem::Tag::STRING, "stdH" });
    headers.push_back({ csv::elem::Tag::STRING, "minH" });
    headers.push_back({ csv::elem::Tag::STRING, "maxH" });

    headers.push_back({ csv::elem::Tag::STRING, "cvH" });
    headers.push_back({ csv::elem::Tag::STRING, "madH" });
    headers.push_back({ csv::elem::Tag::STRING, "q1H" });
    headers.push_back({ csv::elem::Tag::STRING, "q2H" });
    headers.push_back({ csv::elem::Tag::STRING, "q3H" });
    headers.push_back({ csv::elem::Tag::STRING, "iqrH" });
    headers.push_back({ csv::elem::Tag::STRING, "qcdH" });

    headers.push_back({ csv::elem::Tag::STRING, "CbH" });

    headers.push_back({ csv::elem::Tag::STRING, "corrAT" });
    headers.push_back({ csv::elem::Tag::STRING, "corrTH" });
    headers.push_back({ csv::elem::Tag::STRING, "corrAH" });

    writer->writeRow(headers);
}

void StatisticsTool::toCSV(csv::CSVWriter* writer) {
    std::vector<csv::elem> stats;

    stats.push_back({ DECIMAL_TAG, m_error }); // Error

    stats.push_back({ DECIMAL_TAG, m_phiAnisotropy });   // r

    stats.push_back({ DECIMAL_TAG, m_thetaMean });         // mean(T)
    stats.push_back({ DECIMAL_TAG, m_thetaStd });          // std(T)
    stats.push_back({ DECIMAL_TAG, m_thetaRange.lower });  // min(T)
    stats.push_back({ DECIMAL_TAG, m_thetaRange.upper });  // max(T)

    stats.push_back({ DECIMAL_TAG, m_thetaCV });           // cv(T)
    stats.push_back({ DECIMAL_TAG, m_thetaMAD });          // mad(T)
    stats.push_back({ DECIMAL_TAG, m_thetaQ1 });           // q1(T)
    stats.push_back({ DECIMAL_TAG, m_thetaMedian });       // q2(T)
    stats.push_back({ DECIMAL_TAG, m_thetaQ3 });           // q3(T)
    stats.push_back({ DECIMAL_TAG, m_thetaIQR });          // iqr(T)
    stats.push_back({ DECIMAL_TAG, m_thetaQCD });          // qcd(T)

    stats.push_back({ DECIMAL_TAG, m_thetaSkewness });     // skewness
    stats.push_back({ DECIMAL_TAG, m_thetaKurtosis });     // kurtosis

    stats.push_back({ DECIMAL_TAG, m_areaTotal });         // total(A)
    stats.push_back({ DECIMAL_TAG, m_areaMean });          // mean(A)
    stats.push_back({ DECIMAL_TAG, m_areaStd });           // std(A)
    stats.push_back({ DECIMAL_TAG, m_areaRange.lower });   // min(A)
    stats.push_back({ DECIMAL_TAG, m_areaRange.upper });   // max(A)

    stats.push_back({ DECIMAL_TAG, m_areaCV });            // cv(A)
    stats.push_back({ DECIMAL_TAG, m_areaMAD });           // mad(A)
    stats.push_back({ DECIMAL_TAG, m_areaQ1 });            // q1(A)
    stats.push_back({ DECIMAL_TAG, m_areaMedian });        // q2(A)
    stats.push_back({ DECIMAL_TAG, m_areaQ3 });            // q3(A)
    stats.push_back({ DECIMAL_TAG, m_areaIQR });           // iqr(A)
    stats.push_back({ DECIMAL_TAG, m_areaQCD });           // qcd(A)

    stats.push_back({ DECIMAL_TAG, m_heightMean });        // mean(H)
    stats.push_back({ DECIMAL_TAG, m_heightStd });         // std(H)
    stats.push_back({ DECIMAL_TAG, m_heightRange.lower }); // min(H)
    stats.push_back({ DECIMAL_TAG, m_heightRange.upper }); // max(H)

    stats.push_back({ DECIMAL_TAG, m_heightCV });          // cv(H)
    stats.push_back({ DECIMAL_TAG, m_heightMAD });         // mad(H)
    stats.push_back({ DECIMAL_TAG, m_heightQ1 });          // q1(H)
    stats.push_back({ DECIMAL_TAG, m_heightMedian });      // q2(H)
    stats.push_back({ DECIMAL_TAG, m_heightQ3 });          // q3(H)
    stats.push_back({ DECIMAL_TAG, m_heightIQR });         // iqr(H)
    stats.push_back({ DECIMAL_TAG, m_heightQCD });         // qcd(H)

    stats.push_back({ DECIMAL_TAG, m_heightCb });          // C_b

    stats.push_back({ DECIMAL_TAG, m_ATCorr });            // Corr(A ; T)
    stats.push_back({ DECIMAL_TAG, m_THCorr });            // Corr(T ; H)
    stats.push_back({ DECIMAL_TAG, m_AHCorr });            // Corr(A ; H)

    writer->writeRow(stats);
}


void StatisticsTool::print() {
    Console::out << *this << std::endl;
}


std::ostream& operator<<(std::ostream& output, const StatisticsTool& stats) {
    output << Console::timePad << "Statistics: "  << stats.mesh->name << std::endl;
    output << Console::timePad << "-- Error " << std::endl;
    output << Console::timePad << "   - SMAPE: " << stats.m_error << std::endl;
    output << Console::timePad << "-- Thetas " << std::endl;
    output << Console::timePad << "   - Mean    : " << stats.m_thetaMean     << std::endl;
    output << Console::timePad << "   - Std     : " << stats.m_thetaStd      << std::endl;
    output << Console::timePad << "   - Range   : " << stats.m_thetaRange    << std::endl;
    output << Console::timePad << "   - Cv      : " << stats.m_thetaCV       << std::endl;
    output << Console::timePad << "   - MAD     : " << stats.m_thetaMAD      << std::endl;
    output << Console::timePad << "   - Q1      : " << stats.m_thetaQ1       << std::endl;
    output << Console::timePad << "   - Q2      : " << stats.m_thetaMedian   << std::endl;
    output << Console::timePad << "   - Q3      : " << stats.m_thetaQ3       << std::endl;
    output << Console::timePad << "   - IQR     : " << stats.m_thetaIQR      << std::endl;
    output << Console::timePad << "   - QCD     : " << stats.m_thetaQCD      << std::endl;
    output << Console::timePad << "   - Skewness: " << stats.m_thetaSkewness << std::endl;
    output << Console::timePad << "   - Kurtosis: " << stats.m_thetaKurtosis << std::endl;
    output << Console::timePad << "-- Areas " << std::endl;
    output << Console::timePad << "   - Total   : " << stats.m_areaTotal  << std::endl;
    output << Console::timePad << "   - Mean    : " << stats.m_areaMean   << std::endl;
    output << Console::timePad << "   - Std     : " << stats.m_areaStd    << std::endl;
    output << Console::timePad << "   - Range   : " << stats.m_areaRange  << std::endl;
    output << Console::timePad << "   - Cv      : " << stats.m_areaCV     << std::endl;
    output << Console::timePad << "   - MAD     : " << stats.m_areaMAD    << std::endl;
    output << Console::timePad << "   - Q1      : " << stats.m_areaQ1     << std::endl;
    output << Console::timePad << "   - Q2      : " << stats.m_areaMedian << std::endl;
    output << Console::timePad << "   - Q3      : " << stats.m_areaQ3     << std::endl;
    output << Console::timePad << "   - IQR     : " << stats.m_areaIQR    << std::endl;
    output << Console::timePad << "   - QCD     : " << stats.m_areaQCD    << std::endl;
    output << Console::timePad << "-- Heights " << std::endl;
    output << Console::timePad << "   - Mean    : " << stats.m_heightMean   << std::endl;
    output << Console::timePad << "   - Std     : " << stats.m_heightStd    << std::endl;
    output << Console::timePad << "   - Range   : " << stats.m_heightRange  << std::endl;
    output << Console::timePad << "   - Cv      : " << stats.m_heightCV     << std::endl;
    output << Console::timePad << "   - MAD     : " << stats.m_heightMAD    << std::endl;
    output << Console::timePad << "   - Q1      : " << stats.m_heightQ1     << std::endl;
    output << Console::timePad << "   - Q2      : " << stats.m_heightMedian << std::endl;
    output << Console::timePad << "   - Q3      : " << stats.m_heightQ3     << std::endl;
    output << Console::timePad << "   - IQR     : " << stats.m_heightIQR    << std::endl;
    output << Console::timePad << "   - QCD     : " << stats.m_heightQCD    << std::endl;
    output << Console::timePad << "   - Cb      : " << stats.m_heightCb     << std::endl;
    output << Console::timePad << "-- Correlations " << std::endl;
    output << Console::timePad << "  -- AT      : " << stats.m_ATCorr << std::endl;
    output << Console::timePad << "  -- TH      : " << stats.m_THCorr << std::endl;
    output << Console::timePad << "  -- AH      : " << stats.m_AHCorr << std::endl;

    return output;
}