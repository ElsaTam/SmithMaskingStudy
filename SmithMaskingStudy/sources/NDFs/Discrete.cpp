/*
    MIT License

    Copyright (c) 2020 Arthur Cavalier, Mickaël Ribardière and Benjamin Bringier

    Permission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal in
    the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is furnished to do
    so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
    DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
    ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.s
*/
#include "NDFs/Microfacet.hpp"

#include <algorithm>
#include <unordered_set>
#include <random>
#include "utils/math/math.h"
#include "utils/console.h"
#include "utils/params.h"
#include "tools/logger.h"


float Discrete::getTheta(size_t i) const {
    if (i < 0) {
        Console::err << "Error: trying to access theta at index " << i << " (over " << thetas.size() << ")" << std::endl;
        return thetas[0]; // error
    }
    if (i >= thetas.size()) {
        Console::err << "Error: trying to access theta at index " << i << " (over " << thetas.size() << ")" << std::endl;
        return thetas[thetas.size() - 1]; // error
    }
    return thetas[i];
}

float Discrete::getPhi(size_t i) const {
    i = i % phis.size();
    return phis[i];
}

float Discrete::getNextTheta(size_t i) const {
    if (i >= thetas.size() - 1) return 2 * thetas[i] - thetas[i - 1]; // assume constant dTheta
    return thetas[i + 1];
}

float Discrete::getNextPhi(size_t i) const {
    i = i % phis.size();
    if (i == phis.size() - 1) return phis[0] + m_2_pi;
    return phis[i + 1];
}

float Discrete::getCenterTheta(size_t i) const {
    if (i >= thetas.size() - 1) return 0;
    return 0.5f * (getNextTheta(i) + getTheta(i));
}

float Discrete::getCenterPhi(size_t i) const {
    if (i >= phis.size()) {
        Console::err << "Error: trying to access phi at index " << i << " (over " << phis.size() << ")" << std::endl;
        return 0;
    }
    return 0.5f * (getNextPhi(i) + getPhi(i));
}

void Discrete::findIndex(const vec3sc& w, size_t& thetaIndex, size_t& phiIndex) const
{
    vec3sc w_spherical = Conversion::cartesian_to_polar(w); // Spherical coordinates as (theta,phi,radius)
    scal theta = w_spherical[0]; // n.theta() in [0, pi/2]
    scal phi = w_spherical[1] >= phiEnd() ? phiStart() + w_spherical[1] - phiEnd() : w_spherical[1]; // n.phi() in [-pi, pi]

    scal phiStep = (phiEnd() - phiStart()) / (scal)phiSize();
    scal thetaStep = (thetaEnd() - thetaStart()) / (scal)thetaSize();
    phiIndex = (int)floor((phi - phiStart()) / phiStep) % phiSize();
    thetaIndex = floor((theta - thetaStart()) / thetaStep);

    /*float theta = acos(w.z); // acos in range [0,pi]
    float phi = atan2(w.y, w.x); // atan2 in range [-pi,+pi]
    phiIndex = 0, thetaIndex = 0;
    while (thetaIndex < thetas.size() && getTheta(thetaIndex) <= theta) ++thetaIndex;
    --thetaIndex;
    while (phiIndex < phis.size() && getPhi(phiIndex) <= phi)   ++phiIndex;
    --phiIndex;*/
}

scal Discrete::findInterpolatedValue(const std::vector<std::vector<float>>& V, const vec3sc& w) const
{
    // x ---> phi   (-pi vers pi  )  / colonnes (gauche vers droite)
    // y ---> theta ( 0  vers pi/2)  / lignes   (haut vers bas)
    float theta = acos(w.z); // acos in range [0,pi]
    float phi = atan2(w.y, w.x); // atan2 in range [-pi,+pi]
    size_t tIdx, pIdx;
    findIndex(w, tIdx, pIdx);

    size_t cLR = pIdx; // Center Left-Right
    size_t L = (cLR - (size_t)1) % phis.size(); // left,  previous phi
    size_t R = (cLR + (size_t)1) % phis.size(); // right, next phi
    size_t cTB = tIdx; // Center Top-Bottom
    size_t T, B;
    if (cTB == thetas.size() - 1) { // w almost lying on the plane, interpolate only on phi
        T = B = cTB;
    }
    else if (cTB == 0) { // we are almost at the normal, we need to swap some indices
        // TODO
        T = B = cTB;
    }
    else {
        T = cTB - 1;
        B = cTB + 1;
    }

    // find the quarter in which w lies
    vec2sc C = { getCenterPhi(cLR), getCenterTheta(cTB) }; // theta/phi coordinates for the value at [pIdx][tIdx]
    /*
             L        cLR        R
         --------- --------- ---------
        |         |         |         |
    T   |         |         |         | T
        |         |         |         |
         --------- --------- ---------
        |         |    .    |         |
    cTB |         |. . C . .|         | cTB
        |         |    .    |         |
         --------- --------- ---------
        |         |         |         |
    B   |         |         |         | B
        |         |         |         |
         --------- --------- ---------
             L        cLR        R
    */

    // indices of the theta and phi surrounding C
    size_t pL = C.x > phi   ? L : cLR;
    size_t pR = C.x > phi   ? cLR : R;
    size_t tT = C.y > theta ? T : cTB;
    size_t tB = C.y > theta ? cTB : B;

    // phi1 < phi < phi2  &&  theta1 < theta < theta2
    float phi1 = getCenterPhi(pL);
    float phi2 = getCenterPhi(pR);
    float theta1 = getCenterTheta(tB);
    float theta2 = getCenterTheta(tT);

    float BI = Maths::bilinearInterpolation(
        V[pL][tB], V[pL][tT], V[pR][tB], V[pR][tT],
        phi1, phi2,
        theta1, theta2,
        phi, theta);
    //float BI_0 = V[pIdx][tIdx]; // for debug purpose
    return BI;
}

scal Discrete::findValue(const std::vector<std::vector<float>>& V, const vec3sc& w) const
{
    size_t tIdx, pIdx;
    findIndex(w, tIdx, pIdx);
    return V[pIdx][tIdx];
}

void Discrete::addNormal(const vec3sc& n, const scal weight)
{
    // Find the correct cell
    size_t phiIdx = -1, thetaIdx = -1;
    findIndex(n, thetaIdx, phiIdx);
    // Add the normal to the distribution histogram
    //if (thetaIdx >= thetaSize()) {
    //    Console::err << n << std::endl;
    //    Console::err << "Phi index : " << phiIdx << std::endl;
    //    Console::err << "Theta index : " << thetaIdx << std::endl;
    //}
    D_values[phiIdx][thetaIdx] += weight;
}

Discrete::Discrete(const TriangleMesh& mesh, const MicrofacetDistribution* _matchingNDF, scal borderPercentage)
    : matchingNDF(_matchingNDF)
{
    Console::out << Console::timeStamp << "Building Discrete NDF from mesh..." << std::endl;
    m_type = NDFTypes::DISCRETE;
    m_meso_normal = mesh.meso_normal;

    scal phiStep;
    scal thetaStep;

    // -------- Step 0 : initialisation
    //          * Angle ranges
    //          * Tabulation size
    //          * Distribution vectors
    {
        phiStep = (phiEnd() - phiStart()) / (scal)phiSize();
        thetaStep = (thetaEnd() - thetaStart()) / (scal)thetaSize();

        D_values.resize(phiSize());
        pdf_values.resize(phiSize());
        for (int p = 0; p < phiSize(); ++p) {
            const scal phi = phiStart() + p * phiStep;
            phis.push_back(phi);
            D_values[p].resize(thetaSize(), 0);
            pdf_values[p].resize(thetaSize(), 0);
        }
        for (int t = 0; t < thetaSize(); ++t) {
            const scal theta = thetaStart() + t * thetaStep;
            thetas.push_back(theta);
        }
        m_average_normal = vec3sc(0, 0, 0);
        m_D_integral = 0;
    }

    // -------- Step 1 :
    //          * Collect normals in the hemisphere cells
    {
        for (int faceID = 0; faceID < mesh.index.size(); ++faceID)
        {
            // If we are too close to the edge of the surface, discard the point for the minimal visibility
            // (if optixLaunchParams.sideEffect.borderPercentage > 0)
            const vec3i& vIdx = mesh.index[faceID];
            const vec3sc& A = mesh.vertex[vIdx[0]];
            const vec3sc& B = mesh.vertex[vIdx[1]];
            const vec3sc& C = mesh.vertex[vIdx[2]];
            const vec3sc surfPos = (A + B + C)/(scal)3.;
            if (mesh.bounds.closest_distance(surfPos).x < borderPercentage * mesh.bounds.span().x / 2.f
                || mesh.bounds.closest_distance(surfPos).y < borderPercentage * mesh.bounds.span().y / 2.f)
            {
                continue;
            }

            const scal a = mesh.area[faceID];
            const vec3sc& n = mesh.triangle_normal[faceID];
            addNormal(n, a);
            m_average_normal += n * a;
        }
    }

    // -------- Step 2 :
    //          * Weighting by cell sizes to obtain D
    //          * Find the normalization factor: 1 / \int D(w)*cos(theta)
    scal normalizationFactor = 0;
    {
        scal dTheta = (thetaEnd() - thetaStart()) / thetaSize();
        scal dPhi = (phiEnd() - phiStart()) / phiSize();
        for (size_t p = 0; p < phiSize(); ++p) {
            const scal phi = phiStart() + p * phiStep;
            for (size_t t = 0; t < thetaSize(); ++t) {
                // PSurf[i][j] must be weighted by the patch surface
                // sphereSurface  = \int_{0}^{m_pi} sin\theta d\theta * \int_{0}^{m_2_pi} d\phi
                //                = [-cos(m_pi) + cos(0)] * [m_2_pi - 0]
                //                = [1 + 1] * m_2_pi
                //                = m_4_pi
                // hemisphereSurf = sphereSurface / 2
                //                = m_2_pi
                // patchSurface   = \int_{\theta_{min}}^{\theta_{max}} sin\theta d\theta * \int_{\phi_{min}}^{\phi_{max}} d\phi
                //                = [-cos(\theta_{max}) + cos(\theta_{min})] * [\phi_{max} - \phi_{min}]
                // patchFraction  = patchSurface / hemisphereSurf
                //                = patchSurface / m_2_pi
                scal patchThetaMin = thetaStart() + t * thetaStep;
                scal patchThetaMax = thetaStart() + (t + 1) * thetaStep;
                scal patchSurface = (-cos(patchThetaMax) + cos(patchThetaMin)) * phiStep;
                scal patchFraction = patchSurface / m_2_pi;

                scal coneVolume = sin(getCenterTheta(t)) * dTheta * dPhi;

                D_values[p][t] = D_values[p][t] / patchSurface;
                pdf_values[p][t] = D_values[p][t] * cos(getCenterTheta(t));

                // Increase the normalization factor
                // N.B.: at this point, D[i][j] =  D
                //       but since we are working on a surface element, we need to add a factor
                normalizationFactor += D_values[p][t] * cos(getCenterTheta(t)) * sin(getCenterTheta(t)); // (double) Riemann sum
            }
        }
        normalizationFactor *= dTheta * dPhi;
    }

    // -------- Step 3 :
    //          * Normalize
    //          * Write in the Tabulation structure
    {
        for (int p = 0; p < phiSize(); ++p) {
            std::vector<float> pdf_phi;
            for (int t = 0; t < thetaSize(); ++t) {
                D_values[p][t] = D_values[p][t] / normalizationFactor;
                pdf_values[p][t] = pdf_values[p][t] / normalizationFactor;
            }
        }
        m_average_normal /= mesh.surfaceArea;
        m_average_normal = normalize(m_average_normal);
        Console::light << Console::timePad << "Average micronormals = " << m_average_normal << std::endl;
    }

    // -------- Step 4 :
    //          * Build the CDF (for sampling)
    {
        cdf.reserve(phiSize() * thetaSize() + 1); // one more value for the first entry = 0
        cdf.push_back({ -1, -1, 0 });
        for (int p = 0; p < phiSize(); ++p) {
            for (int t = 0; t < thetaSize(); ++t) {
                //scal patchThetaMin = thetaStart() + t * thetaStep;
                //scal patchThetaMax = thetaStart() + (t + 1) * thetaStep;
                //scal patchSurface = (-cos(patchThetaMax) + cos(patchThetaMin)) * phiStep;
                //scal patchFraction = patchSurface / m_2_pi;
                if (pdf_values[p][t] > 0)
                    cdf.push_back({ t, p, std::get<2>(cdf[cdf.size() - 1]) + pdf_values[p][t] });
                m_D_integral += D_values[p][t] * sin(getCenterTheta(t));
            }
        }
        scal dTheta = (thetaEnd() - thetaStart()) / thetaSize();
        scal dPhi = (phiEnd() - phiStart()) / phiSize();
        m_D_integral *= dTheta * dPhi;

        cdf_sum = std::get<2>(cdf[cdf.size() - 1]);
        Console::light << Console::timePad << "Cdf sum = " << cdf_sum << std::endl;
        if (cdf_sum != 1) {
            for (size_t i = 1; i < cdf.size(); ++i)
                std::get<2>(cdf[i]) /= cdf_sum;
            std::get<2>(cdf[cdf.size() - 1]) = 1.0f;
        }
    }

    Console::out << *this << std::endl;

    /*for (int p = 0; p < phiSize(); ++p) {
       Console::out << Console::timePad
            << "phi("       << p << ") = " << getPhi(p) << " / "
            << "nextPhi("   << p << ") = " << getNextPhi(p) << " / "
            << "centerPhi(" << p << ") = " << getCenterPhi(p) << " / " << std::endl;
    }*/
}

Discrete::Discrete(const MicrofacetDistribution* analyticNDF, int samples)
    : matchingNDF(analyticNDF)
{
    Console::out << Console::timeStamp << "Building Discrete NDF from distribution..." << std::endl;
    m_type = NDFTypes::DISCRETE;

    size_t nPhi, nTheta;
    scal phiStep;
    scal thetaStep;

    {
        nPhi = phiSize(); nTheta = thetaSize();
        phiStep = (phiEnd() - phiStart()) / (scal)nPhi;
        thetaStep = (thetaEnd() - thetaStart()) / (scal)nTheta;

        D_values.resize(nPhi);
        pdf_values.resize(nPhi);
        for (int p = 0; p < nPhi; ++p) {
            const scal phi = phiStart() + p * phiStep;
            phis.push_back(phi);
            D_values[p].resize(nTheta, 0);
            pdf_values[p].resize(nTheta, 0);
        }
        for (int t = 0; t < nTheta; ++t) {
            const scal theta = thetaStart() + t * thetaStep;
            thetas.push_back(theta);
        }
    }

    {
        std::unordered_set<int> addedNormals;

        std::default_random_engine rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<scal> rng(0.f, 1.f);
        for (int i = 0; i < samples; ++i)
        {
            scal s1 = rng(gen);
            scal s2 = rng(gen);
            const vec3sc& n = analyticNDF->sample(s1, s2);
            addNormal(n, 1. / n.z);
        }
    }

    scal normalizationFactor = 0;
    {
        scal dTheta = (thetaEnd() - thetaStart()) / thetaSize();
        scal dPhi = (phiEnd() - phiStart()) / phiSize();
        for (size_t p = 0; p < nPhi; ++p) {
            const scal phi = phiStart() + p * phiStep;
            for (size_t t = 0; t < nTheta; ++t) {
                scal patchThetaMin = thetaStart() + t * thetaStep;
                scal patchThetaMax = thetaStart() + (t + 1) * thetaStep;
                scal patchSurface = (-cos(patchThetaMax) + cos(patchThetaMin)) * phiStep;
                scal patchFraction = patchSurface / m_2_pi;

                scal coneVolume = sin(getCenterTheta(t)) * dTheta * dPhi;

                D_values[p][t] = D_values[p][t] / patchSurface;
                pdf_values[p][t] = D_values[p][t] * cos(getCenterTheta(t));
                normalizationFactor += D_values[p][t] * cos(getCenterTheta(t)) * coneVolume;
            }
        }
    }

    {
        for (int p = 0; p < nPhi; ++p) {
            std::vector<float> pdf_phi;
            for (int t = 0; t < nTheta; ++t) {
                D_values[p][t] = D_values[p][t] / normalizationFactor;
                pdf_values[p][t] = pdf_values[p][t] / normalizationFactor;
            }
        }
    }

    {
        cdf.reserve(nPhi * nTheta + 1); // one more value for the first entry = 0
        cdf.push_back({ -1, -1, 0 });
        for (int p = 0; p < nPhi; ++p) {
            for (int t = 0; t < nTheta; ++t) {
                scal patchThetaMin = thetaStart() + t * thetaStep;
                scal patchThetaMax = thetaStart() + (t + 1) * thetaStep;
                scal patchSurface = (-cos(patchThetaMax) + cos(patchThetaMin)) * phiStep;
                scal patchFraction = patchSurface / m_2_pi;
                if (pdf_values[p][t] > 0)
                    cdf.push_back({ t, p, std::get<2>(cdf[cdf.size() - 1]) + pdf_values[p][t] });
            }
        }
        cdf_sum = std::get<2>(cdf[cdf.size() - 1]);
        Console::light << Console::timePad << "Cdf sum = " << cdf_sum << std::endl;
        if (cdf_sum != 1) {
            for (size_t i = 1; i < cdf.size(); ++i)
                std::get<2>(cdf[i]) /= cdf_sum;
            std::get<2>(cdf[cdf.size() - 1]) = 1.0f;
        }
    }
}

/**
* @brief Normal distribution evaluation function
*
* This function must be statisfy :
* \f$ \int_{\Omega_+}D(\omega_h) cos(\theta_h) d\omega_h = 1 \f$
*
* pdf(wh) = dot(wh, wg) * D(wh)
* D(wh) = pdf(wh) / dot(wh, wg)
*
* @param wh        The microfacet normal (half-vector)
* @return scal     \f$  D(\omega_h) \f$
*/
scal Discrete::D(const vec3sc& wh) const
{
    return findInterpolatedValue(D_values, wh);
}

/**
* @brief Compute the pdf for a given direction.
*
* @param u1        the first uniformly distributed sample used to compute w
* @param u2        the second uniformly distributed sample used to compute w
* @param w         the given direction
* @return scal     pdf(w)
*/
scal Discrete::pdf(const vec3sc& w, float u1, float u2) const
{
    return findInterpolatedValue(pdf_values, w);
}

scal Discrete::gk_pos(const vec3sc& wo) const
{
    scal integral = 0;
    scal dTheta = (thetaEnd() - thetaStart()) / thetaSize();
    scal dPhi = (phiEnd() - phiStart()) / phiSize();
    for (int p = 0; p < phiSize(); ++p) {
        for (int t = 0; t < thetaSize(); ++t) {
            vec3sc wm = Conversion::polar_to_cartesian(getCenterTheta(t), getCenterPhi(p));
            scal pm = D_values[p][t] / m_D_integral;
            integral += std::max((scal)0, dot(wo, wm)) * pm * sin(getCenterTheta(t));
        }
    }
    integral *= dTheta * dPhi;
    return integral;
}

scal Discrete::gk_neg(const vec3sc& wo) const
{
    scal integral = 0;
    scal dTheta = (thetaEnd() - thetaStart()) / thetaSize();
    scal dPhi = (phiEnd() - phiStart()) / phiSize();
    for (int p = 0; p < phiSize(); ++p) {
        for (int t = 0; t < thetaSize(); ++t) {
            vec3sc wm = Conversion::polar_to_cartesian(getCenterTheta(t), getCenterPhi(p));
            scal pm = D_values[p][t] / m_D_integral;
            integral += std::min((scal)0, dot(wo, wm)) * pm * sin(getCenterTheta(t));
        }
    }
    integral *= dTheta * dPhi;
    return integral;
}

/**
* G_1(w) = \frac{\cos\theta g(n)}{g(k)}
*/
scal Discrete::G1_ashikhmin(const vec3sc& w) const
{
    return dot(w, m_meso_normal) * gk_pos(m_meso_normal) / gk_pos(w);
}

/**
* G_1(w) = \frac{\cos\theta}{\int_\Omega (w \cdot w_h) D(w_h) dw_h}
*/
scal Discrete::G1_heitz(const vec3sc& w) const
{
    scal integral = 0;
    scal dTheta = (thetaEnd() - thetaStart()) / thetaSize();
    scal dPhi = (phiEnd() - phiStart()) / phiSize();
    for (int p = 0; p < phiSize(); ++p) {
        for (int t = 0; t < thetaSize(); ++t) {
            vec3sc wm = Conversion::polar_to_cartesian(getCenterTheta(t), getCenterPhi(p));
            integral += std::max((scal)0, dot(w, wm)) * D_values[p][t] * sin(getCenterTheta(t));
        }
    }
    integral *= dTheta * dPhi;
    return dot(w, m_meso_normal) / integral;
}

/**
* @brief Smith's shadowing-masking term for a single direction
*
* @param w         a direction vector
* @param wh        the microfacet normal
* @return scal    \f$ G(\omega,\omega_h) \f$
*/
scal Discrete::G1(const vec3sc& wo, const vec3sc& wh) const
{
    if (dot(wo, wh) <= 0) return 0;
    return G1_ashikhmin(wo);
}

/**
* @brief Correlated Smith's shadowing-masking term
*
* @param wi        the incident direction
* @param wo        the outgoing direction
* @param wh        the microfacet normal
* @return scal    \f$ G(\omega_i, \omega_o, \omega_h) \f$
*/
scal Discrete::GAFcorrelated(const vec3sc& wi, const vec3sc& wo, const vec3sc& wh) const
{
    if (matchingNDF)
        return matchingNDF->GAFcorrelated(wi, wo, wh);
    else {
        Console::err << "[Error] Discrete::GAFcorrelated not implemented." << std::endl;
        exit(EXIT_FAILURE);
    }
}


/**
* @brief Sample the microfacet normal distribution
*
* @param u1        an uniformly distributed sample
* @param u2        an uniformly distributed sample
* @return vec3sc    the sampled normal
*/
vec3sc Discrete::sample(float u1, float u2) const
{
    std::vector<std::tuple<int, int, float>>::const_iterator entry = std::lower_bound(
        cdf.begin(),
        cdf.end(),
        u1,
        [](const std::tuple<int, int, float>& info, double value) {
            return std::get<2>(info) < value;
        }
    );
    size_t index = std::min(cdf.size() - 2, (size_t)std::max((ptrdiff_t)0, entry - cdf.begin() - 1));
    // Handle a rare corner-case where a entry has probability 0 but is sampled nonetheless
    while (index < cdf.size() - 1 && std::get<2>(cdf[index + 1]) - std::get<2>(cdf[index]) == 0) {
        ++index;
    }
    int tIdx = std::get<0>(cdf[index]);
    int pIdx = std::get<1>(cdf[index]);
    return Conversion::polar_to_cartesian(getCenterTheta(tIdx), getCenterPhi(pIdx));
}

/**
* @brief Sample the microfacet normal distribution
*
* @param u1        an uniformly distributed sample
* @param u2        an uniformly distributed sample
* @param _pdf       out parameter, for the pdf of the sampled normal
* @param _D         out parameter, for the distribution of the sampled normal
* @return vec3sc    the sampled normal
*/
vec3sc Discrete::sample(float u1, float u2, scal& _pdf, scal& _D) const
{
    std::vector<std::tuple<int, int, float>>::const_iterator entry = std::lower_bound(
        cdf.begin(),
        cdf.end(),
        u1,
        [](const std::tuple<int, int, float>& info, double value) {
            return std::get<2>(info) < value;
        }
    );
    size_t index = std::min(cdf.size() - 2, (size_t)std::max((ptrdiff_t)0, entry - cdf.begin() - 1));
    _pdf = std::get<2>(cdf[index + 1]) - std::get<2>(cdf[index]);
    // Handle a rare corner-case where a entry has probability 0 but is sampled nonetheless
    while (index < cdf.size() - 1 && _pdf == 0) {
        ++index;
        _pdf = std::get<2>(cdf[index + 1]) - std::get<2>(cdf[index]);
    }
    int tIdx = std::get<0>(cdf[index]);
    int pIdx = std::get<1>(cdf[index]);
    vec3sc w = Conversion::polar_to_cartesian(getCenterTheta(tIdx), getCenterPhi(pIdx));

    //scal patchThetaMin = thetaStart() + tIdx * ((thetaEnd() - thetaStart()) / (scal)thetas.size());
    //scal patchThetaMax = thetaStart() + (tIdx + 1) * ((thetaEnd() - thetaStart()) / (scal)thetas.size());
    //scal patchSurface = (-cos(patchThetaMax) + cos(patchThetaMin)) * ((phiEnd() - phiStart()) / (scal)phis.size());
    //scal patchFraction = patchSurface / m_2_pi;
    _pdf *= cdf_sum;
    //_D = _pdf / w.z / patchFraction;
    _D = D(w);
    //_pdf = pdf(w);
    return w;
}

void Discrete::toCSV(const std::string& output) const {
    csv::CSVWriter* writer_distrib = new csv::CSVWriter(output);

    scal thetaStep = (thetaEnd() - thetaStart()) / (scal)thetaSize();
    scal phiStep   = (phiEnd()   - phiStart())   / (scal)phiSize();
    
    std::vector<csv::elem> thetas;
    thetas.push_back({ csv::elem::Tag::STRING, "" });
    for (int j = 0; j < thetaSize(); ++j) { thetas.push_back({ csv::elem::Tag::SCAL, thetaStart() + j * thetaStep }); }
    writer_distrib->writeRow(thetas);

    for (int i = 0; i < phiSize(); ++i) {
        std::vector<csv::elem> vector_d;
        vector_d.push_back({ csv::elem::Tag::SCAL, phiStart() + i * phiStep });

        scal phi = phiStart() + (i + 0.5) * phiStep;

        for (int j = 0; j < thetaSize(); ++j) {
            scal theta = thetaStart() + (j + 0.5) * thetaStep;
            vector_d.push_back({ csv::elem::Tag::SCAL, D(Conversion::polar_to_cartesian(theta, phi)) });
        }

        writer_distrib->writeRow(vector_d);
    }

    writer_distrib->close();
    delete writer_distrib;
}

std::ostream& operator<<(std::ostream& output, const Discrete& ndf) {
    output << Console::timePad << "Discrete distribution: shape (nPhi, nTheta) = (" << ndf.phiSize() << ", " << ndf.thetaSize() << ")" << std::endl;
    output << Console::timePad << Console::indent << "Phi:   from " << ndf.phiStart() << " to " << ndf.phiEnd() << std::endl;
    output << Console::timePad << Console::indent << "Theta: from " << ndf.thetaStart() << " to " << ndf.thetaEnd() << std::endl;
    output << std::endl;

    // Print distribution
    int precision = 9;

    if (ndf.phiSize() <= 100 && ndf.thetaSize() <= 10)
    {
        // Print thetas (first row)
        output << Console::timePad << std::string(precision, ' ') << " | ";
        for (int t = 0; t < ndf.thetaSize(); ++t) {
            output << " ";
            output << Console::fixed_size(ndf.getTheta(t), precision);
            output << " ";
        }
        output << std::endl;

        // Header separation
        output << Console::timePad;
        for (int i = 0; i < ndf.thetaSize() + 1; ++i) {
            output << std::string(precision + 2, '-');
        }
        output << std::endl;

        // Rows
        for (int p = 0; p < ndf.phiSize(); ++p) {
            // theta
            output << Console::timePad << Console::fixed_size(ndf.getPhi(p), precision) << " | ";
            // values
            for (int t = 0; t < ndf.thetaSize(); ++t) {
                output << " ";
                output << Console::fixed_size(ndf[p][t], precision);
                output << " ";
            }
            output << std::endl;
        }
    }
    else {
        std::string separator = "  ...  ";

        // Print thetas (first row)
        output << Console::timePad << std::string(precision, ' ') << " | ";
        for (size_t t = 0; t < 4; ++t) {
            output << " ";
            output << Console::fixed_size(ndf.getTheta(t), precision);
            output << " ";
        }
        output << " " << separator << " ";
        for (size_t t = ndf.thetaSize() - 4; t < ndf.thetaSize(); ++t) {
            output << " ";
            output << Console::fixed_size(ndf.getTheta(t), precision);
            output << " ";
        }
        output << std::endl;

        // Header separation
        output << Console::timePad;
        for (int i = 0; i < 10; ++i) {
            output << std::string(precision, '-');
            if (i != 9) output << "--";
        }
        output << std::endl;

        // First 4 rows
        for (int p = 0; p < std::min(4, (int)ndf.phiSize()); ++p) {
            // theta
            output << Console::timePad << Console::fixed_size(ndf.getPhi(p), precision) << " | ";
            // first 4 thetas
            for (size_t t = 0; t < std::min(4, (int)ndf.thetaSize()); ++t) {
                output << " ";
                output << Console::fixed_size(ndf[p][t], precision);
                output << " ";
            }
            // space
            output << " " << separator << " ";
            // last 4 thetas
            for (size_t t = ndf.thetaSize() - std::min(4, (int)ndf.thetaSize()); t < ndf.thetaSize(); ++t) {
                output << " ";
                output << Console::fixed_size(ndf[p][t], precision);
                output << " ";
            }
            output << std::endl;
        }

        // Separator
        output << Console::timePad;
        output << " " << separator << "  |  ";
        for (int i = 1; i < 10; ++i) {
            if (i != 5) output << " " << separator << "   ";
            else        output << separator << "  ";
        }
        output << std::endl;

        // Last 4 rows
        for (size_t p = ndf.phiSize() - std::min(4, (int)ndf.phiSize()); p < ndf.phiSize(); ++p) {
            // theta
            output << Console::timePad << Console::fixed_size(ndf.getPhi(p), precision) << " | ";
            // first 4 thetas
            for (size_t t = 0; t < std::min(4, (int)ndf.thetaSize()); ++t) {
                output << " ";
                output << Console::fixed_size(ndf[p][t], precision);
                output << " ";
            }
            // space
            output << " " << separator << " ";
            // last 4 thetas
            for (size_t t = ndf.thetaSize() - std::min(4, (int)ndf.thetaSize()); t < ndf.thetaSize(); ++t) {
                output << " ";
                output << Console::fixed_size(ndf[p][t], precision);
                output << " ";
            }
            output << std::endl;
        }
    }

    output << std::endl;
    return output;
}