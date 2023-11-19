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
#include "HDFs/HeightsDiscrete.hpp"

#include <algorithm>
#include <unordered_set>
#include <random>
#include <cmath>
#include "utils/console.h"
#include "utils/params.h"
#include "utils/math/math.h"
#include "tools/logger.h"



Histogram::Histogram()
    : heightMin(0), heightMax(0), hStep(0)
{
    this->pdf.resize(this->size(), 0);
    this->cdf.resize(this->size(), 0);
}

Histogram::Histogram(const vec3sc& w, const TriangleMesh& mesh, scal borderPercentage)
    : heightMin(mesh.bounds.lower.z), heightMax(mesh.bounds.upper.z)
{
    this->hStep = mesh.bounds.span().z / this->size();

    this->pdf.resize(this->size(), 0);
    this->cdf.resize(this->size(), 0);

    std::default_random_engine rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<scal> rng(0.f, 1.f);

    for (int faceID = 0; faceID < mesh.index.size(); ++faceID) {
        const vec3i& idx = mesh.index[faceID];

        /* Get the projection weighting factor */
        const vec3sc& n  = mesh.triangle_normal[faceID];
        const scal weight = dot(w, n);
        if (weight <= 0.) {
            continue;
        }

        /* Uniform sampling from unit-square */
        const vec3sc& v1 = mesh.vertex[idx[0]];
        const vec3sc& v2 = mesh.vertex[idx[1]];
        const vec3sc& v3 = mesh.vertex[idx[2]];

        /* Get number of samples for this face */
        const scal A = mesh.area[faceID];
        const int nSamples = (A - 0) / (mesh.maxArea - 0) * (maxSamples - minSamples) + minSamples;

        for (int i = 0; i < nSamples; ++i) {
            scal s1 = rng(gen);
            scal s2 = rng(gen);
            /* Mapping to unit-triangle */
            scal t1 = 0.5f * s1;
            scal t2 = 0.5f * s2;
            scal offset = t2 - t1;
            if (offset > 0) { t2 += offset; }
            else { t1 -= offset; }

            /* Mapping to arbitrary triangle */
            vec3sc v(t1 * v1 + t2 * v2 + (1.f - t1 - t2) * v3);

            /* Check if not too close from the border */
            if (mesh.bounds.closest_distance(v).x < Parameters::get()->currentParams()->sideEffectParams.borderPercentage * mesh.bounds.span().x / 2.f
                || mesh.bounds.closest_distance(v).y < Parameters::get()->currentParams()->sideEffectParams.borderPercentage * mesh.bounds.span().y / 2.f)
            {
                continue;
            }

            /* Add correct value in the histogram */
            this->pdf[this->getIndex(v3.z)] += weight;
        }
    }

    this->cdf[0] = this->pdf[0];
    for (int hi = 1; hi < this->pdf.size(); ++hi)
    {
        this->cdf[hi] = this->cdf[hi - 1] + this->pdf[hi];
    }
}

float Histogram::sample(float u1) const
{
    u1 *= this->cdf[this->cdf.size() - 1];

    std::vector<scal>::const_iterator entry = std::lower_bound(
        this->cdf.begin(),
        this->cdf.end(),
        u1,
        [](const scal& info, double value) {
            return info < value;
        }
    );
    size_t hi = std::min(cdf.size() - 2, (size_t)std::max((ptrdiff_t)0, entry - cdf.begin() - 1));
    // Handle a rare corner-case where a entry has probability 0 but is sampled nonetheless
    while (hi < cdf.size() - 1 && cdf[hi + 1] - cdf[hi] == 0) {
        ++hi;
    }

    return getHeight(hi);
}

std::ostream& operator<<(std::ostream& output, const Histogram& hist) {
    output << Console::timePad << "Histogram: size (nHeights) = " << hist.size() << std::endl;
    output << Console::timePad << Console::indent << "from " << hist.heightMin << " to " << hist.heightMax << std::endl;
    output << std::endl;

    // Print distribution
    int precision = 9;

    std::string separator = "  ...  ";

    // first 4 pdfs
    for (size_t hi = 0; hi < std::min(4, (int)hist.size()); ++hi) {
        output << " ";
        output << Console::fixed_size(hist.pdf[hi], precision);
        output << " ";
    }
    // space
    output << " " << separator << " ";
    // last 4 pdfs
    for (size_t hi = hist.size() - std::min(4, (int)hist.size()); hi < hist.size(); ++hi) {
        output << " ";
        output << Console::fixed_size(hist.pdf[hi], precision);
        output << " ";
    }

    output << std::endl;
    return output;
}












float HeightsDiscrete::getTheta(size_t i) const {
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

float HeightsDiscrete::getPhi(size_t i) const {
    i = i % phis.size();
    return phis[i];
}

float HeightsDiscrete::getNextTheta(size_t i) const {
    if (i >= thetas.size() - 1) return 2 * thetas[i] - thetas[i - 1]; // assume constant dTheta
    return thetas[i + 1];
}

float HeightsDiscrete::getNextPhi(size_t i) const {
    i = i % phis.size();
    if (i == phis.size() - 1) return phis[0] + m_2_pi;
    return phis[i + 1];
}

float HeightsDiscrete::getCenterTheta(size_t i) const {
    if (i >= thetas.size() - 1) return 0;
    return 0.5f * (getNextTheta(i) + getTheta(i));
}

float HeightsDiscrete::getCenterPhi(size_t i) const {
    if (i >= phis.size()) {
        Console::err << "Error: trying to access phi at index " << i << " (over " << phis.size() << ")" << std::endl;
        return 0;
    }
    return 0.5f * (getNextPhi(i) + getPhi(i));
}

void HeightsDiscrete::findIndex(const vec3sc& w, size_t& thetaIndex, size_t& phiIndex) const
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

scal HeightsDiscrete::findValue(const std::vector<std::vector<float>>& V, const vec3sc& w) const
{
    size_t tIdx, pIdx;
    findIndex(w, tIdx, pIdx);
    return V[pIdx][tIdx];
}



HeightsDiscrete::HeightsDiscrete(const TriangleMesh& mesh, scal borderPercentage)
{
    Console::out << Console::timeStamp << "Building Discrete HDF from mesh..." << std::endl;

    scal phiStep = (phiEnd() - phiStart()) / (scal)phiSize();
    scal thetaStep = (thetaEnd() - thetaStart()) / (scal)thetaSize();

    for (int p = 0; p < phiSize(); ++p) {
        const scal phi = phiStart() + p * phiStep;
        phis.push_back(phi);
    }
    for (int t = 0; t < thetaSize(); ++t) {
        const scal theta = thetaStart() + t * thetaStep;
        thetas.push_back(theta);
    }

    H_values.resize(phiSize());
    for (int p = 0; p < phiSize(); ++p) {
        Console::out << Console::timePad << p+1 << "/" << phiSize() << std::endl;
        for (int t = 0; t < thetaSize(); ++t) {
            const vec3sc dir = Conversion::polar_to_cartesian(getTheta(t), getPhi(p));
            H_values[p].push_back(Histogram(dir, mesh, borderPercentage));
        }
    }

    Console::out << *this << std::endl;
}


/**
* @brief Compute the pdf for a given direction.
*
* @param u1        the first uniformly distributed sample used to compute w
* @param u2        the second uniformly distributed sample used to compute w
* @param w         the given direction
* @return scal     pdf(w)
*/
scal HeightsDiscrete::pdf(const vec3sc& w, float h) const
{
    const vec3sc w_polar = Conversion::cartesian_to_polar(w);
    float theta = w_polar[0]; // acos in range [0,pi]
    float phi = w_polar[1]; // atan2 in range [-pi,+pi]
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
    size_t pL = C.x > phi ? L : cLR;
    size_t pR = C.x > phi ? cLR : R;
    size_t tT = C.y > theta ? T : cTB;
    size_t tB = C.y > theta ? cTB : B;

    // phi1 < phi < phi2  &&  theta1 < theta < theta2
    float phi1 = getCenterPhi(pL);
    float phi2 = getCenterPhi(pR);
    float theta1 = getCenterTheta(tB);
    float theta2 = getCenterTheta(tT);

    scal hLB = H_values[pL][tB].getPdf(h);
    scal hLT = H_values[pL][tT].getPdf(h);
    scal hRB = H_values[pR][tB].getPdf(h);
    scal hRT = H_values[pR][tT].getPdf(h);
    float BI = Maths::bilinearInterpolation(
        hLB, hLT, hRB, hRT,
        phi1, phi2,
        theta1, theta2,
        phi, theta);
    //float BI_0 = H_values[pIdx][tIdx]->getPdf(h); // for debug purpose
    return BI;
}


/**
* @brief Sample the microfacet normal distribution
*
* @param u1        an uniformly distributed sample
* @param u2        an uniformly distributed sample
* @return vec3sc    the sampled normal
*/
scal HeightsDiscrete::sample(const vec3sc& w, float u1) const
{
    size_t tIdx, pIdx;
    findIndex(w, tIdx, pIdx);

    return H_values[pIdx][tIdx].sample(u1);
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
scal HeightsDiscrete::sample(const vec3sc& w, float u1, scal& _pdf) const
{
    size_t tIdx, pIdx;
    findIndex(w, tIdx, pIdx);
    const scal h = H_values[pIdx][tIdx].sample(u1);
    _pdf = H_values[pIdx][tIdx].getPdf(h);
    return h;
}




std::ostream& operator<<(std::ostream& output, const HeightsDiscrete& hdf) {
    output << Console::timePad << "Heights Discrete distribution: shape (nPhi, nTheta) = (" << hdf.phiSize() << ", " << hdf.thetaSize() << ")" << std::endl;
    output << Console::timePad << Console::indent << "Phi:   from " << hdf.phiStart() << " to " << hdf.phiEnd() << std::endl;
    output << Console::timePad << Console::indent << "Theta: from " << hdf.thetaStart() << " to " << hdf.thetaEnd() << std::endl;
    output << Console::timePad << Console::indent << hdf[0][0] << std::endl;
    output << std::endl;
    return output;
}