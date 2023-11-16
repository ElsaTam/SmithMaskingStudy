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

#pragma once

#include "shapes/TriangleMesh.h"

using namespace gdt;

class Histogram
{
private:
    scal heightMin;
    scal heightMax;
    scal hStep;
    std::vector<scal> pdf;
    std::vector<scal> cdf;

    int minSamples = 3;
    int maxSamples = 20;

    scal getHeight(int hi) const {
        return this->heightMin + this->hStep * (hi + 0.5);
    }
    int getIndex(scal h) const {
        return int((h - std::fmod(h, this->hStep)) / this->hStep);
    }

public:
    Histogram();
    Histogram(const vec3sc& w, const TriangleMesh& mesh, scal borderPercentage);
    ~Histogram() { ; }

    inline scal& operator[](size_t hIndex) { return pdf[hIndex]; }
    inline const scal& operator[](size_t hIndex) const { return pdf[hIndex]; }

    float sample(float u1) const;
    scal getPdf(scal h) const {
        const int hi = this->getIndex(h);
        return this->pdf[hi] / this->cdf[this->cdf.size() - 1];
    }
    const std::vector<scal>& getHist() const { return this->pdf; }

    size_t size() const { return 100; }

    friend std::ostream& operator<<(std::ostream& output, const Histogram& hist);
};

class HeightsDiscrete
{
private:
    std::vector<float> thetas; // thetas
    std::vector<float> phis; // phis
    std::vector<std::vector<Histogram>> H_values; // H_values[phi][theta]

    float getTheta(size_t i) const;
    float getPhi(size_t i) const;
    float getNextTheta(size_t i) const;
    float getNextPhi(size_t i) const;
    float getCenterTheta(size_t i) const;
    float getCenterPhi(size_t i) const;
    void findIndex(const vec3sc& w, size_t& thetaIndex, size_t& phiIndex) const;
    scal findValue(const std::vector<std::vector<float>>& V, const vec3sc& w) const;


public:
    HeightsDiscrete(const TriangleMesh& mesh, scal borderPercentage = 0);
    ~HeightsDiscrete() { ; }

    scal sample(const vec3sc& w, float u1) const; // return a random height for given direction w
    scal sample(const vec3sc& w, float u1, scal& _pdf) const;
    scal pdf(const vec3sc& w, float h) const;

    inline static size_t phiSize()   { return 400; } // 400
    inline static size_t thetaSize() { return 100; } // 100

    inline static float phiStart()   { return -m_pi; }
    inline static float phiEnd()     { return m_pi; }
    inline static float thetaStart() { return 0; }
    inline static float thetaEnd()   { return m_pi_2; }

    inline const std::vector<scal>& hist(const vec3sc& w) const {
        size_t tIdx, pIdx;
        findIndex(w, tIdx, pIdx);
        return H_values[pIdx][tIdx].getHist();
    }

    inline std::vector<Histogram>& operator[](size_t phiIndex) { return H_values[phiIndex]; }
    inline const std::vector<Histogram>& operator[](size_t phiIndex) const { return H_values[phiIndex]; }

    friend std::ostream& operator<<(std::ostream& output, const HeightsDiscrete& hdf);
};