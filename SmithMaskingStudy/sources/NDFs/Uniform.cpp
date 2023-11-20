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
#include "gdt/math/conversion.h"
#include "utils/console.h"

Uniform::Uniform()
{
    m_type = NDFTypes::UNIFORM;
}

/**
* @brief Normal distribution evaluation function
*
* This function must be statisfy :
* \f$ \int_{\Omega_+}D(\omega_h) cos(\theta_h) d\omega_h = 1 \f$
*
* @param wh        The microfacet normal (half-vector)
* @return scal     \f$  D(\omega_h) \f$
*/
scal Uniform::D(const vec3sc& wh) const
{
    return 1.f / m_pi;
}

/**
* @brief Compute the pdf for a given direction.
*
* pdf = w.z / m_pi
*
* @param u1        the first uniformly distributed sample used to compute w
* @param u2        the second uniformly distributed sample used to compute w
* @param w         the given direction
* @return scal     pdf(w)
*/
scal Uniform::pdf(const vec3sc& w, float u1, float u2) const
{
    return w.z / m_pi;
}

/**
* @brief Smith's shadowing-masking term for a single direction
*
* @param w         a direction vector
* @param wh        the microfacet normal
* @return scal    \f$ G(\omega,\omega_h) \f$
*/
scal Uniform::G1(const vec3sc& w, const vec3sc& wh) const
{
    Console::print(OutLevel::ERR, "Uniform::G1 not implemented.");
    exit(EXIT_FAILURE);
}

/**
* @brief Correlated Smith's shadowing-masking term
*
* @param wi        the incident direction
* @param wo        the outgoing direction
* @param wh        the microfacet normal
* @return scal    \f$ G(\omega_i, \omega_o, \omega_h) \f$
*/
scal Uniform::GAFcorrelated(const vec3sc& wi, const vec3sc& wo, const vec3sc& wh) const
{
    Console::print(OutLevel::ERR, "Uniform::GAFcorrelated not implemented.");
    exit(EXIT_FAILURE);
}


/**
* @brief Sample the microfacet normal distribution
*
* @param u1        an uniformly distributed sample
* @param u2        an uniformly distributed sample
* @return vec3sc    the sampled normal
*/
vec3sc Uniform::sample(float u1, float u2) const
{
    scal phi = u2 * m_2_pi;
    scal sin_theta = std::sqrt(1.f - u1);
    scal cos_theta = std::sqrt(u1);
    return vec3sc(std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, cos_theta);
}