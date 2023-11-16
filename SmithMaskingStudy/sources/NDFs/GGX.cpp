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

GGX::GGX(scal alpha_x, scal alpha_y)
    :m_alpha_x(alpha_x),m_alpha_y(alpha_y)
{
    m_type = NDFTypes::GGX;
}

/**
 * @brief GGX lambda function
 *  
 * \Lambda(\omega) = \frac{ -1 + \sqrt{1+\alpha(\phi)^2\tan(\theta)^2} }{2}
 * with
 * \f$ \alpha(\phi) = \sqrt{ \cos(\phi)^2\alpha_x^2 + \sin(\phi)^2\alpha_y^2 } \f$
 * 
 * @param omega 
 * @return scal 
 */
scal GGX::lambda(const vec3sc & omega) const
{
    scal alpha_phi = Frame3sc::projected_roughness(omega,m_alpha_x,m_alpha_y);
    scal nu = 1.f / (alpha_phi * Frame3sc::tan_theta(omega));
    scal lambda = 0.5f * (-1.f + sqrt(1.f + 1.f / (nu*nu))); 
    return(lambda);
}

/**
 * @brief Normal distribution evaluation function
 * 
 * This function must be statisfy :
 * \f$ \int_{\Omega_+}D(\omega_h) cos(\theta_h) d\omega_h = 1 \f$
 * 
 * \f$ D(\omega_h) = \frac{1}{\pi\alpha_x\alpha_y\cos(\theta_h)^4 \big(1+\tan(\theta_h)^2 ( \frac{\cos(\phi)^2}{\alpha_x^2} + \frac{\sin(\phi_h)^2}{\alpha_y^2} ) \big)^2 } \f$
 * 
 * @param wh        The microfacet normal (half-vector)
 * @return scal    \f$  D(\omega_h) \f$   
 */
scal GGX::D(const vec3sc & wh) const
{
    scal cos_theta = Frame3sc::cos_theta(wh);
    if(cos_theta <= 0.f)
        return(0.f);

    scal tan_theta = Frame3sc::tan_theta(wh);
    scal cos_phi = Frame3sc::cos_phi(wh);
    scal sin_phi = Frame3sc::sin_phi(wh);
    scal tan_theta_sqr = tan_theta * tan_theta;
    scal cos_theta_sqr = cos_theta * cos_theta;
    scal exponent = ( cos_phi*cos_phi/(m_alpha_x*m_alpha_x) + sin_phi*sin_phi/(m_alpha_y*m_alpha_y) ) * tan_theta_sqr;
    scal root = (1.f + exponent) * cos_theta_sqr;
    return(
                            1.f
    / //---------------------------------------------------            
        (m_pi * m_alpha_x * m_alpha_y * root * root)
    );
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
scal GGX::pdf(const vec3sc& w, float u1, float u2) const
{
    scal alpha_phi;
    if (Real::equivf(m_alpha_x, m_alpha_y))
    {
        alpha_phi = m_alpha_x * m_alpha_x;
    }
    else
    {
        scal phi = std::atan(m_alpha_y / m_alpha_x * std::tan(m_pi + 2.0f * m_pi * u1));
        phi += m_pi * std::floor(2.0f * u1 + 0.5f);
        scal cos_phi_over_alpha = std::cos(phi) / m_alpha_x;
        scal sin_phi_over_alpha = std::sin(phi) / m_alpha_y;
        alpha_phi = cos_phi_over_alpha * cos_phi_over_alpha + sin_phi_over_alpha * sin_phi_over_alpha;
    }

    scal tan_theta_sqr = u2 / ((1.0f - u2) * alpha_phi);
    scal cos_theta = 1.0f / std::sqrt(1.0f + tan_theta_sqr);
    scal tmp = 1.f + tan_theta_sqr * alpha_phi;
    return 1.f / (m_pi * m_alpha_x * m_alpha_y * cos_theta * cos_theta * cos_theta * tmp * tmp);
}

/**
 * @brief Smith's shadowing-masking term for a single direction
 * 
 * \f$ \text{G1}(\omega) = \frac{1}{1 + \Lambda(\omega)}  \f$
 * 
 * @param w         a direction vector
 * @param wh        the microfacet normal
 * @return scal    \f$ G(\omega,\omega_h) \f$
 */
scal GGX::G1(const vec3sc & w, const vec3sc & wh) const
{
    return( Real::min(1.f/(1.f+lambda(w)),1.f) );
}

/**
 * @brief Correlated Smith's shadowing-masking term
 * 
 * \f$ \text{GAF}(\omega) = \frac{1}{1 + \Lambda(\omega_o) + \Lambda(\omega_i)}  \f$
 * 
 * @param wi        the incident direction
 * @param wo        the outgoing direction
 * @param wh        the microfacet normal
 * @return scal    \f$ G(\omega_i, \omega_o, \omega_h) \f$
 */
scal GGX::GAFcorrelated(const vec3sc & wi, const vec3sc & wo, const vec3sc & wh) const
{
    return( Real::min(1.f/(1.f+lambda(wi)+lambda(wo)),1.f) );
}

/**
 * @brief Sample the GGX microfacet normal distribution
 * 
 * Using a sampled vector \f$ u = (s_1,s_2)^t \f$ 
 * and \f$ A(\phi) = \frac{\cos(\phi)^2}{\alpha_x^2} + \frac{\sin(\phi)^2}{\alpha_y^2} \f$ 
 * 
 * For isotropic GGX :
 * \f$ \phi = 2\pi \s_1  \f$
 * \f$ \theta = \arctan( \alpha \sqrt{\frac{s_2}{1-s_2}}) \f$
 * For anisotropic GGX :
 * \f$ \phi =  \f$
 * \f$ \theta = \arctan \big( \sqrt{\frac{s_2}{(1-s_2)A(\phi)}} \big) \f$
 * 
 * @param u         IN: an uniformly distributed 2D sample
 * @return vec3sc    the sampled normal
 */
vec3sc GGX::sample(float u1, float u2) const
{
    scal theta = 0.f;
    scal phi = 0.f;
    if(Real::equivf(m_alpha_x,m_alpha_y))
    {
        phi = m_2_pi * u1;
        theta = std::atan(m_alpha_x * std::sqrt(u2 / (1.f-u2) ));
    }
    else
    {
        phi = std::atan(m_alpha_y / m_alpha_x * std::tan(m_pi + 2.0f*m_pi*u1));
        phi += m_pi * std::floor( 2.0f*u1 + 0.5f );
        
        scal cos_phi = std::cos(phi);
        scal sin_phi = std::sin(phi);
        scal cos_phi_over_alpha = cos_phi / m_alpha_x;
        scal sin_phi_over_alpha = sin_phi / m_alpha_y;
        scal alpha_phi = cos_phi_over_alpha*cos_phi_over_alpha + sin_phi_over_alpha*sin_phi_over_alpha;
        theta = std::atan( std::sqrt(u2 / ((1.0f-u2)*alpha_phi) ) );
    }
    return Conversion::polar_to_cartesian(theta,phi);
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
vec3sc GGX::sample(float u1, float u2, scal& _pdf, scal& _D) const
{
    scal theta = 0.f;
    scal phi = 0.f;
    scal alpha_phi, cos_phi, sin_phi;
    if (Real::equivf(m_alpha_x, m_alpha_y))
    {
        phi = m_2_pi * u1;
        cos_phi = std::cos(phi);
        sin_phi = std::sin(phi);
        alpha_phi = m_alpha_x * m_alpha_x;
    }
    else
    {
        phi = std::atan(m_alpha_y / m_alpha_x * std::tan(m_pi + 2.0f * m_pi * u1));
        phi += m_pi * std::floor(2.0f * u1 + 0.5f);

        cos_phi = std::cos(phi);
        sin_phi = std::sin(phi);
        scal cos_phi_over_alpha = cos_phi / m_alpha_x;
        scal sin_phi_over_alpha = sin_phi / m_alpha_y;
        alpha_phi = cos_phi_over_alpha * cos_phi_over_alpha + sin_phi_over_alpha * sin_phi_over_alpha;
    }

    scal tan_theta_sqr = u2 / ((1.0f - u2) * alpha_phi);
    scal cos_theta = 1.0f / std::sqrt(1.0f + tan_theta_sqr);
    scal tmp = 1.f + tan_theta_sqr * alpha_phi;
    _pdf = 1.f / (m_pi * m_alpha_x * m_alpha_y * cos_theta * cos_theta * cos_theta * tmp * tmp);

    scal cos_theta_sqr = cos_theta * cos_theta;
    scal exponent = (cos_phi * cos_phi / (m_alpha_x * m_alpha_x) + sin_phi * sin_phi / (m_alpha_y * m_alpha_y)) * tan_theta_sqr;
    scal root = (1.f + exponent) * cos_theta_sqr;
    _D = 1.f / (m_pi * m_alpha_x * m_alpha_y * root * root);

    theta = std::atan(std::sqrt(tan_theta_sqr));
    return Conversion::polar_to_cartesian(theta, phi);
}