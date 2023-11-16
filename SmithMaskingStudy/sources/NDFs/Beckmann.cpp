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

Beckmann::Beckmann(scal alpha_x, scal alpha_y)
    :m_alpha_x(alpha_x), m_alpha_y(alpha_y)
{
    m_type = NDFTypes::BECKMANN;
}

/**
 * @brief Evaluate the beckmann p22
 *
 * The Beckmann slope distribution is a normalized normal distribution using the RMS roughness.
 * \f$ \alpha \f$ and \f$ \sigma \f$ are related by \f$ \sigma = \frac{\alpha}{\sqrt{2}} \f$
 *
 * \f$ P_{22}(\tilde{x_h},\tilde{y_h};\sigma_x,\sigma_y) = \frac{1}{2\pi\sigma_x\sigma_y} \exp(-\frac{1}{2} (\frac{x^2}{\sigma_x^2} + \frac{y^2}{\sigma_y^2})) \f$
 *
 * properties:
 * \int_{\mathbb{R}^2} P_{22}(\tilde{w_h}) d\tilde{w_h} = 1
 *
 * @param x         the x slope
 * @param y         the y slope
 * @return scal    the slope distribution value
 */
scal Beckmann::P22(scal x, scal y) const
{
    scal x_sqr = x * x;
    scal y_sqr = y * y;
    scal sigma_x = m_alpha_x * m_i_sqrt_2;
    scal sigma_y = m_alpha_y * m_i_sqrt_2;
    scal sigma_x_sqr = sigma_x * sigma_x;
    scal sigma_y_sqr = sigma_y * sigma_y;

    return(
        std::exp(-0.5f * ((x_sqr / sigma_x_sqr) + (y_sqr / sigma_y_sqr)))
        / //-------------------------------------------------------------------
        (2.f * m_pi * sigma_x * sigma_y)
        );
}

/**
 * @brief Beckmann lambda function
 *
 * \f$ \Lambda(\omega) = \tan\theta\int^{\infty}_{\cot\theta}(\tilde{x_h} - \cot\theta) \Big( \int_{-\infty}^{+\infty} P_{22}(\tilde{x_h},\tilde{y_h})d\tilde{y_h} \Big) d\tilde{x_h} \f$
 * \f$ \Lambda(\omega) = \frac{\exp(-\nu^2)}{2\nu\sqrt{\pi}} - \frac{\text{erf}(\nu)}{2} \fs
 * \f$ \Lambda(\omega) \approx \begin{cases} \frac{1.0 - 1.259\nu+0.396\nu^2}{3.535\nu + 2.181\nu^2} & \text{if}~~\nu < 1.6 \\ 0 &\text{otherwise}\end{cases} \fs
 * \f$ \nu = \frac{1}{\tan\theta\sqrt{\cos(\phi)^2\alpha_x^2 + \sin(\phi)^2\alpha_y^2}} \f$
 *
 * @param omega
 * @return scal
 */
scal Beckmann::lambda(const vec3sc& omega) const
{
    scal lambda = 0.f;
    scal tan_theta = Frame3sc::tan_theta(omega);
    scal alpha_phi = Frame3sc::projected_roughness(omega, m_alpha_x, m_alpha_y);

    scal nu = 1.0f / (alpha_phi * tan_theta);
    if (nu < 1.6f) { lambda = (1.0f - 1.259f * nu + 0.396f * nu * nu) / (3.535f * nu + 2.181f * nu * nu); }
    return(lambda);
}

/**
 * @brief Normal distribution evaluation function
 *
 * Using the slope distribution relation :
 * \f$ D(\omega_h) = P_{22}(\tilde{\omega_{h}}) sec^{4}(\theta_h) \f$
 * \f$ D(\omega_h) = \frac{ P_{22}(\tilde{\omega_{h}}) } { cos^4(\theta_h) } \f$
 * normals \f$ \omega_h \in \Omega_+ \f$ and slopes \f$ \tilde{\omega_h}  \in \mathbb{R}^2 \f$ are linked by the bijection :
 * \f$ \tilde{\omega_h} = ~(-tan(\theta_h)cos(\phi_h), -tan(\theta_h)sin(\phi_h) )^t = ( \tilde{x_h} , \tilde{y_h} )^t  \f$
 * whose inverse is:
 * \f$ \omega_h =\frac{1}{\sqrt{1+\tilde{x_h}^2+\tilde{y_h}^2}} ( -\tilde{x_h} , -\tilde{y_h}, 1 )^t  \f$
 * in tangent space :
 * \f$ \tilde{x_h} = -\cos(\phi)\tan(\theta) =  \frac{-\cos(\phi)\sin(\theta)}{\cos(\theta)} = \frac{-\omega_h.x}{\omega_h.z} \f$
 * \f$ \tilde{y_h} = -\sin(\phi)\tan(\theta) =  \frac{-\sin(\phi)\sin(\theta)}{\cos(\theta)} = \frac{-\omega_h.y}{\omega_h.z} \f$
 *
 * properties :
 * \f$ \int_{\Omega_+}D(\omega_h) cos(\theta_h) d\omega_h = 1 \f$
 *
 * @param wh        The microfacet normal (half-vector)
 * @return scal    \f$  D(\omega_h) \f$
 */
scal Beckmann::D(const vec3sc& wh) const
{
    scal cos_theta = wh.z;
    if (cos_theta <= 0.)
        return(0.);

    scal slope_x = -(wh.x / wh.z);
    scal slope_y = -(wh.y / wh.z);
    scal cos_2_theta = cos_theta * cos_theta;
    scal cos_4_theta = cos_2_theta * cos_2_theta;
    return(
        this->P22(slope_x, slope_y)
        / //-----------------------------------
        cos_4_theta
        );
}

/**
* @brief Compute the pdf for a given direction.
*
* pdf = (1.0f - u2) / (m_pi * m_alpha_x * m_alpha_y * w.z * w.z * w.z)
*
* @param u1        the first uniformly distributed sample used to compute w
* @param u2        the second uniformly distributed sample used to compute w
* @param w         the given direction
* @return scal     pdf(w)
*/
scal Beckmann::pdf(const vec3sc& w, float u1, float u2) const
{
    return (1.0f - u2) / (m_pi * m_alpha_x * m_alpha_y * w.z * w.z * w.z);
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
scal Beckmann::G1(const vec3sc& w, const vec3sc& wh) const
{
    return(Real::min(1.f / (1.f + lambda(w)), 1.f));
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
scal Beckmann::GAFcorrelated(const vec3sc& wi, const vec3sc& wo, const vec3sc& wh) const
{
    return(Real::min(1.f / (1.f + lambda(wi) + lambda(wo)), 1.f));
}


/**
 * @brief Sample the Beckmann normal distribution
 *
 * Using a sampled vector \f$ u = (s_1,s_2)^t \f$
 *
 * For isotropic Beckmann :
 * \f$ \phi = 2\pi \s_1  \f$
 * \f$ \theta = \arctan(\sqrt{-\alpha^2 \ln(1-s_2)}) \f$
 * For anisotropic Beckmann :
 * \f$ \phi =  \phi = \arctan( \frac{\alpha_y}{\alpha_x} \tan(\pi + 2\pi s_1) \f$
 * \f$ \theta = \arctan( \frac{\ln(s_2)}{ \frac{\cos(\phi)^2}{\alpha_x^2} + \frac{\sin(\phi)^2}{\alpha_y^2} } ) \f$
 *
 * @param u         IN: an uniformly distributed 2D sample
 * @return vec3sc    the sampled normal
 */
vec3sc Beckmann::sample(float u1, float u2) const
{
    scal theta = 0.f;
    scal phi = 0.f;

    if (Real::equivf(m_alpha_x, m_alpha_y))
    {
        scal alpha_sqr = m_alpha_x * m_alpha_x;
        phi = m_2_pi * u1;
        theta = std::atan(std::sqrt(-alpha_sqr * std::log(1.f - u2)));
    }
    else
    {
        phi = std::atan(m_alpha_y / m_alpha_x * std::tan(m_pi + m_2_pi * u1));
        phi += m_pi * std::floor(2.f * u1 + 0.5f);

        scal cos_phi = std::cos(phi);
        scal sin_phi = std::sin(phi);
        scal cos_phi_over_alpha = cos_phi / m_alpha_x;
        scal sin_phi_over_alpha = sin_phi / m_alpha_y;
        scal alpha_phi = cos_phi_over_alpha * cos_phi_over_alpha + sin_phi_over_alpha * sin_phi_over_alpha;
        scal tan_theta_sqr = -std::log(1.f - u2) / alpha_phi;
        theta = std::atan(std::sqrt(tan_theta_sqr));
    }
    return Conversion::polar_to_cartesian(theta, phi);
}