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

#include "gdt/math/frame.h"
#include "shapes/TriangleMesh.h"

using namespace gdt;

enum class NDFTypes {
    BECKMANN,
    GGX,
    UNIFORM,
    DISCRETE,
    UNKNOWN
};

/**
 * @brief Microfacet distribution interface
 */
class MicrofacetDistribution
{
    protected:
        NDFTypes m_type = NDFTypes::UNKNOWN;

    public:      
        virtual ~MicrofacetDistribution()
        {;}

        /**
         * @brief Normal distribution evaluation function
         * 
         * This function must be statisfy :
         * \f$ \int_{\Omega_+}D(\omega_h) cos(\theta_h) d\omega_h = 1 \f$
         * 
         * @param wh        The microfacet normal (half-vector)
         * @return scal     \f$  D(\omega_h) \f$   
         */
        virtual scal D(const vec3sc & wh) const = 0;

        /**
         * @brief Compute the pdf for a given direction.
         *
         * @param u1        the first uniformly distributed sample used to compute w
         * @param u2        the second uniformly distributed sample used to compute w
         * @param w         the given direction
         * @return scal     pdf(w)
        */
        virtual scal pdf(const vec3sc& w, float u1, float u2) const = 0;

        /**
         * @brief Smith's shadowing-masking term for a single direction
         * 
         * @param w         a direction vector
         * @param wh        the microfacet normal
         * @return scal    \f$ G(\omega,\omega_h) \f$
         */
        virtual scal G1(const vec3sc & w, const vec3sc & wh) const = 0;
        
        /**
         * @brief Correlated Smith's shadowing-masking term
         * 
         * @param wi        the incident direction
         * @param wo        the outgoing direction
         * @param wh        the microfacet normal
         * @return scal    \f$ G(\omega_i, \omega_o, \omega_h) \f$
         */
        virtual scal GAFcorrelated(const vec3sc & wi, const vec3sc & wo, const vec3sc & wh) const = 0;
               
        /**
         * @brief Sample the microfacet normal distribution
         * 
         * @param u1        an uniformly distributed sample
         * @param u2        an uniformly distributed sample
         * @return vec3sc    the sampled normal
         */
        virtual vec3sc sample(float u1, float u2) const = 0;

        /**
         * @brief Sample the microfacet normal distribution
         *
         * @param u1        an uniformly distributed sample
         * @param u2        an uniformly distributed sample
         * @param _pdf       out parameter, for the pdf of the sampled normal
         * @param _D         out parameter, for the distribution of the sampled normal
         * @return vec3sc    the sampled normal
         */
        virtual inline vec3sc sample(float u1, float u2, scal& _pdf, scal& _D) const
        {
            vec3sc w = sample(u1, u2);
            _pdf = pdf(w, u1, u2);
            _D = D(w);
            return w;
        }

        inline NDFTypes type() const { return m_type; }
        virtual std::string name() const = 0;
};




/**
 * @brief Beckmann Distribution 
 * From : (1963) The Scattering of Electromagnetic Waves from Rough Surfaces. 
 * By   : Beckmann, P., and A. Spizzichino. 1963.
 */
class Beckmann : public MicrofacetDistribution
{
    private:
        scal m_alpha_x; /**< Roughness in the tangent x direction */
        scal m_alpha_y; /**< Roughness in the tangent y direction */

    public:
        Beckmann(scal alpha_x, scal alpha_y);
        ~Beckmann() {;}

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
        scal P22(scal x, scal y) const;

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
        scal lambda(const vec3sc & omega) const;

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
        scal D(const vec3sc & wh) const override;

        /**
         * @brief Compute the pdf for a given direction.
         *
         * pdf = (1.0f - u2) / (m_pi * 0.5 * 0.5 * w.z * w.z * w.z)
         *
         * @param u1        the first uniformly distributed sample used to compute w
         * @param u2        the second uniformly distributed sample used to compute w
         * @param w         the given direction
         * @return scal     pdf(w)
        */
        scal pdf(const vec3sc& w, float u1, float u2) const override;

        /**
         * @brief Smith's shadowing-masking term for a single direction
         * 
         * \f$ \text{G1}(\omega) = \frac{1}{1 + \Lambda(\omega)}  \f$
         * 
         * @param w         a direction vector
         * @param wh        the microfacet normal
         * @return scal    \f$ G(\omega,\omega_h) \f$
         */
        scal G1(const vec3sc & w, const vec3sc & wh) const override;

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
        scal GAFcorrelated(const vec3sc & wi, const vec3sc & wo, const vec3sc & wh) const override;

        /**
         * @brief Sample the Beckmann normal distribution
         * 
         * Using a sampled vector \f$ u = (s_1,s_2)^t \f$ 
         * 
         * For isotropic Beckmann :
         * \f$ \phi = 2\pi \s_1  \f$
         * \f$ \theta = \arctan(\sqrt{-\alpha^2 \ln(1-s_2)}) \f$
         * For anisotropic Beckmann :
         * \f$ \phi =  \phi = \arctan( \frac{\alpha_y}{\alpha_u} \tan(\pi + 2\pi s_1) \f$
         * \f$ \theta = \arctan( \frac{\ln(s_2)}{ \frac{\cos(\phi)^2}{\alpha_x^2} + \frac{\sin(\phi)^2}{\alpha_y^2} } ) \f$
         * 
         * @param u         IN: an uniformly distributed 2D sample
         * @return vec3sc    the sampled normal
         */
        vec3sc sample(float u1, float u2) const override;

        inline std::string name() const override { return "Beckmann"; }
};


/**
 * @brief GGX or Trowbridge-Reitz Distribution 
 * 
 * From : (1975) Average irregularity representation of a rough ray reflection.
 * By   : Trowbridge, S., and K. P. Reitz.
 * 
 * From : (2007) Microfacet models for refraction through rough surfaces. 
 * By   : Walter, B., S. Marschner, H. Li, and K. Torrance.
 */
class GGX : public MicrofacetDistribution
{
    private:
        scal m_alpha_x; /**< Roughness in the tangent x direction */
        scal m_alpha_y; /**< Roughness in the tangent y direction */
    public:
        GGX(scal alpha_x, scal alpha_y);
        ~GGX() {;}

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
        scal lambda(const vec3sc & omega) const;

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
        scal D(const vec3sc & wh) const override;

        /**
         * @brief Compute the pdf for a given direction.
         *
         * pdf = ???
         *
         * @param u1        the first uniformly distributed sample used to compute w
         * @param u2        the second uniformly distributed sample used to compute w
         * @param w         the given direction
         * @return scal     pdf(w)
        */
        scal pdf(const vec3sc& w, float u1, float u2) const override;

        /**
         * @brief Smith's shadowing-masking term for a single direction
         * 
         * \f$ \text{G1}(\omega) = \frac{1}{1 + \Lambda(\omega)}  \f$
         * 
         * @param w         a direction vector
         * @param wh        the microfacet normal
         * @return scal    \f$ G(\omega,\omega_h) \f$
         */
        scal G1(const vec3sc & w, const vec3sc & wh) const override;

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
        scal GAFcorrelated(const vec3sc & wi, const vec3sc & wo, const vec3sc & wh) const override;

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
        vec3sc sample(float u1, float u2) const override;

        /**
         * @brief Sample the microfacet normal distribution
         *
         * @param u1        an uniformly distributed sample
         * @param u2        an uniformly distributed sample
         * @param _pdf       out parameter, for the pdf of the sampled normal
         * @param _D         out parameter, for the distribution of the sampled normal
         * @return vec3sc    the sampled normal
         */
        vec3sc sample(float u1, float u2, scal& _pdf, scal& _D) const override;

        inline std::string name() const override { return "GGX"; }
};


class Uniform : public MicrofacetDistribution
{
public:
    Uniform();
    ~Uniform() { ; }

    /**
    * @brief Normal distribution evaluation function
    *
    * This function must be statisfy :
    * \f$ \int_{\Omega_+}D(\omega_h) cos(\theta_h) d\omega_h = 1 \f$
    *
    * @param wh        The microfacet normal (half-vector)
    * @return scal     \f$  D(\omega_h) \f$
    */
    scal D(const vec3sc& wh) const override;

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
    scal pdf(const vec3sc& w, float u1 = 0, float u2 = 0) const override;

    /**
    * @brief Smith's shadowing-masking term for a single direction
    *
    * @param w         a direction vector
    * @param wh        the microfacet normal
    * @return scal    \f$ G(\omega,\omega_h) \f$
    */
    scal G1(const vec3sc& w, const vec3sc& wh) const override;

    /**
    * @brief Correlated Smith's shadowing-masking term
    *
    * @param wi        the incident direction
    * @param wo        the outgoing direction
    * @param wh        the microfacet normal
    * @return scal    \f$ G(\omega_i, \omega_o, \omega_h) \f$
    */
    scal GAFcorrelated(const vec3sc& wi, const vec3sc& wo, const vec3sc& wh) const override;

    /**
    * @brief Sample the microfacet normal distribution
    *
    * @param u1        an uniformly distributed sample
    * @param u2        an uniformly distributed sample
    * @return vec3sc    the sampled normal
    */
    vec3sc sample(float u1, float u2) const override;

    inline std::string name() const override { return "Uniform"; }
};


class Discrete : public MicrofacetDistribution
{
private:
    std::vector<float> thetas; // thetas Out
    std::vector<float> phis; // phis Out
    std::vector<std::vector<float>> D_values; // D_values[phi][theta]
    std::vector<std::vector<float>> pdf_values; // pdf_values[phi][theta]
    std::vector<std::tuple<int, int, float>> cdf; // for sampling <thetaIdx, phiIdx, D>
    float cdf_sum;
    const MicrofacetDistribution* matchingNDF;
    vec3sc m_meso_normal;
    vec3sc m_average_normal;
    scal m_D_integral;

    float getTheta(size_t i) const;
    float getPhi(size_t i) const;
    float getNextTheta(size_t i) const;
    float getNextPhi(size_t i) const;
    float getCenterTheta(size_t i) const;
    float getCenterPhi(size_t i) const;
    void findIndex(const vec3sc& w, size_t& thetaIndex, size_t& phiIndex) const;
    scal findInterpolatedValue(const std::vector<std::vector<float>>& V, const vec3sc& w) const;
    scal findValue(const std::vector<std::vector<float>>& V, const vec3sc& w) const;

    void addNormal(const vec3sc& n, const scal weight = 1);

public:
    Discrete(const TriangleMesh& mesh, const MicrofacetDistribution* _matchingNDF = nullptr, scal borderPercentage = 0);
    Discrete(const MicrofacetDistribution* analyticNDF, int samples);
    ~Discrete() { ; }

    /**
    * @brief Normal distribution evaluation function
    *
    * This function must be statisfy :
    * \f$ \int_{\Omega_+}D(\omega_h) cos(\theta_h) d\omega_h = 1 \f$
    *
    * @param wh        The microfacet normal (half-vector)
    * @return scal    \f$  D(\omega_h) \f$
    */
    scal D(const vec3sc& wh) const override;

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
    scal pdf(const vec3sc& w, float u1 = 0, float u2 = 0) const override;

    scal gk_pos(const vec3sc& wo) const;
    scal gk_neg(const vec3sc& wo) const;
    scal G1_ashikhmin(const vec3sc& w) const;
    scal G1_heitz(const vec3sc& w) const;

    /**
    * @brief Smith's shadowing-masking term for a single direction
    *
    * @param w         a direction vector
    * @param wh        the microfacet normal
    * @return scal    \f$ G(\omega,\omega_h) \f$
    */
    scal G1(const vec3sc& w, const vec3sc& wh) const override;

    /**
    * @brief Correlated Smith's shadowing-masking term
    *
    * @param wi        the incident direction
    * @param wo        the outgoing direction
    * @param wh        the microfacet normal
    * @return scal    \f$ G(\omega_i, \omega_o, \omega_h) \f$
    */
    scal GAFcorrelated(const vec3sc& wi, const vec3sc& wo, const vec3sc& wh) const override;

    /**
    * @brief Sample the microfacet normal distribution
    *
    * @param u1        an uniformly distributed sample
    * @param u2        an uniformly distributed sample
    * @return vec3sc    the sampled normal
    */
    vec3sc sample(float u1, float u2) const override;

    /**
     * @brief Sample the microfacet normal distribution
     *
     * @param u1        an uniformly distributed sample
     * @param u2        an uniformly distributed sample
     * @param _pdf       out parameter, for the pdf of the sampled normal
     * @param _D         out parameter, for the distribution of the sampled normal
     * @return vec3sc    the sampled normal
     */
    vec3sc sample(float u1, float u2, scal& _pdf, scal& _D) const override;

    inline std::string name() const override { return "Discrete"; }
    inline static size_t phiSize()   { return 400; } // 400
    inline static size_t thetaSize() { return 100; } // 100
    inline const std::vector<std::vector<float>>& getValues() const { return D_values; }

    inline static float phiStart()   { return -m_pi; }
    inline static float phiEnd()     { return m_pi; }
    inline static float thetaStart() { return 0; }
    inline static float thetaEnd()   { return m_pi_2; }

    inline std::vector<float>& operator[](size_t phiIndex) { return D_values[phiIndex]; }
    inline const std::vector<float>& operator[](size_t phiIndex) const { return D_values[phiIndex]; }

    friend std::ostream& operator<<(std::ostream& output, const Discrete& mesh);
};