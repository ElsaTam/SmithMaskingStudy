#include "tests/tests_distribution.h"
#include "utils/math/math.h"
#include "utils/console.h"

float Test::numerical_integration_ndf(MicrofacetDistribution* ndf)
{
    float integral = 0.f;
    const float d_theta = 0.001f;
    const float d_phi = 0.001f;
    for (float theta = 0.f; theta <= m_pi_2; theta += d_theta)
    {
        for (float phi = 0.f; phi <= m_2_pi; phi += d_phi)
        {
            vec3sc omega_h = Conversion::polar_to_cartesian(theta, phi);
            integral += ndf->D(omega_h) * cos(theta) * sin(theta);
        }
    }
    integral *= (d_theta * d_phi);
    Console::out << Console::timePad << "\\int[D(w) cos(theta) dw] = " << integral << std::endl;
    return(integral);
}

float Test::weak_white_furnace_test(MicrofacetDistribution* ndf, const float theta_o, const float phi_o)
{
    Console::light << Console::timePad << "weak_white_furnace_test(" << theta_o << ", " << phi_o << ")" << std::endl;
    vec3sc omega_o = Conversion::polar_to_cartesian(theta_o, phi_o);

    float integral = 0.f;
    const float d_theta = 0.1f;
    const float d_phi = 0.1f;
    for (float theta = 0.f; theta <= m_pi; theta += d_theta)
    {
        Console::light << Console::timePad << "-- theta = " << theta << std::endl;
        for (float phi = 0.f; phi <= 2.f * m_pi; phi += d_phi)
        {
            vec3sc omega_i = Conversion::polar_to_cartesian(theta, phi);
            vec3sc omega_h(omega_o + omega_i);
            normalize(omega_h);

            float   d_omega_i = sinf(theta);
            if (omega_h.z <= 0.f)
                continue;

            float D = ndf->D(omega_h);
            float G1 = ndf->G1(omega_o, omega_h);
            integral += ((D * G1) / fabs(4.f * omega_o.z)) * d_omega_i;
        }
    }
    integral *= d_theta * d_phi;
    Console::out << Console::timePad << "\\int[ (D(wh) * G1((" << theta_o << ", " << phi_o << "), wh)) / | 4 cos(" << theta_o << ") | dw] = " << integral << std::endl;
    return(integral);
}




bool equal_sc(float a, float b, float epsilon)
{
    return fabs(a - b) < epsilon;
}

bool equal_sc(double a, double b, double epsilon)
{
    return abs(a - b) < epsilon;
}

void TEST_CASE_MESH(const TriangleMesh& mesh)
{
    Discrete ndf(mesh);

    Console::out << Console::timeStamp << "Testing Discrete Distribution" << std::endl;
    equal_sc(Test::numerical_integration_ndf(&ndf), 1.f, 0.01f);
    equal_sc(Test::weak_white_furnace_test(&ndf, m_pi / 7.f, m_pi / 3.f), 1.f, 0.01f);
    equal_sc(Test::weak_white_furnace_test(&ndf, m_pi / 7.f, m_pi / 4.f), 1.f, 0.01f);
    equal_sc(Test::weak_white_furnace_test(&ndf, m_pi / 7.f, m_pi / 5.f), 1.f, 0.01f);
    equal_sc(Test::weak_white_furnace_test(&ndf, m_pi / 7.f, m_pi / 6.f), 1.f, 0.01f);
    equal_sc(Test::weak_white_furnace_test(&ndf, m_pi / 7.f, m_pi / 7.f), 1.f, 0.01f);

    equal_sc(Test::weak_white_furnace_test(&ndf, m_pi / 3.f, m_pi / 7.f), 1.f, 0.01f);
    equal_sc(Test::weak_white_furnace_test(&ndf, m_pi / 4.f, m_pi / 7.f), 1.f, 0.01f);
    equal_sc(Test::weak_white_furnace_test(&ndf, m_pi / 5.f, m_pi / 7.f), 1.f, 0.01f);
    equal_sc(Test::weak_white_furnace_test(&ndf, m_pi / 6.f, m_pi / 7.f), 1.f, 0.01f);
    equal_sc(Test::weak_white_furnace_test(&ndf, m_pi / 7.f, m_pi / 7.f), 1.f, 0.01f);
}