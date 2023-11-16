#pragma once

//#include "doctest.h"
#include "NDFs/microfacet.hpp"

namespace Test {

    float numerical_integration_ndf(MicrofacetDistribution* ndf);

    float weak_white_furnace_test(MicrofacetDistribution* ndf, const float theta_o, const float phi_o);

}

void TEST_CASE_MESH(const TriangleMesh& mesh);