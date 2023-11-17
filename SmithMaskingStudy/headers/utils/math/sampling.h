#pragma once

#include "gdt/gdt.h"
#include "gdt/math/vec.h"


namespace gdt {

    inline __both__ vec3sc uniformHemisphereSampling(scal r1, scal r2) {
        scal phi = 2.0 * m_pi * r1;
        scal sinTheta = sqrt(1.0 - r2 * r2);
        scal x = cos(phi) * sinTheta;
        scal y = sin(phi) * sinTheta;
        vec3sc dir(x, y, r2);

        return normalize(dir);
    }

};