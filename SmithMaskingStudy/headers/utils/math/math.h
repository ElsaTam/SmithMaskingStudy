#pragma once

#include <vector>
#include "gdt/math/vec.h"

namespace Conversion
{
    /**
     * @brief Convert Polar coordinates to cartesian coordinates
     *
     * @param theta     The zenithal angle
     * @param phi       The azimuthal angle
     * @return vec3sc    The cartesian vector
     */
    gdt::vec3sc polar_to_colormap(scal theta, scal phi, scal value);

    /**
     * @brief Convert Polar coordinates to cartesian coordinates
     *
     * @param theta     The zenithal angle
     * @param phi       The azimuthal angle
     * @return vec3sc    The cartesian vector
     */
    gdt::vec3sc polar_to_cartesian(scal theta, scal phi, scal radius = 1);

    /**
     * @brief Convert Cartesian coordinates to Spherical coordinates
     *
     * @param v         an arbitraty vector in cartesian coordinates
     * @return vec3sc    Spherical coordinates as (theta,phi,radius)
     */
    gdt::vec3sc cartesian_to_polar(const gdt::vec3sc& v);

    /** Angle Conversion : Degrees -> Radians */
    scal degree_to_radian(scal angle);

    /**  Angle Conversion : Radians -> Degrees */
    scal radian_to_degree(scal angle);

    gdt::vec3sc XYZToRGB(const gdt::vec3sc xyz);

    gdt::vec3sc RGBToXYZ(const gdt::vec3sc rgb);
}

namespace Geometry {
    /*! rotate a vec3 along a normal */
    __both__ gdt::vec3sc rotateAlongNormal(const gdt::vec3sc& w, const gdt::vec3sc& to, const gdt::vec3sc& from = gdt::vec3sc(0, 0, 1));
}

namespace Functions
{
    // https://github.com/lakshayg/erfinv
    // Returns a floating point number y such that std::erf(y)
    // is close to x. The current implementation is quite accurate
    // when x is away from +1.0 and -1.0. As x approaches closer
    // to those values, the error in the result increases.
    long double erfinv(long double x);

    // Refine the result of erfinv by performing Newton-Raphson
    // iteration nr_iter number of times. This method works well
    // when the value of x is away from 1.0 and -1.0
    long double erfinv(long double x, int nr_iter);
}

namespace Maths
{
    float bilinearInterpolation(scal q11, scal q12, scal q21, scal q22, scal x1, scal x2, scal y1, scal y2, scal x, scal y);
}