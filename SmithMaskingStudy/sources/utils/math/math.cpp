#include "utils/math/math.h"

#include <cmath>
#include <limits>
#include <Eigen/Core>
#include <Eigen/Dense>

gdt::vec3sc Conversion::polar_to_colormap(scal theta, scal phi, scal value)
{
    return gdt::vec3sc(std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), value);
}

gdt::vec3sc Conversion::polar_to_cartesian(scal theta, scal phi, scal radius)
{
    return gdt::vec3sc(radius * std::cos(phi) * std::sin(theta), radius * std::sin(phi) * std::sin(theta), radius * std::cos(theta));
}

gdt::vec3sc Conversion::cartesian_to_polar(const gdt::vec3sc& v)
{
    scal radius = length(v);
    scal theta = std::acos(v.z / radius);
    scal phi = std::atan2(v.y, v.x);
    return gdt::vec3sc(theta, phi, radius);
}

scal Conversion::degree_to_radian(scal angle)
{
    return angle * (m_pi / 180.f);
}

scal Conversion::radian_to_degree(scal angle)
{
    return angle * (180.f / m_pi);
}

gdt::vec3sc Conversion::XYZToRGB(const gdt::vec3sc xyz)
{
    return gdt::vec3sc(
        3.240479f * xyz.x - 1.537150f * xyz.y - 0.498535f * xyz.z,
        -0.969256f * xyz.x + 1.875991f * xyz.y + 0.041556f * xyz.z,
        0.055648f * xyz.x - 0.204043f * xyz.y + 1.057311f * xyz.z
    );
}

gdt::vec3sc Conversion::RGBToXYZ(const gdt::vec3sc rgb)
{
    return gdt::vec3sc(
        0.412453f * rgb.x + 0.357580f * rgb.y + 0.180423f * rgb.z,
        0.212671f * rgb.x + 0.715160f * rgb.y + 0.072169f * rgb.z,
        0.019334f * rgb.x + 0.119193f * rgb.y + 0.950227f * rgb.z
    );
}


__both__ gdt::vec3sc Geometry::rotateAlongNormal(const gdt::vec3sc& w, const gdt::vec3sc& to, const gdt::vec3sc& from) {
    if (from == to) return w;

    Eigen::Matrix<scal, 3, 1> originalDir;
    originalDir << w.x, w.y, w.z;

    gdt::vec3sc v = cross(from, to);

    scal s = length(v);
    scal c = dot(from, to);

    Eigen::Matrix<scal, 3, 3> V;
    V << 0, -v.z, v.y,
        v.z, 0, -v.x,
        -v.y, v.x, 0;

    Eigen::Matrix<scal, 3, 3> I = Eigen::Matrix<scal, 3, 3>::Identity();

    Eigen::Matrix<scal, 3, 3> R;
    R = I + V + V * V * (1. / (1. + c));

    Eigen::Matrix<scal, 3, 1> orientedDir;
    orientedDir = R * originalDir;

    return { orientedDir(0), orientedDir(1), orientedDir(2) };
}

long double Functions::erfinv(long double x)
{
    if (x < -1 || x > 1) {
        return std::numeric_limits<long double>::quiet_NaN();
    } else if (x == 1.0) {
        return std::numeric_limits<long double>::infinity();
    } else if (x == -1.0) {
        return -std::numeric_limits<long double>::infinity();
    }

    const long double LN2 = 6.931471805599453094172321214581e-1L;

    const long double A0 = 1.1975323115670912564578e0L;
    const long double A1 = 4.7072688112383978012285e1L;
    const long double A2 = 6.9706266534389598238465e2L;
    const long double A3 = 4.8548868893843886794648e3L;
    const long double A4 = 1.6235862515167575384252e4L;
    const long double A5 = 2.3782041382114385731252e4L;
    const long double A6 = 1.1819493347062294404278e4L;
    const long double A7 = 8.8709406962545514830200e2L;

    const long double B0 = 1.0000000000000000000e0L;
    const long double B1 = 4.2313330701600911252e1L;
    const long double B2 = 6.8718700749205790830e2L;
    const long double B3 = 5.3941960214247511077e3L;
    const long double B4 = 2.1213794301586595867e4L;
    const long double B5 = 3.9307895800092710610e4L;
    const long double B6 = 2.8729085735721942674e4L;
    const long double B7 = 5.2264952788528545610e3L;

    const long double C0 = 1.42343711074968357734e0L;
    const long double C1 = 4.63033784615654529590e0L;
    const long double C2 = 5.76949722146069140550e0L;
    const long double C3 = 3.64784832476320460504e0L;
    const long double C4 = 1.27045825245236838258e0L;
    const long double C5 = 2.41780725177450611770e-1L;
    const long double C6 = 2.27238449892691845833e-2L;
    const long double C7 = 7.74545014278341407640e-4L;

    const long double D0 = 1.4142135623730950488016887e0L;
    const long double D1 = 2.9036514445419946173133295e0L;
    const long double D2 = 2.3707661626024532365971225e0L;
    const long double D3 = 9.7547832001787427186894837e-1L;
    const long double D4 = 2.0945065210512749128288442e-1L;
    const long double D5 = 2.1494160384252876777097297e-2L;
    const long double D6 = 7.7441459065157709165577218e-4L;
    const long double D7 = 1.4859850019840355905497876e-9L;

    const long double E0 = 6.65790464350110377720e0L;
    const long double E1 = 5.46378491116411436990e0L;
    const long double E2 = 1.78482653991729133580e0L;
    const long double E3 = 2.96560571828504891230e-1L;
    const long double E4 = 2.65321895265761230930e-2L;
    const long double E5 = 1.24266094738807843860e-3L;
    const long double E6 = 2.71155556874348757815e-5L;
    const long double E7 = 2.01033439929228813265e-7L;

    const long double F0 = 1.414213562373095048801689e0L;
    const long double F1 = 8.482908416595164588112026e-1L;
    const long double F2 = 1.936480946950659106176712e-1L;
    const long double F3 = 2.103693768272068968719679e-2L;
    const long double F4 = 1.112800997078859844711555e-3L;
    const long double F5 = 2.611088405080593625138020e-5L;
    const long double F6 = 2.010321207683943062279931e-7L;
    const long double F7 = 2.891024605872965461538222e-15L;

    long double abs_x = std::abs(x);

    if (abs_x <= 0.85L) {
        long double r =  0.180625L - 0.25L * x * x;
        long double num = (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0);
        long double den = (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + B0);
        return x * num / den; 
    }

    long double r = sqrt(LN2 - log(1.0L - abs_x));

    long double num, den;
    if (r <= 5.0L) {
        r = r - 1.6L;
        num = (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0);
        den = (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + D0);
    } else {
        r = r - 5.0L;
        num = (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0);
        den = (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + F0);
    }

    if (x < 0L) {
        return -num / den;
    } else {
        return num / den;
    }
}

long double Functions::erfinv(long double x, int nr_iter)
{
    long double k = 0.8862269254527580136490837416706L; // 0.5 * sqrt(pi)
    long double y = erfinv(x);
    while (nr_iter-- > 0) {
        y -= k * (std::erf(y) - x) / std::exp(-y * y);
    }
    return y;
}



/* https://en.wikipedia.org/wiki/Bilinear_interpolation
*
*    |    |          |                  |
* y2 |-- q12 ------- O --------------- q22 --
*    |    |          |                  |
*  y |----+--------- P -----------------+----
*    |    |          |                  |
*    |    |          |                  |
*    |    |          |                  |
* y1 |-- q11 ------- O --------------- q21 --
*    |    |          |                  |
*    +---------------------------------------
*        x1          x                 x2
*/
scal Maths::bilinearInterpolation(scal q11, scal q12, scal q21, scal q22, scal x1, scal x2, scal y1, scal y2, scal x, scal y)
{
    float x2x1, y2y1, x2x, y2y, yy1, xx1;
    x2x1 = x2 - x1;
    if (x2x1 == 0) { // same x, interpolate only on y
        float w1 = (y - y1) / (y2 - y1);
        return w1 * q11 + (1 - w1) * q12;
    }
    y2y1 = y2 - y1;
    if (y2y1 == 0) { // same y, interpolate only on x
        float w1 = (x - x1) / (x2 - x1);
        return w1 * q11 + (1 - w1) * q21;
    }
    x2x = x2 - x;
    y2y = y2 - y;
    yy1 = y - y1;
    xx1 = x - x1;
    float d = x2x1 * y2y1;
    float w11 = x2x * y2y / d;
    float w12 = x2x * yy1 / d;
    float w21 = xx1 * y2y / d;
    float w22 = xx1 * yy1 / d;
    return w11 * q11 + w12 * q12 + w21 * q21 + w22 * q22;
}
