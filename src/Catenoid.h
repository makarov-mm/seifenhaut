#pragma once
#include <cmath>

// MSVC does not define M_PI unless _USE_MATH_DEFINES is set before <cmath>,
// and by the time this header is included that may already be too late.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Analytic reference solution for the catenoid soap film between two coaxial
// rings of radius R separated by distance h.
//
// The minimal surface of revolution is r(z) = c*cosh(z/c) with the boundary
// condition R = c*cosh(h/(2c)). Substituting u = h/(2c) gives
//     h/R = 2u / cosh(u),
// which has solutions only while h/R <= rhoMax = 2*u* / cosh(u*), where u* is
// the root of coth(u) = u. Below rhoMax there are two branches; the stable
// (physical) one is u in (0, u*].
namespace catenoid {

// Root of coth(u) = u, i.e. cosh(u) = u*sinh(u).
inline double uStar() {
    static const double v = [] {
        double lo = 1.0, hi = 1.5;
        for (int i = 0; i < 200; ++i) {
            double m = 0.5 * (lo + hi);
            (std::cosh(m) - m * std::sinh(m) > 0.0) ? lo = m : hi = m;
        }
        return 0.5 * (lo + hi);
    }();
    return v;
}

// Maximum h/R for which a catenoid exists (~1.32549).
inline double rhoMax() {
    const double u = uStar();
    return 2.0 * u / std::cosh(u);
}

// Stable-branch u for a given rho = h/R (rho must be < rhoMax).
inline double solveU(double rho) {
    double lo = 1e-12, hi = uStar();
    for (int i = 0; i < 200; ++i) {
        double m = 0.5 * (lo + hi);
        (2.0 * m / std::cosh(m) < rho) ? lo = m : hi = m;
    }
    return 0.5 * (lo + hi);
}

// Neck (waist) radius c for rings of radius R at separation h.
inline double neckRadius(double h, double R) {
    const double u = solveU(h / R);
    return h / (2.0 * u);
}

// Profile radius at axial coordinate z (z measured from the mid-plane).
inline double radiusAt(double z, double c) { return c * std::cosh(z / c); }

// Exact catenoid area for rings of radius R at separation h.
inline double area(double h, double R) {
    const double u = solveU(h / R);
    const double c = h / (2.0 * u);
    return M_PI * c * c * (std::sinh(2.0 * u) + 2.0 * u);
}

// Separation beyond which two flat disks have less area than the catenoid
// (the film becomes metastable): solve area(h,R) == 2*pi*R^2.
inline double rhoMetastable() {
    static const double v = [] {
        double lo = 0.5, hi = rhoMax() - 1e-9;
        for (int i = 0; i < 200; ++i) {
            double m = 0.5 * (lo + hi);
            (area(m, 1.0) < 2.0 * M_PI) ? lo = m : hi = m;
        }
        return 0.5 * (lo + hi);
    }();
    return v;
}

} // namespace catenoid
