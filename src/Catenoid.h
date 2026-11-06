#pragma once
#include <cmath>
#include <numbers>

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

// Root of f on [lo, hi], where f(lo) > 0 > f(hi) (or the reverse): 200
// halvings take any bracket well past double precision.
inline double bisect(auto&& f, double lo, double hi) {
    const bool rising = f(lo) < 0.0;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        ((f(mid) < 0.0) == rising ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

// Root of coth(u) = u, i.e. cosh(u) = u*sinh(u).
inline double uStar() {
    static const double v = bisect([](double u) { return std::cosh(u) - u * std::sinh(u); }, 1.0, 1.5);
    return v;
}

// Maximum h/R for which a catenoid exists (~1.32549).
inline double rhoMax() {
    const double u = uStar();
    return 2.0 * u / std::cosh(u);
}

// Stable-branch u for a given rho = h/R (rho must be < rhoMax).
inline double solveU(double rho) {
    return bisect([rho](double u) { return 2.0 * u / std::cosh(u) - rho; }, 1e-12, uStar());
}

// Neck (waist) radius c for rings of radius R at separation h.
inline double neckRadius(double h, double R) {
    return h / (2.0 * solveU(h / R));
}

// Profile radius at axial coordinate z (z measured from the mid-plane).
inline double radiusAt(double z, double c) { return c * std::cosh(z / c); }

// Exact catenoid area for rings of radius R at separation h.
inline double area(double h, double R) {
    const double u = solveU(h / R);
    const double c = h / (2.0 * u);
    return std::numbers::pi * c * c * (std::sinh(2.0 * u) + 2.0 * u);
}

// Separation beyond which two flat disks have less area than the catenoid
// (the film becomes metastable): solve area(h,R) == 2*pi*R^2.
inline double rhoMetastable() {
    static const double v = bisect([](double rho) { return area(rho, 1.0) - 2.0 * std::numbers::pi; },
                                   0.5, rhoMax() - 1e-9);
    return v;
}

} // namespace catenoid
