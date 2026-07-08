// Headless verification of the SoapFilm solver against the analytic catenoid.
// Build:  g++ -O2 -std=c++17 -I../src verify.cpp ../src/SoapFilm.cpp -o verify
#include "../src/Catenoid.h"
#include "../src/SoapFilm.h"
#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char* what) {
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok) ++failures;
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("Analytic constants:\n");
    std::printf("  u*            = %.6f  (expected 1.199679)\n", catenoid::uStar());
    std::printf("  (h/R)_max     = %.6f  (expected 1.325487)\n", catenoid::rhoMax());
    std::printf("  (h/R)_meta    = %.6f  (expected ~1.0554)\n", catenoid::rhoMetastable());
    check(std::fabs(catenoid::uStar() - 1.199679) < 1e-5, "u* root of coth(u)=u");
    check(std::fabs(catenoid::rhoMax() - 1.325487) < 1e-5, "critical ratio");

    // 1. Equilibrium shape at h/R = 0.8 vs analytic profile ------------------
    {
        const double h = 0.8;
        SoapFilm film(48, 64, 1.0, h);
        const int steps = film.relaxUntil(1e-8, 60000);
        const double c = catenoid::neckRadius(h, 1.0);
        const auto rr = film.patches()[0].rowRadii(48, 64);
        double maxErr = 0.0;
        for (int i = 0; i <= 48; ++i) {
            const double z = -h / 2 + h * i / 48.0;
            const double ra = catenoid::radiusAt(z, c);
            maxErr = std::max(maxErr, std::fabs(rr[i] - ra) / ra);
        }
        const double A = film.totalArea();
        const double Aa = catenoid::area(h, 1.0);
        std::printf("\nEquilibrium at h/R = 0.8 (%d steps):\n", steps);
        std::printf("  neck radius   = %.5f  (analytic %.5f)\n", film.neckRadius(), c);
        std::printf("  max profile rel. error = %.2e\n", maxErr);
        std::printf("  area          = %.6f  (analytic %.6f, rel err %.2e)\n",
                    A, Aa, std::fabs(A - Aa) / Aa);
        check(maxErr < 5e-3, "profile matches c*cosh(z/c) within 0.5%");
        check(std::fabs(A - Aa) / Aa < 5e-3, "area matches analytic within 0.5%");
        check(!film.pinched(), "film intact well below the limit");
    }

    // 2. Quasi-static ramp: measure the critical separation ------------------
    {
        SoapFilm film(32, 48, 1.0, 1.20);
        film.relaxUntil(1e-8, 60000);
        double hCrit = -1.0;
        for (double h = 1.20; h < 1.45; h += 0.005) {
            film.setSeparation(h);
            film.relaxUntil(3e-8, 12000);
            if (film.pinched()) { hCrit = h; break; }
        }
        std::printf("\nQuasi-static ramp:\n");
        std::printf("  measured critical h/R = %.4f  (analytic %.4f)\n",
                    hCrit, catenoid::rhoMax());
        check(hCrit > 0, "film does collapse");
        check(std::fabs(hCrit - catenoid::rhoMax()) < 0.02,
              "critical ratio within 1.5% of 1.32549");

        // 3. Goldschmidt end state: two flat disks ----------------------------
        film.relaxUntil(1e-9, 40000);
        const double A = film.totalArea();
        std::printf("\nPost-collapse (Goldschmidt) state:\n");
        std::printf("  patches       = %zu\n", film.patches().size());
        std::printf("  total area    = %.6f  (2*pi = %.6f, rel err %.2e)\n",
                    A, 2 * M_PI, std::fabs(A - 2 * M_PI) / (2 * M_PI));
        check(film.patches().size() == 2, "two separate films after pinch");
        check(std::fabs(A - 2 * M_PI) / (2 * M_PI) < 5e-3,
              "halves flatten to disks (area -> 2*pi*R^2)");

        // Flatness: every vertex of each half must end up in its wire ring
        // plane. This is the direct test for the needle artifact — a tail
        // left over from the pinch that never finishes retracting.
        double flatErr = 0.0;
        for (const auto& p : film.patches()) {
            double x0 = 0.0;
            int npin = 0;
            for (size_t i = 0; i < p.X.size(); ++i)
                if (p.pinned[i]) { x0 += p.X[i].x; ++npin; }
            x0 /= npin;
            for (const auto& v : p.X)
                flatErr = std::max(flatErr, std::fabs(v.x - x0));
        }
        std::printf("  max deviation from ring plane = %.2e\n", flatErr);
        check(flatErr < 1e-2, "disks are flat (no residual needle)");
    }

    std::printf("\n%s (%d failure%s)\n", failures ? "VERIFICATION FAILED" : "ALL CHECKS PASSED",
                failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
