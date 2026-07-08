#include "SoapFilm.h"
#include "Catenoid.h"
#include <algorithm>
#include <cmath>

double Vec3::norm() const { return std::sqrt(norm2()); }

// ---------------------------------------------------------------- FilmPatch

FilmPatch FilmPatch::tube(int rows, int cols, double x0, double x1,
                          const std::vector<double>& rowRadius) {
    FilmPatch p;
    p.X.reserve((rows + 1) * cols);
    for (int i = 0; i <= rows; ++i) {
        const double t = double(i) / rows;
        const double x = x0 + (x1 - x0) * t;
        const double r = rowRadius[i];
        for (int j = 0; j < cols; ++j) {
            const double a = 2.0 * M_PI * j / cols;
            p.X.emplace_back(x, r * std::cos(a), r * std::sin(a));
        }
    }
    p.pinned.assign(p.X.size(), 0);
    for (int j = 0; j < cols; ++j) {
        p.pinned[j] = 1;
        p.pinned[rows * cols + j] = 1;
    }
    p.tris.reserve(rows * cols * 2);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j) {
            const int jn = (j + 1) % cols;
            const int a = i * cols + j, b = i * cols + jn;
            const int c = (i + 1) * cols + j, d = (i + 1) * cols + jn;
            p.tris.push_back({a, b, d});
            p.tris.push_back({a, d, c});
        }
    return p;
}

double FilmPatch::step() {
    const size_t nv = X.size();
    grad_.assign(nv, Vec3{});
    vnorm_.assign(nv, Vec3{});
    diag_.assign(nv, 0.0);

    // Exact gradient of total surface area. For a triangle (a,b,c) with unit
    // normal n, d(Area)/da = 0.5 * n x (c - b) — bounded by half the length
    // of the opposite edge even as the triangle degenerates, so unlike the
    // clamped-cotangent form it always remains a true descent direction.
    //
    // The diagonal preconditioner D_i (a Jacobi-style scaling assembled from
    // |cot|-clamped edge weights) only rescales the step per vertex; it does
    // not change the fixed point grad = 0, which is the discrete minimal
    // surface (coordinates harmonic w.r.t. Laplace-Beltrami).
    double areaNow = 0.0; // |(b-a)x(c-a)| is twice the triangle area
    for (const auto& t : tris) {
        const int a = t[0], b = t[1], c = t[2];
        const Vec3 ab = X[b] - X[a];
        const Vec3 ac = X[c] - X[a];
        const Vec3 bc = X[c] - X[b];
        Vec3 n = ab.cross(ac);
        const double ln = n.norm();
        areaNow += 0.5 * ln;
        vnorm_[a] += n; // |n| = 2*area, so this is an area-weighted normal
        vnorm_[b] += n;
        vnorm_[c] += n;
        if (ln > 1e-14) {
            n = n * (1.0 / ln);
            grad_[a] += n.cross(bc) * 0.5;
            grad_[b] += ac.cross(n) * 0.5; // = 0.5 * n x (a - c)
            grad_[c] += n.cross(ab) * 0.5;
        }
        auto cotw = [&](const Vec3& u, const Vec3& v) {
            const double cr = u.cross(v).norm();
            return std::clamp(u.dot(v) / std::max(cr, 1e-14), -20.0, 20.0);
        };
        const double wa = std::fabs(0.5 * cotw(X[b] - X[a], X[c] - X[a]));
        const double wb = std::fabs(0.5 * cotw(X[a] - X[b], X[c] - X[b]));
        const double wc = std::fabs(0.5 * cotw(X[a] - X[c], X[b] - X[c]));
        diag_[b] += wa; diag_[c] += wa; // angle at a weights edge (b,c)
        diag_[a] += wb; diag_[c] += wb;
        diag_[a] += wc; diag_[b] += wc;
    }

    // Backtracking safeguard, deferred by one step: areaNow is the result of
    // the PREVIOUS accepted move, measured for free during the gradient pass
    // (|n| is twice the triangle area). Surface tension can only ever lower
    // the area, so if the previous move raised it, that move is undone and
    // retried smaller. The scheme stays monotone — divergence is impossible
    // by construction — at half the per-step cost of a separate area pass.
    //
    // The rejection threshold sits above the FP noise floor of summing a few
    // thousand triangle areas (~1e-12 relative); with a 1e-14 threshold the
    // noise triggers spurious rejections, and since rejection halves omega
    // while acceptance only multiplies it back gently, the ratchet pins
    // omega to its floor and freezes the flow solid.
    if (lastArea_ >= 0.0 && areaNow > lastArea_ * (1.0 + 1e-11)) {
        X = xPrev_;
        omega_ = std::max(0.5 * omega_, 1e-6);
        return 1.0; // gradient is stale for the restored state; try again
    }
    lastArea_ = areaNow;
    omega_ = std::min(omegaBase_, omega_ * 1.5);

    xPrev_ = X;
    double maxDisp2 = 0.0;
    for (size_t i = 0; i < nv; ++i) {
        if (pinned[i]) continue;
        Vec3 d = grad_[i] * (-1.0 / std::max(diag_[i], 1e-12));
        // Project the motion onto the vertex normal. In the continuum the
        // area gradient is purely normal; the discrete gradient carries a
        // spurious tangential component that does next to nothing for the
        // area but slowly drags vertices across the surface until edges
        // collapse (observed: fan edges shrinking to 1e-13 near the disk
        // centre, whose degenerate-triangle noise then jams the
        // backtracking safeguard). Normal-only motion leaves the mesh
        // parametrisation alone, and the equilibria are unchanged: a
        // vanishing normal gradient is exactly zero mean curvature.
        const double nl = vnorm_[i].norm();
        if (nl > 1e-14) {
            const Vec3 nh = vnorm_[i] * (1.0 / nl);
            d = nh * d.dot(nh);
        }
        // Step cap: a fixed global scale only. An edge-proportional trust
        // region (a relic from before the backtracking safeguard was sound)
        // starves exactly the vertices that must move the most: the apex of
        // a retracting tip sits on ever-shrinking fan edges, its allowance
        // shrinks with them, and the needle never finishes retracting. With
        // every displacement proportional to omega, the area-monotone
        // backtracking alone is a complete safeguard: any step that would
        // misbehave is undone and retried smaller.
        const double lim = 0.05;
        const double dn = d.norm();
        if (dn > lim) d = d * (lim / dn);
        d = d * omega_;
        X[i] += d;
        maxDisp2 = std::max(maxDisp2, d.norm2());
    }
    return std::sqrt(maxDisp2);
}

double FilmPatch::implicitStep(double dt) {
    const size_t nv = X.size();

    // Assemble per-triangle edge weights (positively clamped cotangents,
    // duplicates across shared edges simply add) and lumped vertex masses.
    eI_.clear(); eJ_.clear(); eW_.clear();
    mass_.assign(nv, 1e-12);
    auto edgeW = [&](int i, int j, int k) {
        const Vec3 u = X[i] - X[k];
        const Vec3 v = X[j] - X[k];
        const double cr = u.cross(v).norm();
        const double cot = u.dot(v) / std::max(cr, 1e-14);
        eI_.push_back(i); eJ_.push_back(j);
        eW_.push_back(std::clamp(0.5 * cot, 1e-9, 10.0));
    };
    for (const auto& t : tris) {
        const int a = t[0], b = t[1], c = t[2];
        edgeW(b, c, a); edgeW(c, a, b); edgeW(a, b, c);
        const double A =
            0.5 * (X[b] - X[a]).cross(X[c] - X[a]).norm();
        mass_[a] += A / 3.0; mass_[b] += A / 3.0; mass_[c] += A / 3.0;
    }

    // Matrix-free apply of (M + dt*L) on the free vertices, with pinned
    // vertices treated as identity rows (Dirichlet data lives in the RHS).
    auto applyA = [&](const std::vector<Vec3>& in, std::vector<Vec3>& out) {
        for (size_t i = 0; i < nv; ++i)
            out[i] = pinned[i] ? in[i] : in[i] * mass_[i];
        for (size_t e = 0; e < eW_.size(); ++e) {
            const int i = eI_[e], j = eJ_[e];
            const double w = dt * eW_[e];
            const Vec3 d = (in[i] - in[j]) * w;
            if (!pinned[i]) out[i] += d;
            if (!pinned[j]) out[j] -= d;
        }
    };

    // RHS: M x for free rows, current position for pinned rows.
    cgB_.assign(nv, Vec3{});
    for (size_t i = 0; i < nv; ++i)
        cgB_[i] = pinned[i] ? X[i] : X[i] * mass_[i];

    // Conjugate gradient on the SPD system, warm-started from X.
    xPrev_ = X;
    cgR_.assign(nv, Vec3{}); cgP_.assign(nv, Vec3{}); cgAp_.assign(nv, Vec3{});
    applyA(X, cgAp_);
    double rr = 0.0, bb = 0.0;
    for (size_t i = 0; i < nv; ++i) {
        cgR_[i] = cgB_[i] - cgAp_[i];
        cgP_[i] = cgR_[i];
        rr += cgR_[i].dot(cgR_[i]);
        bb += cgB_[i].dot(cgB_[i]);
    }
    const double tol2 = std::max(bb, 1e-30) * 1e-20;
    for (int it = 0; it < 200 && rr > tol2; ++it) {
        applyA(cgP_, cgAp_);
        double pAp = 0.0;
        for (size_t i = 0; i < nv; ++i) pAp += cgP_[i].dot(cgAp_[i]);
        if (pAp <= 0.0) break;
        const double alpha = rr / pAp;
        double rr2 = 0.0;
        for (size_t i = 0; i < nv; ++i) {
            X[i] += cgP_[i] * alpha;
            cgR_[i] -= cgAp_[i] * alpha;
            rr2 += cgR_[i].dot(cgR_[i]);
        }
        const double beta = rr2 / rr;
        rr = rr2;
        for (size_t i = 0; i < nv; ++i)
            cgP_[i] = cgR_[i] + cgP_[i] * beta;
    }

    resetDescent(); // explicit-descent cache is stale after an implicit move
    double maxDisp2 = 0.0;
    for (size_t i = 0; i < nv; ++i)
        maxDisp2 = std::max(maxDisp2, (X[i] - xPrev_[i]).norm2());
    return std::sqrt(maxDisp2);
}

double FilmPatch::totalArea() const {
    double A = 0.0;
    for (const auto& t : tris)
        A += 0.5 * (X[t[1]] - X[t[0]]).cross(X[t[2]] - X[t[0]]).norm();
    return A;
}

std::vector<double> FilmPatch::rowRadii(int rows, int cols) const {
    std::vector<double> r(rows + 1, 0.0);
    for (int i = 0; i <= rows; ++i) {
        double s = 0.0;
        for (int j = 0; j < cols; ++j) {
            const Vec3& v = X[i * cols + j];
            s += std::sqrt(v.y * v.y + v.z * v.z);
        }
        r[i] = s / cols;
    }
    return r;
}

// ----------------------------------------------------------------- SoapFilm

SoapFilm::SoapFilm(int rows, int cols, double ringRadius, double separation)
    : rows_(rows), cols_(cols), R_(ringRadius), h_(separation) {
    // Freshly dipped film: start from a straight cylinder (that is what the
    // wire frame drags out of the soap bath) and let the flow find the
    // catenoid on its own.
    std::vector<double> rr(rows + 1, R_);
    patches_.push_back(
        FilmPatch::tube(rows_, cols_, -h_ / 2, h_ / 2, rr));
    pinched_ = false;
}

void SoapFilm::setSeparation(double h) {
    const double s = h / h_;
    for (auto& p : patches_) {
        for (size_t i = 0; i < p.X.size(); ++i) p.X[i].x *= s;
        p.resetDescent();
    }
    h_ = h;
}

double SoapFilm::relax(int n) {
    // After the pinch the task changes character: the two caps must diffuse
    // flat, which is stiff for explicit stepping. One unconditionally
    // stable implicit step does the work of thousands of explicit ones, so
    // the work budget n is exchanged at a rough cost ratio of 50:1 (a CG
    // solve costs about fifty explicit sweeps).
    if (pinched_) {
        double d = 0.0;
        const int m = std::max(1, n / 50);
        for (int it = 0; it < m; ++it) {
            d = 0.0;
            for (auto& p : patches_)
                d = std::max(d, p.implicitStep(0.05));
        }
        return d;
    }
    double d = 0.0;
    for (int it = 0; it < n; ++it) {
        d = 0.0;
        for (auto& p : patches_) d = std::max(d, p.step());
        if (!pinched_) {
            const auto rr = patches_[0].rowRadii(rows_, cols_);
            int kmin = 1;
            for (int i = 1; i < rows_; ++i)
                if (rr[i] < rr[kmin]) kmin = i;
            if (rr[kmin] < pinchThreshold) {
                doPinch(kmin);
                d = 1.0; // topology changed; not converged
            }
        }
    }
    return d;
}

int SoapFilm::relaxUntil(double tol, int maxSteps) {
    int steps = 0;
    double lastA = totalArea();
    while (steps < maxSteps) {
        const double d = relax(1);
        ++steps;
        if (d < tol) break;
        // The backtracking line search makes the raw displacement jitter, so
        // also stop when the area itself has stagnated: surface tension has
        // no other objective, a stationary area means a stationary film.
        if (steps % 200 == 0) {
            const double A = totalArea();
            if (!pinched_ && std::fabs(lastA - A) < 1e-12 * A) break;
            if (pinched_ && std::fabs(lastA - A) < 1e-13 * A) break;
            lastA = A;
        }
    }
    return steps;
}

double SoapFilm::neckRadius() const {
    if (pinched_) return 0.0;
    const auto rr = patches_[0].rowRadii(rows_, cols_);
    return *std::min_element(rr.begin(), rr.end());
}

double SoapFilm::totalArea() const {
    double A = 0.0;
    for (const auto& p : patches_) A += p.totalArea();
    return A;
}

void SoapFilm::doPinch(int atRow) {
    // Cut the tube at the thin row and resample each half onto a clean
    // structured disk mesh that follows the current axisymmetric profile
    // (row-mean radius vs row-mean x). Visually this is the same surface,
    // but the elements are well-shaped: the badly stretched slivers that
    // accumulate near the neck during the collapse would otherwise shrink
    // the trust region to nothing and leave a needle at each disk centre
    // that never finishes retracting.
    const FilmPatch t = patches_[0];
    const auto rr = t.rowRadii(rows_, cols_);
    std::vector<double> rx(rows_ + 1, 0.0);
    for (int i = 0; i <= rows_; ++i) {
        double sx = 0.0;
        for (int j = 0; j < cols_; ++j) sx += t.X[i * cols_ + j].x;
        rx[i] = sx / cols_;
    }

    auto buildHalf = [&](bool leftSide) {
        // Profile from the wire ring towards the neck, radius forced to be
        // monotone decreasing so it can be used as an interpolation table.
        std::vector<double> pr, px;
        double rmin = 1e30;
        auto push = [&](int i) {
            rmin = std::min(rmin, rr[i]);
            pr.push_back(rmin);
            px.push_back(rx[i]);
        };
        if (leftSide)
            for (int i = 0; i <= atRow; ++i) push(i);
        else
            for (int i = rows_; i >= atRow; --i) push(i);

        auto xAt = [&](double rho) {
            if (rho >= pr.front()) return px.front();
            for (size_t i = 1; i < pr.size(); ++i)
                if (rho >= pr[i]) {
                    const double f =
                        (rho - pr[i]) / std::max(pr[i - 1] - pr[i], 1e-12);
                    return px[i] + (px[i - 1] - px[i]) * f;
                }
            return px.back();
        };

        const double rTip = pr.back();
        (void)rTip; // depth of the funnel is all that survives the snap
        // The deep funnel left by the cut is one half of a catenoid — a
        // near-minimal surface whose relaxation rate vanishes (all the area
        // gradient sits in the thin tip, where clamped-cotangent slivers
        // also inflate the preconditioner). A real film crosses this region
        // on the fast inertial timescale (Taylor-Culick retraction), which
        // the quasi-static model does not contain. So the snap-back is
        // shortcut: the half is initialised as a smooth shallow dome —
        // cosine profile, zero slope at the axis, depth a fraction of the
        // funnel depth. No thin tip means a healthy normal gradient at
        // every vertex, and the dome settles briskly into the flat
        // Goldschmidt disk.
        constexpr double kSnap = 0.35;
        const double xRing = px.front();
        const double depth = kSnap * (px.back() - xRing);
        auto xCap = [&](double rho) {
            const double t = std::min(rho / R_, 1.0);
            return xRing + depth * 0.5 * (1.0 + std::cos(M_PI * t));
        };
        // Ring radii graded geometrically so that the radial spacing matches
        // the angular spacing (2*pi*rho/cols) on every ring — triangles stay
        // near-isotropic all the way to the centre. A uniform radial grid
        // would produce 5:1 slivers near the fan, and with signed cotangent
        // weights (no discrete maximum principle) such slivers support a
        // spurious non-flat critical point of the discrete area functional:
        // the flow freezes on a shallow bump with a genuinely zero gradient.
        std::vector<double> rings{R_};
        const double q = 1.0 - 2.0 * M_PI / cols_;
        const double rhoMin = 0.02 * R_;
        while (rings.back() * q > rhoMin) rings.push_back(rings.back() * q);
        const int M = int(rings.size()) - 1;
        FilmPatch p;
        p.X.reserve((M + 1) * cols_ + 1);
        for (int jr = 0; jr <= M; ++jr) {
            const double rho = rings[jr];
            const double x = xCap(rho);
            for (int j = 0; j < cols_; ++j) {
                const double a = 2.0 * M_PI * j / cols_;
                p.X.emplace_back(x, rho * std::cos(a), rho * std::sin(a));
            }
        }
        p.pinned.assign(p.X.size(), 0);
        for (int j = 0; j < cols_; ++j) p.pinned[j] = 1; // wire ring
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < cols_; ++j) {
                const int jn = (j + 1) % cols_;
                const int a = i * cols_ + j, b = i * cols_ + jn;
                const int c = (i + 1) * cols_ + j, d = (i + 1) * cols_ + jn;
                p.tris.push_back({a, b, d});
                p.tris.push_back({a, d, c});
            }
        p.X.emplace_back(xCap(0.0), 0.0, 0.0); // centre closes the dome
        p.pinned.push_back(0);
        const int apexIdx = int(p.X.size()) - 1;
        for (int j = 0; j < cols_; ++j) {
            const int jn = (j + 1) % cols_;
            p.tris.push_back({M * cols_ + j, M * cols_ + jn, apexIdx});
        }
        return p;
    };

    FilmPatch a = buildHalf(true);
    FilmPatch b = buildHalf(false);
    patches_.clear();
    patches_.push_back(std::move(a));
    patches_.push_back(std::move(b));
    pinched_ = true;
}
