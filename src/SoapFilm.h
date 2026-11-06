#pragma once
#include "Vec3.h"

#include <array>
#include <cstdint>
#include <vector>

// One connected soap film patch: a triangle mesh whose boundary vertices are
// pinned to the supporting wire rings. Relaxation is monotone gradient
// descent on total surface area (surface tension acting alone), so the film
// always descends towards a discrete minimal surface and can never blow up.
//
// The axis of symmetry of the initial configuration is x; the rings live in
// planes x = const.
class FilmPatch {
public:
    using Triangle = std::array<int, 3>;

    std::vector<Vec3> X;            // vertex positions
    std::vector<Triangle> tris;     // triangle indices, CCW outward
    std::vector<std::uint8_t> pinned; // 1 = boundary vertex, fixed

    // Build a tube (cylinder topology) between rings at x0 and x1 with the
    // radius of row i given by rowRadius[i]. rows = axial segment count,
    // cols = angular vertex count.
    static FilmPatch tube(int rows, int cols, double x0, double x1,
                          const std::vector<double>& rowRadius);

    // One relaxation step towards the minimal surface: exact area gradient
    // (cross-product form, bounded even on degenerate triangles), diagonal
    // preconditioner, and a backtracking safeguard that reverts any step
    // that would increase area. Returns the max vertex displacement.
    double step();

    // One implicit mean-curvature-flow step (implicit fairing): solves
    // (M + dt*L) x' = M x with a matrix-free conjugate gradient, pinned
    // vertices held as Dirichlet data. Unconditionally stable, so dt can be
    // large — this is what makes the post-pinch flattening fast: the dome
    // must diffuse flat from the boundary inwards, and explicit stepping
    // crosses that diffusion time scale tens of thousands of steps at a
    // time. L uses positively clamped cotangent weights (required for SPD);
    // on the near-isotropic resampled disk the clamping bias is negligible,
    // and the flat disk is an exact fixed point for any weights.
    // Returns the max vertex displacement.
    double implicitStep(double dt);

    // Forget cached line-search state (call after externally moving X).
    void resetDescent() { lastArea_ = -1.0; omega_ = kOmegaBase; }

    [[nodiscard]] double totalArea() const;

    // Mean distance of each row's vertices from the x-axis, for structured
    // tubes (rows*cols layout) only.
    [[nodiscard]] std::vector<double> rowRadii(int rows, int cols) const;

private:
    static constexpr double kOmegaBase = 0.9; // nominal damping
    static constexpr double kStepCap = 0.05;  // max displacement per step, in R

    std::vector<Vec3> grad_;    // scratch: area gradient
    std::vector<Vec3> vnorm_;   // scratch: area-weighted vertex normals
    std::vector<double> diag_;  // scratch: Jacobi diagonal (step scaling)
    std::vector<Vec3> xPrev_;   // scratch: positions before the trial step
    // scratch for the implicit step / conjugate gradient
    std::vector<int> eI_, eJ_;
    std::vector<double> eW_, mass_;
    std::vector<Vec3> cgR_, cgP_, cgAp_, cgB_;
    double omega_ = kOmegaBase; // current damping (adapted by backtracking)
    double lastArea_ = -1.0;    // area accepted at the previous step
};

// The whole simulation: one tube that may pinch off into two capped disks.
class SoapFilm {
public:
    // rows = axial segments, cols = angular resolution.
    SoapFilm(int rows, int cols, double ringRadius, double separation);

    // Move the rings to a new separation. Interior x-coordinates are
    // stretched affinely; radii are kept, then relaxation takes over.
    void setSeparation(double h);
    [[nodiscard]] double separation() const { return h_; }
    [[nodiscard]] double ringRadius() const { return R_; }

    // Advance the relaxation by n units of work (explicit steps before the
    // pinch, implicit steps at a 50:1 exchange rate after it). Detects and
    // performs the pinch-off. Returns the max displacement of the last step.
    double relax(int n);

    // Relax until maxDisplacement < tol or maxSteps exhausted.
    // Returns the number of steps taken.
    int relaxUntil(double tol, int maxSteps);

    [[nodiscard]] bool pinched() const { return pinched_; }
    [[nodiscard]] double neckRadius() const; // min row radius (tube state only)
    [[nodiscard]] double totalArea() const;
    [[nodiscard]] const std::vector<FilmPatch>& patches() const { return patches_; }

    // Neck radius below which the film is considered pinched. Physically a
    // real film ruptures at finite thickness; numerically we cut while the
    // neck triangles are still well-shaped. The measured critical h/R is
    // insensitive to this value because the collapse, once started, is
    // irreversible: no equilibrium exists beyond (h/R)_max.
    double pinchThreshold = 0.05;

private:
    void doPinch(int atRow);
    FilmPatch buildCap(double xRing, double xNeck) const;

    int rows_, cols_;
    double R_, h_;
    bool pinched_ = false;
    std::vector<FilmPatch> patches_;
};
