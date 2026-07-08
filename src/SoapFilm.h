#pragma once
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

// MSVC does not define M_PI unless _USE_MATH_DEFINES is set before <cmath>.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Minimal 3D vector — deliberately no external math library.
struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    double norm2() const { return dot(*this); }
    double norm() const;
};

// One connected soap film patch: a triangle mesh whose boundary vertices are
// pinned to the supporting wire rings. Relaxation is monotone gradient
// descent on total surface area (surface tension acting alone), so the film
// always descends towards a discrete minimal surface and can never blow up.
//
// The axis of symmetry of the initial configuration is x; the rings live in
// planes x = const.
class FilmPatch {
public:
    std::vector<Vec3> X;                    // vertex positions
    std::vector<std::array<int, 3>> tris;   // triangle indices, CCW outward
    std::vector<uint8_t> pinned;            // 1 = boundary vertex, fixed

    // Build a tube (cylinder topology) between rings at x0 and x1 with radius
    // profile r(i) supplied by the caller. rows = axial segment count,
    // cols = angular vertex count.
    static FilmPatch tube(int rows, int cols, double x0, double x1,
                          const std::vector<double>& rowRadius);

    // One relaxation step towards the minimal surface: exact area gradient
    // (cross-product form, bounded even on degenerate triangles), diagonal
    // preconditioner, and a backtracking safeguard that reverts any step
    // that would increase area. Returns max vertex displacement.
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
    // Returns max vertex displacement.
    double implicitStep(double dt);

    // Forget cached line-search state (call after externally moving X).
    void resetDescent() { lastArea_ = -1.0; omega_ = omegaBase_; }

    double totalArea() const;

    // Mean distance of row vertices from the x-axis, for structured tubes
    // only (rows*cols layout). Returns per-row means.
    std::vector<double> rowRadii(int rows, int cols) const;

private:
    std::vector<Vec3> grad_;    // scratch: area gradient
    std::vector<Vec3> vnorm_;   // scratch: area-weighted vertex normals
    std::vector<double> diag_;  // scratch: Jacobi diagonal (step scaling)
    std::vector<Vec3> xPrev_;   // scratch: positions before the trial step
    // scratch for the implicit step / conjugate gradient
    std::vector<int> eI_, eJ_;
    std::vector<double> eW_, mass_;
    std::vector<Vec3> cgR_, cgP_, cgAp_, cgB_;
    double omegaBase_ = 0.9;    // nominal damping
    double omega_ = 0.9;        // current damping (adapted by backtracking)
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
    double separation() const { return h_; }
    double ringRadius() const { return R_; }

    // Advance n explicit MCF steps. Detects and performs pinch-off.
    // Returns max displacement of the last step.
    double relax(int n);

    // Relax until maxDisplacement < tol or maxSteps exhausted.
    // Returns the number of steps taken.
    int relaxUntil(double tol, int maxSteps);

    bool pinched() const { return pinched_; }
    double neckRadius() const;             // min row radius (tube state only)
    double totalArea() const;
    const std::vector<FilmPatch>& patches() const { return patches_; }

    // Neck radius below which the film is considered pinched. Physically a
    // real film ruptures at finite thickness; numerically we cut while the
    // neck triangles are still well-shaped. The measured critical h/R is
    // insensitive to this value because the collapse, once started, is
    // irreversible: no equilibrium exists beyond (h/R)_max.
    double pinchThreshold = 0.05;

private:
    void doPinch(int atRow);

    int rows_, cols_;
    double R_, h_;
    bool pinched_ = false;
    std::vector<FilmPatch> patches_;
};
