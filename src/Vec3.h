#pragma once
#include <cmath>

// Minimal 3D vector — deliberately no external math library.
struct Vec3 {
    double x = 0, y = 0, z = 0;

    constexpr Vec3() = default;
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    constexpr Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    constexpr Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }

    [[nodiscard]] friend constexpr Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
    [[nodiscard]] friend constexpr Vec3 operator-(Vec3 a, const Vec3& b) { return a -= b; }
    [[nodiscard]] friend constexpr Vec3 operator*(Vec3 a, double s) { return a *= s; }
    [[nodiscard]] friend constexpr Vec3 operator*(double s, Vec3 a) { return a *= s; }
    [[nodiscard]] friend constexpr Vec3 operator-(const Vec3& a) { return {-a.x, -a.y, -a.z}; }

    [[nodiscard]] constexpr double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    [[nodiscard]] constexpr Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    [[nodiscard]] constexpr double norm2() const { return dot(*this); }
    [[nodiscard]] double norm() const { return std::sqrt(norm2()); }
};
