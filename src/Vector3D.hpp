#pragma once

#include "precompiled_headers.hpp"

namespace GDM {

template <typename T = float>
    requires std::is_floating_point_v<T>
class Vector3D {
  public:
    Vector3D();
    Vector3D(T x, T y, T z);
    Vector3D(std::initializer_list<T> list);
    Vector3D(Vector3D const & other);
    Vector3D(Vector3D && other);
    ~Vector3D();

    Vector3D & operator=(Vector3D const & other);
    Vector3D & operator=(Vector3D && other) noexcept;
    friend std::ostream & operator<<(std::ostream & os, Vector3D const & vec) {
        os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
        return os;
    }

    Vector3D operator+(Vector3D const & other) const;
    Vector3D operator-(Vector3D const & other) const;
    Vector3D operator*(T scalar) const;
    Vector3D operator/(T scalar) const;

    Vector3D & operator+=(Vector3D const & other);
    Vector3D & operator-=(Vector3D const & other);
    Vector3D & operator*=(T scalar);
    Vector3D & operator/=(T scalar);

    bool operator==(Vector3D const & other) const;
    bool operator!=(Vector3D const & other) const;

    auto operator[](int index) const -> T;

    auto X() const -> T;
    auto Y() const -> T;
    auto Z() const -> T;
    auto length() const -> T;
    auto lengthSquared() const -> T;
    auto dot(Vector3D const & other) const -> T;
    auto cross(Vector3D const & other) const -> Vector3D;
    auto normalized() const -> Vector3D;
    auto normalize() -> void;

  private:
    /// The x-coordinate of the Vector
    T x;
    /// The y-coordinate of the Vector
    T y;
    /// The z-coordinate of the Vector
    T z;
};

/// Default constructor, initializes the Vector to (0, 0, 0)
/// @tparam T Floating point type
template <typename T>
    requires std::is_floating_point_v<T>
Vector3D<T>::Vector3D() : x(0), y(0), z(0) {
}

/// Constructor, initializes the Vector to (x, y, z)
/// @tparam T Floating point type
/// @param x X coordinate
/// @param y Y coordinate
/// @param z Z coordinate
template <typename T>
    requires std::is_floating_point_v<T>
Vector3D<T>::Vector3D(T x, T y, T z) : x(x), y(y), z(z) {
}

/// Constructor, initializes the Vector to (x, y, z)
/// @tparam T Floating point type
/// @param list Initializer list
template <typename T>
    requires std::is_floating_point_v<T>
Vector3D<T>::Vector3D(std::initializer_list<T> list) {
    assert(list.size() == 3);
    auto it = list.begin();
    x = *it;
    y = *(++it);
    z = *(++it);
}

/// Copy constructor
/// @tparam T Floating point type
/// @param other Vector to copy
template <typename T>
    requires std::is_floating_point_v<T>
Vector3D<T>::Vector3D(Vector3D const & other) : x(other.x), y(other.y), z(other.z) {
}

/// Move constructor
/// @tparam T The type used to store the coordinates
/// @param other The Vector to move
template <typename T>
    requires std::is_floating_point_v<T>
Vector3D<T>::Vector3D(Vector3D && other) : x(other.x), y(other.y), z(other.z) {
}

/// Destructor
/// @tparam T Floating point type
template <typename T>
    requires std::is_floating_point_v<T>
Vector3D<T>::~Vector3D() {
}

/// Copy assignment operator
/// @tparam T Floating point type
/// @param other Vector to copy
/// @return Reference to this Vector
template <typename T>
    requires std::is_floating_point_v<T>
auto Vector3D<T>::operator=(Vector3D const & other) -> Vector3D<T> & {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
}

/// Move assignment operator
/// @tparam T Floating point type
/// @param other Vector to move
/// @return Reference to this Vector
template <typename T>
    requires std::is_floating_point_v<T>
auto Vector3D<T>::operator=(Vector3D && other) noexcept -> Vector3D<T> & {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
}

/// Addition operator
/// @tparam T Floating point type
/// @param other Vector to add
/// @return Vector sum
template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::operator+(Vector3D const & other) const -> Vector3D<T> {
    return Vector3D(x + other.x, y + other.y, z + other.z);
}

/// Subtraction operator
/// @tparam T Floating point type
/// @param other Vector to subtract
/// @return Vector difference
template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::operator-(Vector3D const & other) const -> Vector3D<T> {
    return Vector3D(x - other.x, y - other.y, z - other.z);
}

/// Multiplication operator
/// @tparam T Floating point type
/// @param scalar Scalar to multiply by
template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::operator*(T scalar) const -> Vector3D<T> {
    return Vector3D(x * scalar, y * scalar, z * scalar);
}

template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::operator/(T scalar) const -> Vector3D<T> {
    return Vector3D(x / scalar, y / scalar, z / scalar);
}

/// Addition assignment operator
/// @tparam T Floating point type
/// @param other Vector to add
/// @return Reference to this Vector
template <typename T>
    requires std::is_floating_point_v<T>
auto Vector3D<T>::operator+=(Vector3D const & other) -> Vector3D<T> & {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

/// Subtraction assignment operator
/// @tparam T The type used to store the coordinates
/// @param other Vector to subtract
/// @return Reference to this Vector
template <typename T>
    requires std::is_floating_point_v<T>
auto Vector3D<T>::operator-=(Vector3D const & other) -> Vector3D<T> & {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

/// Multiplication assignment operator
/// @tparam The type used to store the coordinates
/// @param scalar Scalar to multiply by
/// @return A reference to this Vector
template <typename T>
    requires std::is_floating_point_v<T>
auto Vector3D<T>::operator*=(T scalar) -> Vector3D<T> & {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}

/// Division assignment operator
/// @tparam The type used to store the coordinates
/// @param scalar Scalar to divide by
/// @return A reference to this Vector
template <typename T>
    requires std::is_floating_point_v<T>
auto Vector3D<T>::operator/=(T scalar) -> Vector3D<T> & {
    x /= scalar;
    y /= scalar;
    z /= scalar;
    return *this;
}

/// Equality operator
/// @tparam T The type used to store the coordinates
/// @param other Vector to compare to
/// @return True if the Vectors are equal, false otherwise
template <typename T>
    requires std::is_floating_point_v<T>
auto Vector3D<T>::operator==(Vector3D const & other) const -> bool {
    return x == other.x && y == other.y && z == other.z;
}

/// Inequality operator
/// @tparam T The type used to store the coordinates
/// @param other Vector to compare to
/// @return True if the Vectors are not equal, false otherwise
template <typename T>
    requires std::is_floating_point_v<T>
auto Vector3D<T>::operator!=(Vector3D const & other) const -> bool {
    return !(*this == other);
}

/// Index operator
/// @tparam T The type used to store the coordinates
/// @param index Index of the coordinate to return
/// @return The coordinate at the given index
template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::operator[](int index) const -> T {
    assert(index >= 0 && index < 3);
    switch (index) {
    case 0:
        return x;
    case 1:
        return y;
    default:
        return z;
    }
}

/// Returns the X coordinate
/// @tparam T The type used to store the coordinates
/// @return The X coordinate
template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::X() const -> T {
    return x;
}

/// Returns the Y coordinate
/// @tparam T The type used to store the coordinates
/// @return The Y coordinate
template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::Y() const -> T {
    return y;
}

/// Returns the Z coordinate
/// @tparam T The type used to store the coordinates
/// @return The Z coordinate
template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::Z() const -> T {
    return z;
}

/// Returns the length of the Vector
/// @tparam T The type used to store the coordinates
/// @return The length of the Vector
template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::length() const -> T {
    return std::sqrt(x * x + y * y + z * z);
}

/// Returns the squared length of the Vector
/// @tparam T The type used to store the coordinates
/// @return The squared length of the Vector
template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::lengthSquared() const -> T {
    return x * x + y * y + z * z;
}

/// Returns the dot product of this Vector and another Vector
/// @tparam T The type used to store the coordinates
/// @param other The other Vector
/// @return The dot product of this Vector and another Vector
template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::dot(Vector3D const & other) const -> T {
    return x * other.x + y * other.y + z * other.z;
}

/// Returns the cross product of this Vector and another Vector
/// @tparam T The type used to store the coordinates
/// @param other The other Vector
/// @return The cross product of this Vector and another Vector
template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::cross(Vector3D const & other) const -> Vector3D<T> {
    return Vector3D(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
}

/// Returns a normalized copy of this Vector
/// @tparam T The type used to store the coordinates
/// @return A normalized copy of this Vector
template <typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] auto Vector3D<T>::normalized() const -> Vector3D<T> {
    float len = length();
    return Vector3D(x / len, y / len, z / len);
}

/// Normalizes this Vector
/// @tparam T The type used to store the coordinates
/// @return A reference to this Vector
template <typename T>
    requires std::is_floating_point_v<T>
auto Vector3D<T>::normalize() -> void {
    T len = length();
    x /= len;
    y /= len;
    z /= len;
}

} // namespace GDM
