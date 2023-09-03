#pragma once

#include "precompiled_headers.hpp"

#include "Vector3D.hpp"

namespace GDM {

/// @brief Represents a matrix of type T with m rows and n coloumns
/// @tparam T The type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
template <typename T = float, std::size_t m = 3, std::size_t n = 3>
  requires std::is_floating_point_v<T>
class Matrix {

public:
  // Represents a row of the matrix
  class Row;

  explicit Matrix();
  Matrix(const Matrix<T, m, n> &mat);
  explicit Matrix(const std::vector<T> values);
  Matrix(const std::initializer_list<T> list);
  Matrix(Matrix<T, m, n> &&mat);
  Matrix(std::array<std::array<T, n>, m> values);
  ~Matrix();

  friend std::ostream &operator<<(std::ostream &os,
                                  const Matrix<T, m, n> &mat) {
    for (std::size_t i = 0; i < m; i++) {
      os << "(";
      for (std::size_t j = 0; j < n; j++) {
        os << mat.m_matrix[i][j];
        if(j != n -1) {
          os << ", ";
        }
      }
      os << ")\n";
    }
    return os;
  }

  inline auto operator()(std::size_t i, std::size_t j) -> T &;
  inline auto operator()(std::size_t i, std::size_t j) const -> T &;
  auto operator[](std::size_t row) -> Row { return Row(*this, row); }
  auto operator[](std::size_t row) const -> Row & { return Row(*this, row); }
  auto operator==(const Matrix<T, m, n> &mat) const -> bool;
  auto operator!=(const Matrix<T, m, n> &mat) const -> bool;
  auto operator=(const Matrix<T, m, n> &mat) -> Matrix<T, m, n> &;
  auto operator=(Matrix<T, m, n> &&mat) -> Matrix<T, m, n> &;

  auto operator+(const Matrix<T, m, n> &mat) const -> Matrix<T, m, n>;
  auto operator-(const Matrix<T, m, n> &mat) const -> Matrix<T, m, n>;
  auto operator*(const T &scalar) const -> Matrix<T, m, n>;
  auto operator*(const Vector3D<T> &vec) const -> Vector3D<T>;
  auto operator/(const T &scalar) const -> Matrix<T, m, n>;

  auto operator+=(const Matrix<T, m, n> &mat) -> Matrix<T, m, n> &;
  auto operator-=(const Matrix<T, m, n> &mat) -> Matrix<T, m, n> &;
  auto operator*=(const T &scalar) -> Matrix<T, m, n> &;
  auto operator*=(const Matrix<T, m, n> &mat) -> Matrix<T, m, n> &;
  auto operator/=(const T &scalar) -> Matrix<T, m, n> &;

  inline auto operator!() -> Matrix<T, m, n>;

  inline constexpr auto getRows() -> std::size_t;
  inline constexpr auto getColoumns() -> std::size_t;

  auto determinant() const -> T;
  auto cofactor(std::size_t i, std::size_t j) const -> T;
  auto minor(std::size_t i, std::size_t j) const -> T;
  auto inverse() -> Matrix<T, m, n>;
  auto transpose() const -> Matrix<T, m, n>;
  auto adjoint() const -> Matrix<T, m, n>;
  auto identity() const -> Matrix<T, m, n>;
  /// Implementation of the matrix row
  class Row {
  public:
    Matrix &parent;
    std::size_t row;
    Row(Matrix &parent, std::size_t row) : parent(parent), row(row) {}
    auto operator[](std::size_t col) -> T & { return parent(row, col); }
  };

private:
  /// The values stored in the matrix
  std::array<std::array<T, n>, m> m_matrix;
};

/// @brief Construct a new Matrix object
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
Matrix<T, m, n>::Matrix() {
  static_assert(m > 0 && n > 0);
}

/// @brief Construct a new Matrix object
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param mat The matrix to copy
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
Matrix<T, m, n>::Matrix(const Matrix<T, m, n> &mat) {
  static_assert(m > 0 && n > 0);
  m_matrix = mat.m_matrix;
}

/// @brief Construct a new Matrix object
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param mat The matrix to move
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
Matrix<T, m, n>::Matrix(Matrix<T, m, n> &&mat) {
  static_assert(m > 0 && n > 0);
  m_matrix = std::move(mat.m_matrix);
}

/**
 * @brief Construct a new Matrix object
 *
 * @tparam T the type of value used to build the matrix
 * @tparam m the number of rows in the matrix
 * @tparam n the number of columns in the matrix
 * @param values the values for the matrix in row-major order
 */
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
Matrix<T, m, n>::Matrix(const std::vector<T> values) {

  // check that the matrix invariants are satisfied
  static_assert(m > 0 && n > 0);
  // check that the values vector is not empty
  if (values.empty()) {
    throw std::invalid_argument("The provided vector is empty");
  }
  // check that the number of values is equal to rowsCount * columnsCount
  if (m * n != values.size()) {
    throw std::runtime_error("Total number of matrix values does not match the "
                             "rows and coloumns provided");
  }
  // populate the matrix with the provided values
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      this->m_matrix[i][j] = values[i * m + j];
    }
  }
}

/// @brief Construct a new Matrix object from a list of values
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param list The list of values to construct the matrix from
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
Matrix<T, m, n>::Matrix(std::initializer_list<T> list) {

  static_assert(m > 0 && n > 0);
  if (list.size() == 0) {
    throw std::invalid_argument("The provided list is empty");
  }

  if (m * n != list.size()) {
    throw std::runtime_error("Total number of matrix values does not match the "
                             "rows and coloumns provided");
  }
  auto it = list.begin();
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++, std::advance(it, 1)) {
      this->m_matrix[i][j] = *it;
    }
  }
}

/// @brief Construct a new Matrix object from an array of values
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of columns in the matrix
/// @param values The array of values to construct the matrix from
/// @return Matrix<T, m, n> The constructed matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
Matrix<T, m, n>::Matrix(std::array<std::array<T, n>, m> mat) {
  static_assert(m > 0 && n > 0);
  m_matrix = mat;
}

/// Destructor for the matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of columns in the matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
Matrix<T, m, n>::~Matrix() {}

/// @brief Get a reference to the value at the specified row and coloumn
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param i The row index
/// @param j The coloumn index
/// @return T& The value at the specified row and coloumn
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] inline auto Matrix<T, m, n>::operator()(std::size_t i,
                                                      std::size_t j) -> T & {
  assert(i < m && j < n);
  return this->m_matrix[i][j];
}

/// @brief Get a const reference to the value at the specified row and coloumn
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param i The row index
/// @param j The coloumn index
/// @return T& The value at the specified row and coloumn
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] inline auto Matrix<T, m, n>::operator()(std::size_t i,
                                                      std::size_t j) const
    -> T & {
  assert(i < m && j < n);
  return this->m_matrix[i][j];
}

/// @brief Copy assignment operator
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param mat The matrix to copy
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
auto Matrix<T, m, n>::operator=(const Matrix<T, m, n> &mat)
    -> Matrix<T, m, n> & {
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      this->m_matrix[i][j] = mat[i][j];
    }
  }
  return *this;
}

/// @brief Move assignment operator
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param mat The matrix to move
/// @return Matrix<T, m, n> & A reference to the moved matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
auto Matrix<T, m, n>::operator=(Matrix<T, m, n> &&mat) -> Matrix<T, m, n> & {
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      this->m_matrix[i][j] = std::move(mat[i][j]);
    }
  }
  return *this;
}

/// Adds two matrices together
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam b The number of coloumns in the matrix
/// @param mat The matrix to add
/// @return Matrix<T, m, n> The sum of the two matrices
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] auto Matrix<T, m, n>::operator+(const Matrix<T, m, n> &mat) const
    -> Matrix<T, m, n> {
  std::array<std::array<T, n>, m> result;
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      result[i][j] = this->m_matrix[i][j] + mat.m_matrix[i][j];
    }
  }
  return Matrix<T, m, n>(result);
}

/// Subtracts a matrix from this matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam b The number of coloumns in the matrix
/// @param mat The matrix to subtract
/// @return Matrix<T, m, n> The difference of the two matrices
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] auto Matrix<T, m, n>::operator-(const Matrix<T, m, n> &mat) const
    -> Matrix<T, m, n> {
  std::array<std::array<T, n>, m> result;
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      result[i][j] = this->m_matrix[i][j] - mat.m_matrix[i][j];
    }
  }
  return Matrix<T, m, n>(result);
}

/// Multiplies the Matrix by a scalar
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param scalar The scalar to multiply by
/// @return Matrix<T, m, n> The product of the matrix and the scalar
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] auto Matrix<T, m, n>::operator*(const T &scalar) const
    -> Matrix<T, m, n> {
  std::array<std::array<T, n>, m> result;
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      result[i][j] = this->m_matrix[i][j] * scalar;
    }
  }
  return Matrix<T, m, n>(result);
}

/// Multiplies the matrix with a vector
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param vec The vector to multiply by
/// @return Vector3D<T> The product of the matrix and the vector
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] auto Matrix<T, m, n>::operator*(const Vector3D<T> &vec) const
    -> Vector3D<T> {
  static_assert(m == 3 && n == 3, "Matrix must have 3 columns and 3 rows");
  T x = 0;
  T y = 0;
  T z = 0;
  for (std::size_t i = 0; i < m; i++) {
    x += this->m_matrix[0][i] * vec[i];
    y += this->m_matrix[1][i] * vec[i];
    z += this->m_matrix[2][i] * vec[i];
  }
  return Vector3D<T>(x, y, z);
}

/// Divides the Matrix by a scalar
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param scalar The scalar to divide by
/// @return Matrix<T, m, n> The quotient of the matrix and the scalar
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] auto Matrix<T, m, n>::operator/(const T &scalar) const
    -> Matrix<T, m, n> {
  std::array<std::array<T, n>, m> result;
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      result[i][j] = this->m_matrix[i][j] / scalar;
    }
  }
  return Matrix<T, m, n>(result);
}

/// @brief Add a matrix to this matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param mat The matrix to add
/// @return Matrix<T, m, n> & A reference to the moved matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] auto Matrix<T, m, n>::operator+=(const Matrix<T, m, n> &mat)
    -> Matrix<T, m, n> & {
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      m_matrix[i][j] += mat(i, j);
    }
  }
  return *this;
}

/// Subtracts a matrix from this matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param mat The matrix to subtract
/// @return Matrix<T, m, n> & A reference to the moved matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] Matrix<T, m, n> &
Matrix<T, m, n>::operator-=(const Matrix<T, m, n> &mat) {
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      this->m_matrix[i][j] -= mat[i][j];
    }
  }
  return *this;
}

/// @brief Multiply this matrix by a scalar
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param scalar The scalar to multiply by
/// @return Matrix<T, m, n> & A reference to the moved matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] auto Matrix<T, m, n>::operator*=(const T &scalar)
    -> Matrix<T, m, n> & {
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      this->m_matrix[i][j] *= scalar;
    }
  }
  return *this;
}

/// @brief Multiply this matrix by another matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param mat The matrix to multiply by
/// @return Matrix<T, m, n> & A reference to the moved matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] auto Matrix<T, m, n>::operator*=(const Matrix<T, m, n> &mat)
    -> Matrix<T, m, n> & {
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      this->m_matrix[i][j] *= mat[i][j];
    }
  }
  return *this;
}

/// @brief Divide this matrix by a scalar
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param scalar The scalar to divide by
/// @return Matrix<T, m, n> & A reference to the moved matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] auto Matrix<T, m, n>::operator/=(const T &scalar)
    -> Matrix<T, m, n> & {
  return *this = this * scalar;
}

/// Compares two matrices for equality
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param mat The matrix to compare to
/// @return true If the matrices are equal
/// @return false If the matrices are not equal
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] inline auto
Matrix<T, m, n>::operator==(const Matrix<T, m, n> &mat) const -> bool {
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      if (this->m_matrix[i][j] != mat.m_matrix[i][j]) {
        return false;
      }
    }
  }
  return true;
}

/// Creates the inverse of the matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam b The number of coloumns in the matrix
/// @return Matrix<T, m, n> The inverse of the matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] inline auto Matrix<T, m, n>::operator!() -> Matrix<T, m, n> {
  return this->inverse();
}

/// Compares two matrices for inequality
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param mat The matrix to compare to
/// @return true If the matrices are not equal
/// @return false If the matrices are equal
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] inline auto
Matrix<T, m, n>::operator!=(const Matrix<T, m, n> &mat) const -> bool {
  return !(*this == mat);
}

/// @brief Get the number of rows in the matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @return std::size_t The number of rows in the matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] constexpr inline auto Matrix<T, m, n>::getRows() -> std::size_t {
  return m;
}

/// @brief Gets the number of coloumns in the matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @return std::size_t The number of coloumns in the matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] constexpr inline auto Matrix<T, m, n>::getColoumns()
    -> std::size_t {
  return n;
}

/// Calculates the determinant of the matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @return T The determinant of the matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] inline auto Matrix<T, m, n>::determinant() const -> T {
  static_assert(m == n, "Matrix must be square");
  if constexpr (m == 1) {
    return this->m_matrix[0][0];
  } else if constexpr (m == 2) {
    return this->m_matrix[0][0] * this->m_matrix[1][1] -
           this->m_matrix[0][1] * this->m_matrix[1][0];
  } else {
    T det = 0;
    for (std::size_t i = 0; i < n; i++) {
      det += this->m_matrix[0][i] * this->cofactor(0, i);
    }
    return det;
  }
}

/// Calculates the cofactor of the matrix at the specified row and coloumn
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @param i The row index
/// @param j The coloumn index
/// @return T The cofactor of the matrix at the specified row and coloumn
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] inline auto Matrix<T, m, n>::cofactor(std::size_t i,
                                                    std::size_t j) const -> T {
  static_assert(m == n, "Matrix must be square");
  assert(i < m && j < n);
  return std::pow(-1, i + j) * this->minor(i, j);
}

/**
 * Calculates the minor of the matrix at the specified row and column
 *
 * @tparam T The underlying type of the matrix
 * @tparam m The number of rows in the matrix
 * @tparam n The number of columns in the matrix
 * @param i The row index
 * @param j The column index
 *
 * @return T The minor of the matrix at the specified row and column
 */
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] inline auto Matrix<T, m, n>::minor(std::size_t i,
                                                 std::size_t j) const -> T {
  static_assert(m == n, "Matrix must be square");
  assert(i < m && j < n);
  Matrix<T, m - 1, n - 1> minor;
  std::size_t minorRow = 0;
  std::size_t minorCol = 0;
  for (std::size_t row = 0; row < m; row++) {
    for (std::size_t col = 0; col < n; col++) {
      if (row != i && col != j) {
        minor[minorRow][minorCol] = this->m_matrix[row][col];
        minorCol++;
        if (minorCol == minor.getColoumns()) {
          minorCol = 0;
          minorRow++;
        }
      }
    }
  }
  return minor.determinant();
}

/// @brief Calculates the inverse of the matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @return Matrix<T, m, n> The inverse of the matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] inline auto Matrix<T, m, n>::inverse() -> Matrix<T, m, n> {
  static_assert(m == n, "Matrix must be square");
  std::array<std::array<T, n>, m> inverted;
  T det = this->determinant();
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      inverted[j][i] = this->cofactor(i, j) / det;
    }
  }
  return Matrix<T, m, n>(inverted);
}

/// @brief Calculates the transpose of the matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @return Matrix<T, m, n> The transpose of the matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] inline auto Matrix<T, m, n>::transpose() const
    -> Matrix<T, m, n> {
  Matrix<T, m, n> transpose;
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      transpose[j][i] = this->m_matrix[i][j];
    }
  }
  return transpose;
}

/// @brief Calculates the adjoint of the matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @return Matrix<T, m, n> The adjoint of the matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] inline auto Matrix<T, m, n>::adjoint() const -> Matrix<T, m, n> {
  static_assert(m == n, "Matrix must be square");
  Matrix<T, m, n> adjoint;
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      adjoint[i][j] = this->cofactor(i, j);
    }
  }
  return adjoint;
}

/// @brief Calculates the identity matrix of the same size as the matrix
/// @tparam T The underlying type of the matrix
/// @tparam m The number of rows in the matrix
/// @tparam n The number of coloumns in the matrix
/// @return Matrix<T, m, n> The identity matrix of the same size as the matrix
template <typename T, std::size_t m, std::size_t n>
  requires std::is_floating_point_v<T>
[[nodiscard]] auto Matrix<T, m, n>::identity() const -> Matrix<T, m, n> {
  Matrix<T, m, n> identity(*this->m_matrix);
  for (std::size_t i = 0; i < m; i++) {
    for (std::size_t j = 0; j < n; j++) {
      identity[i][j] = (i == j) ? 1 : 0;
    }
  }
  return identity;
}

} // namespace GDM
