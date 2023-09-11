#include <gtest/gtest.h>

#include "../src/matrix.hpp"

using namespace GDM;

/*
 * Testing simple matrix creation
 */
TEST(Matrix, Creation) {
    // Arrange and Act
    Matrix<> matrix({1, 2, 3, 4, 5, 6, 7, 8, 9});

    // Assert
    EXPECT_EQ(3, matrix.getRows());
    EXPECT_EQ(3, matrix.getColoumns());
    EXPECT_EQ(1, matrix(0, 0));
    EXPECT_EQ(2, matrix(0, 1));
    EXPECT_EQ(3, matrix(0, 2));
    EXPECT_EQ(4, matrix(1, 0));
    EXPECT_EQ(5, matrix(1, 1));
    EXPECT_EQ(6, matrix(1, 2));
    EXPECT_EQ(7, matrix(2, 0));
    EXPECT_EQ(8, matrix(2, 1));
    EXPECT_EQ(9, matrix(2, 2));
}

/*
 * Testing simple matrix creation from list
 */
TEST(Matrix, CreationFromList) {
    // Arrange and Act
    Matrix<> matrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    // Assert
    EXPECT_EQ(3, matrix.getRows());
    EXPECT_EQ(3, matrix.getColoumns());
    EXPECT_EQ(1, matrix(0, 0));
    EXPECT_EQ(2, matrix(0, 1));
    EXPECT_EQ(3, matrix(0, 2));
    EXPECT_EQ(4, matrix(1, 0));
    EXPECT_EQ(5, matrix(1, 1));
    EXPECT_EQ(6, matrix(1, 2));
    EXPECT_EQ(7, matrix(2, 0));
    EXPECT_EQ(8, matrix(2, 1));
    EXPECT_EQ(9, matrix(2, 2));
}

/*
 * Testing copy constructor
 */
TEST(Matrix, Copy) {
    // Arrange
    Matrix<> firstMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    // Act
    Matrix<> secondMatrix(firstMatrix);
    firstMatrix = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    // Assert
    EXPECT_EQ(3, secondMatrix.getRows());
    EXPECT_EQ(3, secondMatrix.getColoumns());
    EXPECT_EQ(1, secondMatrix(0, 0));
    EXPECT_EQ(2, secondMatrix(0, 1));
    EXPECT_EQ(3, secondMatrix(0, 2));
    EXPECT_EQ(4, secondMatrix(1, 0));
    EXPECT_EQ(5, secondMatrix(1, 1));
    EXPECT_EQ(6, secondMatrix(1, 2));
    EXPECT_EQ(7, secondMatrix(2, 0));
    EXPECT_EQ(8, secondMatrix(2, 1));
    EXPECT_EQ(9, secondMatrix(2, 2));
}

/*
 * Testing move constructor
 */
TEST(Matrix, Move) {
    // Arrange
    Matrix<> firstMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    // Act
    Matrix<> secondMatrix(std::move(firstMatrix));
    firstMatrix = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    // Assert
    EXPECT_EQ(3, secondMatrix.getRows());
    EXPECT_EQ(3, secondMatrix.getColoumns());
    EXPECT_EQ(1, secondMatrix(0, 0));
    EXPECT_EQ(2, secondMatrix(0, 1));
    EXPECT_EQ(3, secondMatrix(0, 2));
    EXPECT_EQ(4, secondMatrix(1, 0));
    EXPECT_EQ(5, secondMatrix(1, 1));
    EXPECT_EQ(6, secondMatrix(1, 2));
    EXPECT_EQ(7, secondMatrix(2, 0));
    EXPECT_EQ(8, secondMatrix(2, 1));
    EXPECT_EQ(9, secondMatrix(2, 2));
}

/*
 * Testing eqaulity operator
 */
TEST(Matrix, Equality) {
    // Arrange
    Matrix<> firstMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<> secondMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    // Act and Assert
    Matrix<> thirdMatrix = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    ASSERT_TRUE(firstMatrix == secondMatrix);
    ASSERT_FALSE(firstMatrix == thirdMatrix);
}

/*
 * Testing inequality operator
 */
TEST(Matrix, Inequality) {
    // Arrange
    Matrix<> firstMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<> secondMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<> thirdMatrix = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    // Act and Assert
    ASSERT_FALSE(firstMatrix != secondMatrix);
    ASSERT_TRUE(firstMatrix != thirdMatrix);
}

/*
 * Testing addition operator
 */
TEST(Matrix, Addition) {
    // Arrange
    Matrix<> firstMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<> secondMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    // Act
    Matrix<> fourthMatrix = firstMatrix + secondMatrix;

    // Assert
    ASSERT_EQ(2, fourthMatrix(0, 0));
    ASSERT_EQ(4, fourthMatrix(0, 1));
    ASSERT_EQ(6, fourthMatrix(0, 2));
    ASSERT_EQ(8, fourthMatrix(1, 0));
    ASSERT_EQ(10, fourthMatrix(1, 1));
    ASSERT_EQ(12, fourthMatrix(1, 2));
    ASSERT_EQ(14, fourthMatrix(2, 0));
    ASSERT_EQ(16, fourthMatrix(2, 1));
    ASSERT_EQ(18, fourthMatrix(2, 2));
}

/*
 * Testing subtraction operator
 */
TEST(Matrix, Subtraction) {
    // Arrange
    Matrix<> firstMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<> secondMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    // Act
    Matrix<> fourthMatrix = firstMatrix - secondMatrix;

    // Assert
    ASSERT_EQ(0, fourthMatrix(0, 0));
    ASSERT_EQ(0, fourthMatrix(0, 1));
    ASSERT_EQ(0, fourthMatrix(0, 2));
    ASSERT_EQ(0, fourthMatrix(1, 0));
    ASSERT_EQ(0, fourthMatrix(1, 1));
    ASSERT_EQ(0, fourthMatrix(1, 2));
    ASSERT_EQ(0, fourthMatrix(2, 0));
    ASSERT_EQ(0, fourthMatrix(2, 1));
    ASSERT_EQ(0, fourthMatrix(2, 2));
}

/*
 * Testing multiplication operator with a scalar
 */
TEST(Matrix, ScalarMultiplication) {
    // Arrange
    Matrix<> firstMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<> secondMatrix = {2, 4, 6, 8, 10, 12, 14, 16, 18};

    // Act
    Matrix<> fourthMatrix = firstMatrix * 2;

    // Assert
    ASSERT_EQ(2, fourthMatrix(0, 0));
    ASSERT_EQ(4, fourthMatrix(0, 1));
    ASSERT_EQ(6, fourthMatrix(0, 2));
    ASSERT_EQ(8, fourthMatrix(1, 0));
    ASSERT_EQ(10, fourthMatrix(1, 1));
    ASSERT_EQ(12, fourthMatrix(1, 2));
    ASSERT_EQ(14, fourthMatrix(2, 0));
    ASSERT_EQ(16, fourthMatrix(2, 1));
    ASSERT_EQ(18, fourthMatrix(2, 2));
}

TEST(Matrix, VectorMultiplication) {
    // Arrange
    Matrix<> matrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Vector3D<> vec = {1, 2, 3};

    // Act
    Vector3D<> resultVector = matrix * vec;

    // Assert
    ASSERT_EQ(14, resultVector[0]);
    ASSERT_EQ(32, resultVector[1]);
    ASSERT_EQ(50, resultVector[2]);
}

/*
 * Testing multiplication operator with a matrix
 */
TEST(Matrix, MatrixMultiplication) {
    // Arrange
    Matrix<float, 2, 3> firstMatrix = {1, 3, -2, 0, -1, 4};
    Matrix<float, 3, 2> secondMatrix = {2, -2, 1, 5, -3, 4};

    // Act
    Matrix<float, 2, 2> fourthMatrix = firstMatrix * secondMatrix;

    // Assert
    EXPECT_EQ(11, fourthMatrix(0, 0));
    EXPECT_EQ(5, fourthMatrix(0, 1));
    EXPECT_EQ(-13, fourthMatrix(1, 0));
    EXPECT_EQ(11, fourthMatrix(1, 1));
}
/*
 * Testing division operator with a scalar
 */
TEST(Matrix, ScalarDivision) {
    // Arrange
    Matrix<> firstMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<> secondMatrix = {2, 4, 6, 8, 10, 12, 14, 16, 18};

    // Act
    Matrix<> fourthMatrix = secondMatrix / 2;

    // Assert
    ASSERT_EQ(1, fourthMatrix(0, 0));
    ASSERT_EQ(2, fourthMatrix(0, 1));
    ASSERT_EQ(3, fourthMatrix(0, 2));
    ASSERT_EQ(4, fourthMatrix(1, 0));
    ASSERT_EQ(5, fourthMatrix(1, 1));
    ASSERT_EQ(6, fourthMatrix(1, 2));
    ASSERT_EQ(7, fourthMatrix(2, 0));
    ASSERT_EQ(8, fourthMatrix(2, 1));
    ASSERT_EQ(9, fourthMatrix(2, 2));
}

/*
 * Testing determinant calculation
 */
TEST(Matrix, Determinant) {
    // Arrange
    Matrix<> firstMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<> secondMatrix = {1.0f, 2.0f, -1.0f, 2.0f, 1.0f, 2.0f, -1.0f, 2.0f, 1.0f};
    Matrix<> thirdMatrix = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    // Act and Assert
    ASSERT_EQ(0, firstMatrix.determinant());
    /*det A = 1 (cofactor of 1) + 2 (cofactor of 2) + (-1) cofactor of (-1)
     * = 1(-3) + 2(-4) + (-1)5
     * = -3 - 8 - 5
     * = -16
     */
    ASSERT_EQ(-16.0f, secondMatrix.determinant());
    ASSERT_EQ(0, thirdMatrix.determinant());
}

/*
 * Testing matrix inversion
 */
TEST(Matrix, Inversion) {
    // Arrange
    Matrix<> matrix({1.0f, 2.0f, -1.0f, 2.0f, 1.0f, 2.0f, -1.0f, 2.0f, 1.0f});

    // Act
    Matrix<> matrix2 = matrix.inverse();

    // Assert
    EXPECT_EQ(matrix2[0][0], 0.1875f);  // 3 / 16
    EXPECT_EQ(matrix2[0][1], 0.25f);    // 1 / 4
    EXPECT_EQ(matrix2[0][2], -0.3125f); // -5 / 16
    EXPECT_EQ(matrix2[1][0], 0.25f);    // 1 / 4
    EXPECT_EQ(matrix2[1][1], 0.0f);     // 0
    EXPECT_EQ(matrix2[1][2], 0.25f);    // 1 / 4
    EXPECT_EQ(matrix2[2][0], -0.3125f); // -5 / 16
    EXPECT_EQ(matrix2[2][1], 0.25f);    // 1 / 4
    EXPECT_EQ(matrix2[2][2], 0.1875f);  // 3 / 16
}

/*
 * Testing matrix inversion with inversion operator
 */
TEST(Matrix, InversionOperator) {
    // Arrange
    Matrix<> matrix({1.0f, 2.0f, -1.0f, 2.0f, 1.0f, 2.0f, -1.0f, 2.0f, 1.0f});

    // Act
    Matrix<> matrix2 = !matrix;

    // Assert
    EXPECT_EQ(matrix2[0][0], 0.1875f);  // 3 / 16
    EXPECT_EQ(matrix2[0][1], 0.25f);    // 1 / 4
    EXPECT_EQ(matrix2[0][2], -0.3125f); // -5 / 16
    EXPECT_EQ(matrix2[1][0], 0.25f);    // 1 / 4
    EXPECT_EQ(matrix2[1][1], 0.0f);     // 0
    EXPECT_EQ(matrix2[1][2], 0.25f);    // 1 / 4
    EXPECT_EQ(matrix2[2][0], -0.3125f); // -5 / 16
    EXPECT_EQ(matrix2[2][1], 0.25f);    // 1 / 4
    EXPECT_EQ(matrix2[2][2], 0.1875f);  // 3 / 16
}

TEST(Matrix, CanBePrinted) {
    // Arrange
    Matrix<> matrix({1.0f, 2.0f, -1.0f, 2.0f, 1.0f, 2.0f, -1.0f, 2.0f, 1.0f});

    // Act
    std::stringstream ss;
    ss << matrix;

    // Assert
    EXPECT_EQ(ss.str(), "(1, 2, -1)\n(2, 1, 2)\n(-1, 2, 1)\n");
}
