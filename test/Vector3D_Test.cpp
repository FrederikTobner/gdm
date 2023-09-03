#include <gtest/gtest.h>

#include "../src/Vector3D.hpp"

using namespace GDM;

TEST(Vector3D, Constructor) {
    // Arrange and Act
    Vector3D<> v1;
    // Assert
    EXPECT_EQ(v1.X(), 0);
    EXPECT_EQ(v1.Y(), 0);
    EXPECT_EQ(v1.Z(), 0);
}

TEST(Vector3D, ConstructorWithParameters) {
    // Arrange and Act
    Vector3D<> v2(1, 2, 3);

    // Assert
    EXPECT_EQ(v2.X(), 1);
    EXPECT_EQ(v2.Y(), 2);
    EXPECT_EQ(v2.Z(), 3);
}

TEST(Vector3D, CopyConstructor) {
    // Arrange
    Vector3D<> v1(1, 2, 3);

    // Act
    Vector3D<> v2(v1);

    // Assert
    EXPECT_EQ(v2.X(), 1);
    EXPECT_EQ(v2.Y(), 2);
    EXPECT_EQ(v2.Z(), 3);
}

TEST(Vector3D, MoveConstructor) {
    // Arrange
    Vector3D<> v1(1, 2, 3);

    // Act
    Vector3D<> v2(std::move(v1));

    // Assert
    EXPECT_EQ(v2.X(), 1);
    EXPECT_EQ(v2.Y(), 2);
    EXPECT_EQ(v2.Z(), 3);
}

TEST(Vector3D, AssignmentOperator) {
    // Arrange
    Vector3D<> v1(1, 2, 3);

    // Act
    Vector3D<> v2 = v1;

    // Assert
    EXPECT_EQ(v2.X(), 1);
    EXPECT_EQ(v2.Y(), 2);
    EXPECT_EQ(v2.Z(), 3);
}

TEST(Vector3D, MoveAssignment) {
    // Arrange
    Vector3D<> v1(1, 2, 3);

    // Act
    Vector3D<> v2 = std::move(v1);

    // Assert
    EXPECT_EQ(v2.X(), 1);
    EXPECT_EQ(v2.Y(), 2);
    EXPECT_EQ(v2.Z(), 3);
}

TEST(Vector3D, AdditionOperator) {
    // Arrange
    Vector3D<> v1(1, 2, 3);
    Vector3D<> v2(4, 5, 6);

    // Act
    Vector3D<> v3 = v1 + v2;

    // Assert
    EXPECT_EQ(v3.X(), 5);
    EXPECT_EQ(v3.Y(), 7);
    EXPECT_EQ(v3.Z(), 9);
}

TEST(Vector3D, SubtractionOperator) {
    // Arrange
    Vector3D<> v1(1, 2, 3);
    Vector3D<> v2(4, 5, 6);

    // Act
    Vector3D<> v3 = v1 - v2;

    // Assert
    EXPECT_EQ(v3.X(), -3);
    EXPECT_EQ(v3.Y(), -3);
    EXPECT_EQ(v3.Z(), -3);
}

TEST(Vector3D, ScalarMultiplication) {
    // Arrange
    Vector3D<> v1(1, 2, 3);

    // Act
    Vector3D<> v2 = v1 * 2;

    // Assert
    EXPECT_EQ(v2.X(), 2);
    EXPECT_EQ(v2.Y(), 4);
    EXPECT_EQ(v2.Z(), 6);
}

TEST(Vector3D, ScalarDivision) {
    //  Arrange
    Vector3D<> v1(2, 4, 6);

    // Act
    Vector3D<> v2 = v1 / 2;

    // Assert
    EXPECT_EQ(v2.X(), 1);
    EXPECT_EQ(v2.Y(), 2);
    EXPECT_EQ(v2.Z(), 3);
}

TEST(Vector3D, DotProduct) {
    // Arrange
    Vector3D<> v1(1, 2, 3);
    Vector3D<> v2(4, 5, 6);

    // Act
    float dotProduct = v1.dot(v2);

    // Assert
    EXPECT_EQ(dotProduct, 32);
}

TEST(Vector3D, CrossProduct) {
    // Arrange
    Vector3D<> v1(1, 2, 3);
    Vector3D<> v2(4, 5, 6);

    // Act
    Vector3D<> v3 = v1.cross(v2);

    // Assert
    EXPECT_EQ(v3.X(), -3);
    EXPECT_EQ(v3.Y(), 6);
    EXPECT_EQ(v3.Z(), -3);
}

TEST(Vector3D, Length) {
    // Arrange
    Vector3D<> v1(0, 3, 4);

    // Act
    float length = v1.length();

    // Assert
    EXPECT_EQ(length, 5.0);
}

TEST(Vector3D, Normalize) {
    // Arrange
    Vector3D<> v1(0, 3, 4);

    // Act
    v1.normalize();

    // Assert
    EXPECT_EQ(v1.X(), 0);
    EXPECT_EQ(v1.Y(), (float)0.6);
    EXPECT_EQ(v1.Z(), (float)0.8);
    EXPECT_EQ(v1.length(), 1.0);
}

TEST(Vector3D, CanBePrinted) {
    // Arrange
    Vector3D<> v1(1, 2, 3);

    // Act
    std::stringstream ss;
    ss << v1;

    // Assert
    EXPECT_EQ(ss.str(), "(1, 2, 3)");
}
