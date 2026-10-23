#include "3dGeomUtil.h"

#include<deque>
#include<vector>
#include<iostream>

#include "gtest/gtest.h"


constexpr double cgEpsilon = 0.0001;

std::ostream& operator<<(std::ostream &aOut, Plane const &aPlane) {
  aOut << aPlane.mNormal(0) << ' ' << aPlane.mNormal(1) << ' ' << aPlane.mNormal(2) << ' ' << " : " << aPlane.mConstant;
  return aOut;
}

bool eq(Vector const &aV1, Vector const &aV2, double const aEpsilon) {
  return (aV1 - aV2).norm() < aEpsilon;
}

bool eq(Vector const &aV1, Vector const &aV2) {
  return eq(aV1, aV2, cgEpsilon);
}

bool eq(double const aF1, double const aF2, double const aEpsilon) {
  return ::abs(aF1 - aF2) < aEpsilon;
}

bool eq(double const aF1, double const aF2) {
  return eq(aF1, aF2, cgEpsilon);
}

Vertex getPlaneIntersectionNormals(Vertex aPoint,
                                         Vector aDir1,
                                         Vector aDir2,
                                         Vector aDir3) {
  Vector normal1 = aDir1.normalized();
  Plane plane1(normal1, normal1.dot(aPoint));
  Vector normal2 = aDir2.normalized();
  Plane plane2(normal2, normal2.dot(aPoint));
  Vector normal3 = aDir3.normalized();
  Plane plane3(normal3, normal3.dot(aPoint));
  return Plane::intersect(plane1, plane2, plane3);
}

TEST(vector, getAperpendicular) {
  static constexpr double csEpsilon = 1e-10;
  {
    Vector vector{1.0, 0.0, 0.0};
    EXPECT_LT(::abs(util::getAperpendicular(vector).dot(vector)), csEpsilon);
  }
  {
    Vector vector{1.0, 1.0, 0.0};
    vector.normalize();
    EXPECT_LT(::abs(util::getAperpendicular(vector).dot(vector)), csEpsilon);
  }
  {
    Vector vector{1.0, 0.0, 1.0};
    vector.normalize();
    EXPECT_LT(::abs(util::getAperpendicular(vector).dot(vector)), csEpsilon);
  }
  {
    Vector vector{1.0, -1.0, -1.0};
    vector.normalize();
    EXPECT_LT(::abs(util::getAperpendicular(vector).dot(vector)), csEpsilon);
  }
}

TEST(ray, averageErrorSquared) {
  {
    Ray ray({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0});
    std::vector<Vertex> points{};
    EXPECT_EQ(ray.getAverageErrorSquared(points), 0.0);
  }
  {
    Ray ray({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0});
    std::vector<Vertex> points{{2.0, 0.0, 0.0}, {-3.0, 0.0, 0.0}};
    EXPECT_EQ(ray.getAverageErrorSquared(points), 0.0);
  }
  {
    Ray ray({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0});
    std::vector<Vertex> points{{2.0, 1.0, 0.0}, {-3.0, 0.0, 1.0}};
    EXPECT_GT(ray.getAverageErrorSquared(points), 0.0);
  }
}

TEST(planeIntersection, Normals) {
  {
    Vertex common{1.0, 2.0, 3.0};
    auto result = getPlaneIntersectionNormals(common, {1.0, 2.0, 3.0}, {3.0, 1.0, 2.0}, {3.0, 2.0, 1.0});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{3.0, -2.0, 1.0};
    auto result = getPlaneIntersectionNormals(common, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{3.0, -2.0, -1.0};
    auto result = getPlaneIntersectionNormals(common, {1.0, -2.0, 3.0}, {-1.0, 2.0, 3.0}, {1.0, 2.0, -3.0});
    EXPECT_TRUE(eq(common, result));
  }
}

Vertex getPlaneIntersectionProportionPoints(Vertex aCommon,
                                                  double aProportion1, Vertex aOne1,
                                                  double aProportion2, Vertex aOne2,
                                                  double aProportion3, Vertex aOne3) {
  Plane plane1 = Plane::createFrom1proportion2points(aProportion1, aOne1, aOne1 + 1.0 / aProportion1 * (aCommon - aOne1));
  Plane plane2 = Plane::createFrom1proportion2points(aProportion2, aOne2, aOne2 + 1.0 / aProportion2 * (aCommon - aOne2));
  Plane plane3 = Plane::createFrom1proportion2points(aProportion3, aOne3, aOne3 + 1.0 / aProportion3 * (aCommon - aOne3));
  return Plane::intersect(plane1, plane2, plane3);
}

TEST(planeIntersection, Proportion) {
  {
    Vertex common{0.0, 0.0, 0.0};
    auto result = getPlaneIntersectionProportionPoints(common, 0.5, {1.0, 0.0, 0.0}, 0.5, {0.0, 1.0, 0.0}, 0.5, {0.0, 0.0, 1.0});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{0.0, 0.0, 0.0};
    auto result = getPlaneIntersectionProportionPoints(common, 0.5, {1.0, 0.0, 0.0}, 0.2, {0.0, 1.0, 0.0}, 0.1, {0.0, 0.0, 1.0});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{-1.0, 2.0, 3.0};
    auto result = getPlaneIntersectionProportionPoints(common, 0.1, {10.0, 10.0, 0.0}, 0.2, {0.0, 10.0, 10.0}, 0.3, {10.0, 0.0, 10.0});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{-1.0, 2.0, 3.0};
    auto result = getPlaneIntersectionProportionPoints(common, 0.1, {-10.0, 10.0, 0.0}, 0.2, {0.0, -10.0, 10.0}, 0.3, {10.0, 0.0, 10.0});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{-1.0, 2.0, 3.0};
    auto result = getPlaneIntersectionProportionPoints(common, 0.1, {10.0, 10.0, 0.0}, 0.2, {0.0, -10.0, 10.0}, 0.3, {-10.0, 0.0, 10.0});
    EXPECT_TRUE(eq(common, result));
  }
}

Vertex getPlaneIntersectionVertices(Vertex aCommon,
                                         Vertex aOther1,
                                         Vertex aOther2,
                                         Vertex aOther3) {
  Plane plane1 = Plane::createFrom3points(aOther1, aOther2, aCommon);
  Plane plane2 = Plane::createFrom3points(aOther2, aOther3, aCommon);
  Plane plane3 = Plane::createFrom3points(aOther1, aOther3, aCommon);
  return Plane::intersect(plane1, plane2, plane3);
}

TEST(planeIntersection, Vertices) {
  {
    Vertex common{1.0, 2.0, 3.0};
    auto result = getPlaneIntersectionVertices(common, {10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{1.0, 2.0, 3.0};
    auto result = getPlaneIntersectionVertices(common, {-10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{1.0, 2.0, 3.0};
    auto result = getPlaneIntersectionVertices(common, {-10.0, 0.0, 0.0}, {0.0, -10.0, 0.0}, {0.0, 0.0, 10.0});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{1.0, 2.0, 3.0};
    auto result = getPlaneIntersectionVertices(common, {-10.0, 0.0, 0.0}, {0.0, -10.0, 0.0}, {0.0, 0.0, -10.0});
    EXPECT_TRUE(eq(common, result));
  }
}

Vertex getPlaneIntersection1vector2points(Vertex aCommon,
                                                Vertex aOne1, Vector aDir1,
                                                Vertex aOne2, Vector aDir2,
                                                Vertex aOne3, Vector aDir3) {
  Plane plane1 = Plane::createFrom1vector2points(aDir1, aOne1, aCommon);
  Plane plane2 = Plane::createFrom1vector2points(aDir2, aOne2, aCommon);
  Plane plane3 = Plane::createFrom1vector2points(aDir3, aOne3, aCommon);
  return Plane::intersect(plane1, plane2, plane3);
}

TEST(planeIntersection, VectorPoints) {
  {
    Vertex common{1.0, 2.0, -3.0};
    auto result = getPlaneIntersection1vector2points(common, {10.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                                                             {0.0, 10.0, 0.0}, {0.0, 0.0, 1.0},
                                                             {0.0, 0.0, 10.0}, {1.0, 0.0, 0.0});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{1.0, 2.0, -3.0};
    auto result = getPlaneIntersection1vector2points(common, {10.0, 0.0, 0.0}, {0.0, 1.0, 1.0},
                                                             {0.0, 10.0, 0.0}, {1.0, 0.0, -1.0},
                                                             {0.0, 0.0, 10.0}, {1.0, 1.0, 0.0});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{1.0, 2.0, -3.0};
    auto result = getPlaneIntersection1vector2points(common, {10.0, 0.0, 0.0}, {-4.0, 1.0, 1.0},
                                                             {0.0, 10.0, 0.0}, {1.0, -4.0, -1.0},
                                                             {0.0, 0.0, 10.0}, {1.0, 1.0, -4.0});
    EXPECT_TRUE(eq(common, result));
  }
}

Vertex getPlaneIntersection2vectors1point(Vertex aCommon,
                                                Vector aOne1, Vector aOther1,
                                                Vector aOne2, Vector aOther2,
                                                Vector aOne3, Vector aOther3) {
  Plane plane1 = Plane::createFrom2vectors1point(aOne1, aOther1, aCommon);
  Plane plane2 = Plane::createFrom2vectors1point(aOne2, aOther2, aCommon);
  Plane plane3 = Plane::createFrom2vectors1point(aOne3, aOther3, aCommon);
  return Plane::intersect(plane1, plane2, plane3);
}

TEST(planeIntersection, VectorsPoint) {
  {
    Vertex common{1.0, 2.0, -3.0};
    auto result = getPlaneIntersection1vector2points(common, {10.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                                                             {0.0, 10.0, 0.0}, {0.0, 0.0, 1.0},
                                                             {0.0, 0.0, 10.0}, {1.0, 0.0, 0.0});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{1.0, 2.0, -3.0};
    auto result = getPlaneIntersection1vector2points(common, {10.0, 1.0, 0.0}, {1.0, 10.0, 0.0},
                                                             {0.0, 10.0, 1.0}, {0.0, 1.0, 10.0},
                                                             {1.0, 0.0, 10.0}, {10.0, 0.0, 1.0});
    EXPECT_TRUE(eq(common, result));
  }
}

TEST(planeIntersection, Ray) {
  {
    Ray ray(Vertex{1.0, 2.0, -3.0}, Vector{1.0, 1.0, 1.0});
    Plane plane = Plane::createFrom3points(Vertex{10.0, 1.0, 2.0}, Vertex{11.0, 11.1, 2.0}, Vertex{12.0, 1.1, 4.4});
    auto result = plane.intersect(ray);
    EXPECT_TRUE(result.mValid);
  }
  {
    Ray ray(Vertex{1.0, 2.0, -3.0}, Vector{-1.0, 2.0, 3.0});
    Plane plane = Plane::createFrom3points(Vertex{10.0, 1.0, 2.0}, Vertex{11.0, 11.1, 2.0}, Vertex{12.0, 1.1, 4.4});
    auto result = plane.intersect(ray);
    EXPECT_TRUE(result.mValid);
    EXPECT_LT(result.mDistance, 0.0);
  }
  {
    Ray ray(Vertex{1.0, 2.0, -3.0}, Vector{0.0, 2.0, 0.0});
    Plane plane = Plane::createFrom3points(Vertex{10.0, 1.0, 2.0}, Vertex{10.0, 11.1, 2.0}, Vertex{10.0, 1.1, 4.4});
    auto result = plane.intersect(ray);
    EXPECT_FALSE(result.mValid);
  }
  {
    Ray ray(Vertex{1.0, 2.0, -3.0}, Vector{0.0, 2.0, 0.0});
    Plane plane = Plane::createFrom3points(Vertex{10.0, 10.0, 2.0}, Vertex{0.0, 10.0, 2.0}, Vertex{10.0, 10.0, 10.4});
    auto result = plane.intersect(ray);
    EXPECT_TRUE(result.mValid);
    EXPECT_LT((result.mPoint - Vertex{1.0, 10.0, -3.0}).norm(), 0.00001);
    EXPECT_GT(::abs(result.mCosIncidence), 0.9999);
  }
}

TEST(planeProjection, Point) {
  {
    Vertex point{0.0, 0.0, 0.0};
    Plane plane = Plane::createFrom3points({2.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 2.0});
    auto projected = plane.project(point);
    EXPECT_TRUE(eq(projected, Vertex(0.666666, 0.666666, 0.666666)));
  }
  {
    Vertex point{0.0, 0.0, 0.0};
    Plane plane = Plane::createFrom3points({2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}, {2.0, 0.0, 1.0});
    auto projected = plane.project(point);
    EXPECT_TRUE(eq(projected, Vertex(2.0, 0.0, 0.0)));
  }
  {
    Vertex point{1.0, 2.0, 3.0};
    Plane plane = Plane::createFrom3points({3.0, 2.0, 3.0}, {1.0, 4.0, 3.0}, {1.0, 2.0, 5.0});
    auto projected = plane.project(point);
    EXPECT_TRUE(eq(projected, Vertex(1.666666, 2.666666, 3.666666)));
  }
  {
    Vertex point{-1.0, -2.0, 3.0};
    Plane plane = Plane::createFrom3points({1.0, -2.0, 3.0}, {1.0, -3.0, 3.0}, {1.0, -2.0, 4.0});
    auto projected = plane.project(point);
    EXPECT_TRUE(eq(projected, Vertex(1.0, -2.0, 3.0)));
  }
  {
    Vertex point{1.666666, 2.666666, 3.666666};
    Plane plane = Plane::createFrom3points({3.0, 2.0, 3.0}, {1.0, 4.0, 3.0}, {1.0, 2.0, 5.0});
    auto projected = plane.project(point);
    EXPECT_TRUE(eq(projected, point));
  }
}

TEST(planeDistance, Point) {
  {
    Vertex point{0.0, 0.0, 0.0};
    Plane plane = Plane::createFrom3points({2.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 2.0});
    auto dist = ::abs(plane.distance(point));
    EXPECT_TRUE(eq(dist, 1.15468));
  }
  {
    Vertex point{0.0, 0.0, 0.0};
    Plane plane = Plane::createFrom3points({2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}, {2.0, 0.0, 1.0});
    auto dist = ::abs(plane.distance(point));
    EXPECT_TRUE(eq(dist, 2.0));
  }
  {
    Vertex point{1.0, 2.0, 3.0};
    Plane plane = Plane::createFrom3points({3.0, 2.0, 3.0}, {1.0, 4.0, 3.0}, {1.0, 2.0, 5.0});
    auto dist = ::abs(plane.distance(point));
    EXPECT_TRUE(eq(dist, 1.15468));
  }
  {
    Vertex point{-1.0, -2.0, 3.0};
    Plane plane = Plane::createFrom3points({1.0, -2.0, 3.0}, {1.0, -3.0, 3.0}, {1.0, -2.0, 4.0});
    auto dist = ::abs(plane.distance(point));
    EXPECT_TRUE(eq(dist, 2.0));
  }
  {
    Vertex point{1.666666, 2.666666, 3.666666};
    Plane plane = Plane::createFrom3points({3.0, 2.0, 3.0}, {1.0, 4.0, 3.0}, {1.0, 2.0, 5.0});
    auto dist = ::abs(plane.distance(point));
    EXPECT_TRUE(eq(dist, 0.0));
  }
}

TEST(toWhichSide, Points) {
  Vertex triangle0{3.0, 2.0, 5.0};
  Vertex triangle1{1.0, 4.0, 5.0};
  Vertex triangle2{6.0, 5.0, 5.0};
  Vertex start = (triangle0 + triangle1 + triangle2) / 3.0;
  {
    Vertex end = start + Vector{1.0, 0.0, 0.0};
    auto converter = util::getBarycentricInverse(triangle0, triangle1, triangle2);
    EXPECT_EQ(util::toWhichSide(converter * start, converter * end), 2u);
  }
  {
    Vertex end = start + Vector{0.0, 1.0, 0.0};
    auto converter = util::getBarycentricInverse(triangle0, triangle1, triangle2);
    EXPECT_EQ(util::toWhichSide(converter * start, converter * end), 1u);
  }
  {
    Vertex end = start + Vector{-1.0, -1.0, 0.0};
    auto converter = util::getBarycentricInverse(triangle0, triangle1, triangle2);
    EXPECT_EQ(util::toWhichSide(converter * start, converter * end), 0u);
  }
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
