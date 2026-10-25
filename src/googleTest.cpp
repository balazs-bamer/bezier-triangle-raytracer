#include "3dGeomUtil.h"
#include "bezierMesh.h"

#include<set>
#include<deque>
#include<random>
#include<string>
#include<vector>
#include<iostream>
#include<functional>

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

std::string getException(std::function<void()> const &aCode) {
  try {
    aCode();
  }
  catch(char const * const aMessage) {
    return aMessage;
  }
  return "";
}

// Open grid of aCount x aCount squares over [-1, 1] x [-1, 1] on the surface z = aCurvature * (x^2 + y^2)
Mesh makeGrid(int32_t const aCount, double const aCurvature) {
  auto vertex = [aCount, aCurvature](int32_t const aI, int32_t const aJ) {
    double const x = 2.0 * aI / aCount - 1.0;
    double const y = 2.0 * aJ / aCount - 1.0;
    return Vertex{x, y, aCurvature * (x * x + y * y)};
  };
  Mesh result;
  for(int32_t i = 0; i < aCount; ++i) {
    for(int32_t j = 0; j < aCount; ++j) {
      result.push_back(Triangle{vertex(i, j), vertex(i + 1, j), vertex(i + 1, j + 1)});
      result.push_back(Triangle{vertex(i, j), vertex(i + 1, j + 1), vertex(i, j + 1)});
    }
  }
  return result;
}

// Keeps the triangles of a closed solid of revolution above z = 0, and gives their original indices too.
std::pair<Mesh, std::vector<uint32_t>> makeUpperHalf(Mesh const &aClosed) {
  Mesh result;
  std::vector<uint32_t> indices;
  for(uint32_t i = 0u; i < aClosed.size(); ++i) {
    auto const &triangle = aClosed[i];
    if(triangle[0u](2) + triangle[1u](2) + triangle[2u](2) > 0.0) {
      result.push_back(triangle);
      indices.push_back(i);
    }
    else { // Nothing to do
    }
  }
  result.standardizeNormals(Vector{0.0, 0.0, 1.0});
  return std::make_pair(result, indices);
}

uint32_t countBoundaryEdges(Mesh const &aMesh) {
  uint32_t result = 0u;
  for(auto const &neighbours : aMesh.getFace2neighbours()) {
    result += std::count(neighbours.mFellowTriangles.cbegin(), neighbours.mFellowTriangles.cend(), Mesh::csNoNeighbour);
  }
  return result;
}

uint32_t countCoreFaces(Mesh const &aMesh) {
  uint32_t result = 0u;
  for(uint32_t i = 0u; i < aMesh.size(); ++i) {
    result += (aMesh.isCoreFace(i) ? 1u : 0u);
  }
  return result;
}

// Consistent means each edge is traversed in opposite directions by its two triangles.
bool isConsistentlyOriented(Mesh const &aMesh) {
  std::set<std::array<double, 6u>> directedEdges;
  bool result = true;
  for(auto const &triangle : aMesh) {
    for(uint32_t i = 0u; i < 3u; ++i) {
      auto const &start = triangle[i];
      auto const &end = triangle[(i + 1u) % 3u];
      result = result && directedEdges.insert({start(0), start(1), start(2), end(0), end(1), end(2)}).second;
    }
  }
  return result;
}

TEST(openMesh, Grid) {
  int32_t const count = 6;
  auto mesh = makeGrid(count, 0.2);
  mesh.standardizeVertices();
  mesh.standardizeNormals(Vector{0.0, 0.0, -1.0});
  EXPECT_FALSE(mesh.isClosed());
  EXPECT_EQ(countBoundaryEdges(mesh), 4u * count);
  EXPECT_EQ(mesh.getBoundaryVertices().size(), 4u * count);
  EXPECT_EQ(countCoreFaces(mesh), 2u * (count - 2) * (count - 2));
  EXPECT_TRUE(isConsistentlyOriented(mesh));
  for(auto const &triangle : mesh) {
    EXPECT_LT(util::getNormal(triangle)(2), 0.0);
  }
  BezierMesh bezier(mesh);
  EXPECT_EQ(bezier.size(), 3u * countCoreFaces(mesh));
}

TEST(openMesh, SingleTriangle) {
  Mesh mesh;
  mesh.push_back(Triangle{Vertex{0.0, 0.0, 0.0}, Vertex{1.0, 0.0, 0.0}, Vertex{0.0, 1.0, 0.0}});
  mesh.standardizeNormals(Vector{0.0, 0.0, 1.0});
  EXPECT_EQ(countBoundaryEdges(mesh), 3u);
  EXPECT_EQ(mesh.getBoundaryVertices().size(), 3u);
  EXPECT_EQ(countCoreFaces(mesh), 0u);
  EXPECT_GT(util::getNormal(mesh[0u])(2), 0.0);
  BezierMesh bezier(mesh);
  EXPECT_EQ(bezier.size(), 0u);
}

TEST(openMesh, ClosedSphere) {
  Mesh mesh;
  mesh.makeUnitSphere(7, 3);
  mesh.standardizeVertices();
  mesh.standardizeNormals(Vector{-1.0, 0.0, 0.0});
  EXPECT_TRUE(mesh.isClosed());
  EXPECT_EQ(countBoundaryEdges(mesh), 0u);
  EXPECT_EQ(countCoreFaces(mesh), mesh.size());
  EXPECT_TRUE(isConsistentlyOriented(mesh));
}

TEST(openMesh, Errors) {
  auto grid = makeGrid(3, 0.0);
  EXPECT_EQ(getException([&grid](){ grid.standardizeNormals(); }), "Vertex on edge detected.");

  Vertex const a{0.0, 0.0, 0.0};
  Vertex const b{2.0, 0.0, 0.0};
  Mesh nonManifold;
  nonManifold.push_back(Triangle{a, b, Vertex{1.0, 1.0, 0.0}});
  nonManifold.push_back(Triangle{b, a, Vertex{1.0, -1.0, 0.0}});
  nonManifold.push_back(Triangle{a, b, Vertex{1.0, 0.0, 1.0}});
  EXPECT_EQ(getException([&nonManifold](){ nonManifold.standardizeNormals(Vector{0.0, 0.0, 1.0}); }), "Edge belonging to more than two triangles detected.");

  Vertex const middle{1.0, 0.0, 0.0};
  Vertex const below{1.0, -1.0, 0.0};
  Mesh tJunction;
  tJunction.push_back(Triangle{a, b, Vertex{1.0, 1.0, 0.0}});
  tJunction.push_back(Triangle{a, below, middle});
  tJunction.push_back(Triangle{middle, below, b});
  EXPECT_EQ(getException([&tJunction](){ tJunction.standardizeNormals(Vector{0.0, 0.0, 1.0}); }), "Vertex on edge detected.");

  Mesh twoPieces;
  twoPieces.push_back(Triangle{a, b, Vertex{1.0, 1.0, 0.0}});
  twoPieces.push_back(Triangle{Vertex{10.0, 0.0, 0.0}, Vertex{12.0, 0.0, 0.0}, Vertex{11.0, 1.0, 0.0}});
  EXPECT_EQ(getException([&twoPieces](){ twoPieces.standardizeNormals(Vector{0.0, 0.0, 1.0}); }), "Mesh consisting of more than one piece detected.");
}

// Surface normal of a Bezier triangle by central differences, also outside the triangle.
Vector getBezierNormal(BezierTriangle const &aBezier, Vertex const &aBary) {
  double const h = 1e-6;
  Vector const du = (aBezier.interpolate(aBary(0) + h, aBary(1), aBary(2) - h) - aBezier.interpolate(aBary(0) - h, aBary(1), aBary(2) + h)) / (2.0 * h);
  Vector const dv = (aBezier.interpolate(aBary(0), aBary(1) + h, aBary(2) - h) - aBezier.interpolate(aBary(0), aBary(1) - h, aBary(2) + h)) / (2.0 * h);
  return du.cross(dv).normalized();
}

// Largest angle in degrees between the normals of neighbouring Bezier triangles along their common sides,
// first inside the original triangles, second across the original sides.
std::pair<double, double> getMaxNormalJumps(BezierMesh const &aBezier) {
  double inside = 0.0;
  double across = 0.0;
  for(uint32_t i = 0u; i < aBezier.size(); ++i) {
    auto const &bezier = aBezier[i];
    auto const neighbours = bezier.getNeighbours();
    for(double s = 0.05; s < 1.0; s += 0.1) {
      // The same point in own and neighbour barycentric coordinates: the original side is reversed in the neighbour,
      // (V1, centre) is (V0, centre) in the next triangle, and (V0, centre) is (V1, centre) in the previous one.
      std::array<std::pair<Vertex, Vertex>, 3u> const sides = {{ {{1.0 - s, s, 0.0}, {s, 1.0 - s, 0.0}},
                                                                 {{0.0, 1.0 - s, s}, {1.0 - s, 0.0, s}},
                                                                 {{1.0 - s, 0.0, s}, {0.0, 1.0 - s, s}} }};
      for(uint32_t side = 0u; side < 3u; ++side) {
        if(neighbours[side] != Mesh::csNoNeighbour) {
          double const cosine = getBezierNormal(bezier, sides[side].first).dot(getBezierNormal(aBezier[neighbours[side]], sides[side].second));
          double &worst = (side == 0u ? across : inside);
          worst = std::max(worst, ::acos(std::min(1.0, cosine)) * 180.0 / cgPi);
        }
        else { // Nothing to do
        }
      }
    }
  }
  return std::make_pair(inside, across);
}

TEST(openBezierMesh, C1) {
  auto grid = makeGrid(8, 0.3);                                         // Open meshes are height fields along their average normal, exactly C1.
  grid.standardizeVertices();
  grid.standardizeNormals(Vector{0.0, 0.0, 1.0});
  auto const [insideGrid, acrossGrid] = getMaxNormalJumps(BezierMesh(grid));
  EXPECT_LT(insideGrid, 1e-4);
  EXPECT_LT(acrossGrid, 1e-4);

  Mesh sphere;
  sphere.makeUnitSphere(21, 15);
  sphere.standardizeVertices();
  sphere.standardizeNormals();
  auto const [half, half2sphere] = makeUpperHalf(sphere);
  auto const [insideHalf, acrossHalf] = getMaxNormalJumps(BezierMesh(half));
  EXPECT_LT(insideHalf, 1e-4);
  EXPECT_LT(acrossHalf, 1e-4);

  auto const [insideClosed, acrossClosed] = getMaxNormalJumps(BezierMesh(sphere)); // Closed meshes are C1 inside the original triangles,
  EXPECT_LT(insideClosed, 1e-4);                                                    // and only approximately across their sides.
  EXPECT_LT(acrossClosed, 0.2);
}

TEST(openBezierMesh, Intersect) {
  Mesh sphere;
  sphere.makeUnitSphere(21, 15);
  sphere.standardizeVertices();
  sphere.standardizeNormals();
  auto [half, half2sphere] = makeUpperHalf(sphere);
  BezierMesh open(half);
  double zMinCore = std::numeric_limits<double>::max();
  for(uint32_t i = 0u; i < half.size(); ++i) {
    if(open.getBezierIndex(i) != Mesh::csNoNeighbour) {
      for(auto const &vertex : half[i]) {
        zMinCore = std::min(zMinCore, vertex(2));
      }
    }
    else { // Nothing to do
    }
  }

  std::mt19937 generator(5u);
  std::uniform_real_distribution<double> distribution(-1.0, 1.0);
  uint32_t checkedInside = 0u;
  for(uint32_t i = 0u; i < 2000u; ++i) {
    Vector direction{distribution(generator), distribution(generator), distribution(generator)};
    direction.normalize();
    Ray const ray(direction * 3.0, -direction);
    auto const hitOpen = open.intersect(ray);
    if(hitOpen.mWhat == BezierIntersection::What::cIntersect) {        // Never beyond the rim.
      EXPECT_GT(hitOpen.mIntersection.mPoint(2), zMinCore - 0.02);
    }
    else { // Nothing to do
    }
    if(direction(2) > zMinCore + 0.2) {                                // Well inside, must be on the unit sphere.
      ASSERT_EQ(hitOpen.mWhat, BezierIntersection::What::cIntersect);
      EXPECT_LT(::abs(hitOpen.mIntersection.mPoint.norm() - 1.0), 0.005);
      ++checkedInside;
    }
    else { // Nothing to do
    }
  }
  EXPECT_GT(checkedInside, 100u);
}

TEST(openBezierMesh, SplitKeepsRim) {
  Mesh sphere;
  sphere.makeUnitSphere(7, 3);
  sphere.standardizeVertices();
  sphere.standardizeNormals();
  auto [half, half2sphere] = makeUpperHalf(sphere);
  BezierMesh open(half);
  auto refined = open.splitThickBezierTriangles();
  EXPECT_GT(refined.size(), half.size());
  refined.standardizeVertices();
  EXPECT_EQ(getException([&refined](){ refined.standardizeNormals(Vector{0.0, 0.0, 1.0}); }), "");
  EXPECT_EQ(refined.getBoundaryVertices(), half.getBoundaryVertices());
  EXPECT_TRUE(isConsistentlyOriented(refined));
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
