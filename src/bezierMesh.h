#ifndef BEZIER_TRIANGLE_RAYTRACER_BEZIERMESH
#define BEZIER_TRIANGLE_RAYTRACER_BEZIERMESH

#include "mesh.h"
#include "bezierTriangle.h"
#include <functional>
#include <numeric>


class BezierMesh final {
private:
  static constexpr std::array<double, 3u> csSampleRatiosOriginalSide = { 0.25, 0.5, 0.75 };
  static constexpr double                 csBezierHeightPerPerimeterLimit = 0.03;
  static constexpr double                 csSplitBezierInterpolateFactor  = 0.7;   // For triangle splitting, new vertex is computed as
                                                                                   // barycentricSplit * (1.0 - csSBIF) + interpolate(barycentricSplit) * csSBIF
  static constexpr double                 csRimBarycentricBand            = 1.0;   // See intersectBeyondRim
  static constexpr double                 csMinBarycentricOfHit           = -0.2;  // Newton may converge to the polynomial continuation of the patch, which
                                                                                   // BezierTriangle::classifyIntersection lets through when it is outside
                                                                                   // two divider planes. Real hits were not seen below -0.1.

  struct Rim final {           // For Bezier triangles whose neighbour across the original edge is not built.
    Plane  mUnderlyingPlane;   // Same as in the BezierTriangle.
    Matrix mBarycentricInverse;// Same as in the BezierTriangle.
    Plane  mRimPlane;          // The plane between the original neighbours, positive inside.
  };

  std::vector<BezierTriangle>    mMesh;
  typename Mesh::TheMesh         mOriginalTriangles;
  typename Mesh::Face2neighbours mOriginalNeighbours;
  std::vector<uint32_t>          mFace2bezierBase;  // Index of the first of the 3 Bezier triangles of each original triangle, Mesh::csNoNeighbour if not built.
  std::vector<uint32_t>          mBuiltFaces;       // Original triangle index for each 3 consecutive Bezier triangles.
  std::vector<uint32_t>          mBezier2rim;       // Index in mRims for each Bezier triangle, Mesh::csNoNeighbour if not on the rim.
  std::vector<Rim>               mRims;

  using Neighbours          = typename Mesh::Neighbours;
  using MissingFieldsMethod = std::function<void(BezierTriangle &, Vertex const &aOriginalCentroid, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious)>;

public:
  // Builds Bezier triangles for the core triangles of aMesh, see Mesh::isCoreFace. For a closed mesh these are all.
  // The other triangles of an open mesh only serve as neighbours, so the result is exactly what the same triangles would give in a bigger mesh.
  BezierMesh(Mesh const &aMesh);
  // The same, but only for the core triangles where aMask is true.
  BezierMesh(Mesh const &aMesh, std::vector<bool> const &aMask);

  auto size() const { return mMesh.size(); }
  auto cbegin() const { return mMesh.cbegin(); }
  auto cend() const { return mMesh.cend(); }
  auto const &operator[](uint32_t const aI) const { return mMesh[aI]; }
  // Index of the first of the 3 Bezier triangles made of the original triangle, or Mesh::csNoNeighbour if not built.
  uint32_t getBezierIndex(uint32_t const aIndexOriginal) const { return mFace2bezierBase[aIndexOriginal]; }

  Mesh interpolate(int32_t const aDivisor) const;
  std::vector<Vertex> dumpControlPoints() const;
  Mesh splitThickBezierTriangles() const;        // TODO consider automatic mechanism to iterate until either
                                                        // - no more split or under a percentage
                                                        // - absolute height left is below a value
                                                        // - minimal triangle perimeter is below a value or percentage of the original
                                                        // Triangles not built are copied, and split if needed to match their built neighbours.
  BezierIntersection intersect(Ray const &aRay) const;

private:
  void setMissingFields(Mesh const &aMesh, MissingFieldsMethod aMissingFieldsMethod);
  Vertex getSplitVertex(uint32_t const aIndexOriginal, uint32_t const aSide) const;
  static void append2split(Mesh &aResult, Triangle const &aOriginalTriangle, uint8_t const aSplit, std::array<Vertex, 3u> const &aSplitVertices);
  static void append3split(Mesh &aResult, Triangle const &aOriginalTriangle, uint8_t const aSplit, std::array<Vertex, 3u> const &aSplitVertices);
  static void append4split(Mesh &aResult, Triangle const &aOriginalTriangle, std::array<Vertex, 3u> const &aSplitVertices);
  Vertex interpolate(uint32_t const aIndex, double const aBary0, double const aBary1, double const aBary2) const;
  BezierIntersection intersectBeyondRim(Rim const &aRim, BezierTriangle const &aBezier, Ray const &aRay) const;
};

#endif
