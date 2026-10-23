#ifndef BEZIER_TRIANGLE_RAYTRACER_BEZIER
#define BEZIER_TRIANGLE_RAYTRACER_BEZIER

#include "3dGeomUtil.h"


struct BezierIntersection final {
  enum class What : uint32_t {
    cFollowSide0 = 0u,
    cFollowSide1 = 1u,
    cFollowSide2 = 2u,
    cNone        = 3u,
    cIntersect   = 4u
  };

  Intersection mIntersection;
  Vertex       mBarycentric;
  Vector       mNormal;
  What         mWhat;
};

class BezierTriangle final {    // Cubic Bezier triangle
public:
  enum class LimitPlaneIntersection : uint32_t {
    cThis = 0u,
    cNone = 1u
  };

  static constexpr uint32_t csControlPointsSize                       = 10u;
private:
  static constexpr uint32_t csControlIndexOriginalVertex0             = 0u;
  static constexpr uint32_t csControlIndexOriginalVertex1             = 1u;
  static constexpr uint32_t csControlIndexAboveOriginalCentroid       = 2u;
  static constexpr uint32_t csControlIndexOnOriginalSide0             = 3u;
  static constexpr uint32_t csControlIndexOnOriginalSide1             = 4u;
  static constexpr uint32_t csControlIndexOnSideToOriginalCentroid0   = 5u;
  static constexpr uint32_t csControlIndexOnSideToOriginalCentroid1   = 6u;
  static constexpr uint32_t csControlIndexOnSideFromOriginalCentroid0 = 7u;
  static constexpr uint32_t csControlIndexOnSideFromOriginalCentroid1 = 8u;
  static constexpr uint32_t csControlIndexMiddle                      = 9u;

  static constexpr uint32_t csControlIndex300 = csControlIndexOriginalVertex0;
  static constexpr uint32_t csControlIndex030 = csControlIndexOriginalVertex1;
  static constexpr uint32_t csControlIndex003 = csControlIndexAboveOriginalCentroid;
  static constexpr uint32_t csControlIndex210 = csControlIndexOnOriginalSide0;
  static constexpr uint32_t csControlIndex120 = csControlIndexOnOriginalSide1;
  static constexpr uint32_t csControlIndex021 = csControlIndexOnSideToOriginalCentroid0;
  static constexpr uint32_t csControlIndex012 = csControlIndexOnSideToOriginalCentroid1;
  static constexpr uint32_t csControlIndex102 = csControlIndexOnSideFromOriginalCentroid0;
  static constexpr uint32_t csControlIndex201 = csControlIndexOnSideFromOriginalCentroid1;
  static constexpr uint32_t csControlIndex111 = csControlIndexMiddle;

  // Power basis of the surface in u = bary0, v = bary1: S(u, v) = sum A_ij * u^i * v^j for i + j <= 3
  static constexpr uint32_t csPowerCoefficientsSize = 10u;
  static constexpr uint32_t csPowerIndex00 = 0u;
  static constexpr uint32_t csPowerIndex10 = 1u;
  static constexpr uint32_t csPowerIndex20 = 2u;
  static constexpr uint32_t csPowerIndex30 = 3u;
  static constexpr uint32_t csPowerIndex01 = 4u;
  static constexpr uint32_t csPowerIndex11 = 5u;
  static constexpr uint32_t csPowerIndex21 = 6u;
  static constexpr uint32_t csPowerIndex02 = 7u;
  static constexpr uint32_t csPowerIndex12 = 8u;
  static constexpr uint32_t csPowerIndex03 = 9u;

  static constexpr double    csProportionControlOnOriginalSide           = 0.291;
  static constexpr double    csProportionControlOnOriginalVertexCentroid = 0.304;
  static constexpr double    csProportionControlOnOriginalMedian         = 0.2;
  static constexpr double    csHeightSafetyFactor                        = 1.33333333;
  static constexpr double    csOneThird                                  = 1.0 / 3.0;
  static constexpr uint32_t  csRootSearchIterations4bezierEvaluations    = 4u;
  static constexpr int32_t   csHeightSampleDivisor                       = 5;
  static constexpr double    csMaxIntersectionDistanceFromRay            = 0.01;
  static constexpr double    csMinimalRayDistance                        = 1.0;
  static constexpr double    csIntersectionEstimationEpsilon             = 0.000001;
  static constexpr uint32_t  csRootSearchIterations4kajiya               = 12u;
  static constexpr double    csKajiyaBarycentricEpsilon                  = 1e-12;  // Newton stops when the barycentric step is smaller.
  static constexpr double    csKajiyaDeterminantEpsilon                  = 1e-20;  // Below this the ray is regarded parallel to the surface.

  Plane                    mUnderlyingPlane;         // Plane span by mControlPoints[0-2], normal is normalized and points outwards.
  std::array<Plane, 3u>    mNeighbourDividerPlanes;  // For indices 1 and 2, these are about the plane going through each edge in direction the average of the adjacent side normals.
                                                     // Index 0 will be aPlaneBetweenOriginalNeighbours
                                                     // Their normal is such that for all planes any point P belonging to this Bezier triangle has plane.distance(P) >= 0
  std::array<uint32_t, 3u> mNeighbours;              // With new indices after Clough-Tocher split.
  // 0-1 identical to original triangle vertices, named aOriginalCommonVertex0, aOriginalCommonVertex1 in constructor
  // 2 somewhere above the middle of the original triangle center, to be calculated
  // 3-8 somewhere above sides of triangle (0-2): 3i+3, 3i+4 on the same side as (i, (i+1)%3), to be calculated
  // 9 in the middle, to be calculated
  std::array<Vertex, csControlPointsSize>  mControlPoints;
  std::array<Vector, csPowerCoefficientsSize> mPowerCoefficients; // Computed from mControlPoints, see csPowerIndex*. Derivatives are taken from these analytically.
  Matrix                   mBarycentricInverse;      // T = (v1 v2 v3), b = barycentric coefficient column, v = point on plane, v=Tb, b = mBI * v
                                                     // multiplication by the inverse proven to be much faster than solving the linear equation under Eigen.
  double                   mHeightInside;            // Sampled biggest distance of the Bezier triangle measured from the underlying triangle, < 0
  double                   mHeightOutside;           // this is > 0
  Vector                   mBezierDerivativeDirectionVectorA; // We calculate directional derivatives in two perpendicular directions
  Vector                   mBezierDerivativeDirectionVectorB; // and their cross product will yield the Bezier surface normal

public:
  // Vertices in argument are in order such that the normal points to the desired direction.
  // Neighbour i is neighbour of edge (i, i + 1)
  BezierTriangle(Vertex const &aOriginalCommonVertex0, Vertex const &aOriginalCommonVertex1, Vertex const &aOriginalCentroid,
                 Vector const &aAverageNormal0, Vector const &aAverageNormal1, Plane const &aPlaneBetweenOriginalNeighbours,
                 std::array<uint32_t, 3u> const &aNeighbourIndices);

  void setMissingFields1(Vertex const &aOriginalCentroid, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious);
  void setMissingFields2(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &);
  void setMissingFields3(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious);

  Vertex getControlPoint(uint32_t const aI) const { return mControlPoints[aI]; }
  std::array<uint32_t, 3u> getNeighbours() const { return mNeighbours; }           // TODO remove

  Vertex interpolateLinear(double const aBary0, double const aBary1, double const aBary2) const;
  Vertex interpolateLinear(double const aBary0, double const aBary1) const { return interpolateLinear(aBary0, aBary1, 1.0 - aBary0 - aBary1); }
  Vertex interpolateLinear(Vertex const &aBary) const { return interpolateLinear(aBary(0u), aBary(1u), aBary(2u)); }

  Vertex interpolate(double const aBary0, double const aBary1, double const aBary2) const;
  Vertex interpolate(double const aBary0, double const aBary1) const { return interpolate(aBary0, aBary1, 1.0 - aBary0 - aBary1); }
  Vertex interpolate(Vertex const &aBary) const { return interpolate(aBary(0u), aBary(1u), aBary(2u)); }
  Vertex interpolateAboveOriginalCentroid() const { return mControlPoints[csControlIndexAboveOriginalCentroid]; }

  BezierIntersection intersectByBezierEvaluations(Ray const &aRay, LimitPlaneIntersection const aShouldLimitPlaneIntersection) const;
  // Ray is regarded as the intersection of two planes containing it, and Newton's method solves for the barycentric coordinates
  // where the Bezier surface touches both planes. Uses the power basis and its analytical derivatives.
  BezierIntersection intersectKajiya(Ray const &aRay, LimitPlaneIntersection const aShouldLimitPlaneIntersection) const;

private:
  Vector getNormal(Vector const &aBarycentric) const;
  void computePowerCoefficients();
  // Sets aResult.mWhat using mNeighbourDividerPlanes, and mCosIncidence for an own intersection. Needs mPoint and mNormal already set.
  void classifyIntersection(BezierIntersection &aResult, Ray const &aRay) const;

  // Works both for the vector coefficients and for their projections to a plane normal.
  template<typename tCoefficient>
  static void evaluatePower(std::array<tCoefficient, csPowerCoefficientsSize> const &aCoefficients, double const aU, double const aV,
                            tCoefficient &aValue, tCoefficient &aDerivativeU, tCoefficient &aDerivativeV);
};

#endif
