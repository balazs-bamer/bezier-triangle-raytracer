#include "bezierTriangle.h"


BezierTriangle::BezierTriangle(Vertex const &aOriginalCommonVertex0, Vertex const &aOriginalCommonVertex1, Vertex const &aOriginalCentroid,
                                      Vector const &aAverageNormal0, Vector const &aAverageNormal1, Plane const &aPlaneBetweenOriginalNeighbours,
                                      Vector const &aLiftDirection, std::array<uint32_t, 3u> const &aNeighbourIndices)
  : mNeighbours(aNeighbourIndices)
  , mLiftDirection(aLiftDirection) {
  mControlPoints[csControlIndexOriginalVertex0] = aOriginalCommonVertex0;
  mControlPoints[csControlIndexOriginalVertex1] = aOriginalCommonVertex1;

  Plane const commonPlaneVertex0{aAverageNormal0, aOriginalCommonVertex0.dot(aAverageNormal0) };
  Plane const commonPlaneVertex1{aAverageNormal1, aOriginalCommonVertex1.dot(aAverageNormal1) };

  // The thirds of the original side lifted onto the tangent planes of the vertices. The triangle on the other side computes exactly the same.
  mControlPoints[csControlIndexOnOriginalSide0] = liftOnto(commonPlaneVertex0, aOriginalCommonVertex0 + (aOriginalCommonVertex1 - aOriginalCommonVertex0) * csOneThird);
  mControlPoints[csControlIndexOnOriginalSide1] = liftOnto(commonPlaneVertex1, aOriginalCommonVertex1 + (aOriginalCommonVertex0 - aOriginalCommonVertex1) * csOneThird);

  mNeighbourDividerPlanes[0u] = aPlaneBetweenOriginalNeighbours;   // Oriented in setMissingFields1, when the middle control point is known.
}

// The rest of the control points follow the C1 conditions of the Clough-Tocher split. Inside the original triangle each
// control point next to a vertex shared by two split triangles is the average of its three neighbours there.
void BezierTriangle::setMissingFields1(Vertex const &aOriginalCentroid, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious) {
  auto const &vertex0 = mControlPoints[csControlIndex300];
  auto const &vertex1 = mControlPoints[csControlIndex030];
  auto const &side0   = mControlPoints[csControlIndex210];
  auto const &side1   = mControlPoints[csControlIndex120];
  // The summation order is the same as in the neighbour computing the same point, to get exactly the same result.
  mControlPoints[csControlIndex201] = (vertex0 + side0 + aTrianglePrevious.mControlPoints[csControlIndex120]) * csOneThird;
  mControlPoints[csControlIndex021] = (vertex1 + aTriangleNext.mControlPoints[csControlIndex210] + side1) * csOneThird;

  // The derivative across the original side at its middle must lie in the tangent plane there, which is common with the triangle
  // on the other side. Its normal is the one of the plane between the original neighbours, made perpendicular to the common
  // boundary curve. With the common tangent planes at the vertices, this makes the surface C1 across the original side.
  Vector sideNormal = (vertex1 - vertex0).cross(mNeighbourDividerPlanes[0u].mNormal);   // Lies in the plane between the original neighbours.
  if(sideNormal.dot(util::getNormal(vertex0, vertex1, aOriginalCentroid)) < 0.0) {
    sideNormal = -sideNormal;
  }
  else { // Nothing to do
  }
  Vector const tangent = (vertex1 + side1 - side0 - vertex0).normalized();               // Of the boundary curve in its middle.
  Vector const middleNormal = (sideNormal - tangent * sideNormal.dot(tangent)).normalized();
  // The derivative towards the centroid in the middle of the original side is 3/4 * (r201 + 2 * (111 - (210 + 120) / 2) + r021),
  // where r201 = 201 - (300 + 210) / 2 and r021 = 021 - (120 + 030) / 2. It must be perpendicular to middleNormal.
  Vector const r201 = mControlPoints[csControlIndex201] - (vertex0 + side0) / 2.0;
  Vector const r021 = mControlPoints[csControlIndex021] - (side1 + vertex1) / 2.0;
  double const wanted = middleNormal.dot((side0 + side1) / 2.0) - middleNormal.dot(r201 + r021) / 2.0;
  Vertex const flat = (vertex0 + vertex1 + aOriginalCentroid) * csOneThird;
  Vector const lift = (mLiftDirection.isZero() ? middleNormal : mLiftDirection);
  mControlPoints[csControlIndex111] = flat + lift * ((wanted - middleNormal.dot(flat)) / middleNormal.dot(lift));

  mNeighbourDividerPlanes[0u].makeDistancePositive(mControlPoints[csControlIndexMiddle]);
}

void BezierTriangle::setMissingFields2(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious) {
  mControlPoints[csControlIndex102] = (mControlPoints[csControlIndex201] + mControlPoints[csControlIndex111] + aTrianglePrevious.mControlPoints[csControlIndex111]) * csOneThird;
  mControlPoints[csControlIndex012] = (mControlPoints[csControlIndex021] + aTriangleNext.mControlPoints[csControlIndex111] + mControlPoints[csControlIndex111]) * csOneThird;
}

void BezierTriangle::setMissingFields3(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &) {
  mControlPoints[csControlIndexAboveOriginalCentroid] = (mControlPoints[csControlIndexOnSideFromOriginalCentroid0] +
                                                         mControlPoints[csControlIndexOnSideToOriginalCentroid1] +
                                                         aTriangleNext.mControlPoints[csControlIndexOnSideToOriginalCentroid1]) / 3.0;

  mUnderlyingPlane = Plane::createFrom3points(mControlPoints[csControlIndexOriginalVertex0], mControlPoints[csControlIndexOriginalVertex1], mControlPoints[csControlIndexAboveOriginalCentroid]);

  mBarycentricInverse = util::getBarycentricInverse(mControlPoints[csControlIndexOriginalVertex0], mControlPoints[csControlIndexOriginalVertex1], mControlPoints[csControlIndexAboveOriginalCentroid]);

  mHeightInside = 0.0;
  mHeightOutside = 0.0;
  Triangle barycentric{Vertex{{1.0, 0.0, 0.0}}, Vertex{{0.0, 1.0, 0.0}}, Vertex{{0.0, 0.0, 1.0}}};
  util::divide(barycentric, csHeightSampleDivisor, [this](Triangle && aNewBary) {
    for(uint32_t i = 0u; i < 3u; ++i) {
      auto distance = mUnderlyingPlane.distance(interpolate(aNewBary[i]));  // Will count most points 2 or 3 times, but don't care now.
      mHeightInside = std::min(mHeightInside, distance);
      mHeightOutside = std::max(mHeightOutside, distance);
    }
  } );
  mHeightInside *= csHeightSafetyFactor;
  mHeightOutside *= csHeightSafetyFactor;
  mBezierDerivativeDirectionVectorA = Vertex{1.0, 0.0, -1.0}; // Parallel to the side from 2 to 0, no matter how long.
  mBezierDerivativeDirectionVectorB = mBarycentricInverse *      // Perpendicular to the previous one.
      (mControlPoints[csControlIndexAboveOriginalCentroid] - mControlPoints[csControlIndexOriginalVertex0]).cross(mUnderlyingPlane.mNormal);
  computePowerCoefficients();
}

void BezierTriangle::computePowerCoefficients() {
  static constexpr std::array<double, 4u> csFactorials = { 1.0, 1.0, 2.0, 6.0 };
  static constexpr uint32_t cInvalid = csPowerCoefficientsSize;
  // Exponents of bary0, bary1, bary2 and the control point index belonging to them.
  static constexpr std::array<std::array<uint32_t, 4u>, csControlPointsSize> csControlExponents = {{
    { 3u, 0u, 0u, csControlIndex300 }, { 0u, 3u, 0u, csControlIndex030 }, { 0u, 0u, 3u, csControlIndex003 },
    { 2u, 1u, 0u, csControlIndex210 }, { 1u, 2u, 0u, csControlIndex120 }, { 0u, 2u, 1u, csControlIndex021 },
    { 0u, 1u, 2u, csControlIndex012 }, { 1u, 0u, 2u, csControlIndex102 }, { 2u, 0u, 1u, csControlIndex201 },
    { 1u, 1u, 1u, csControlIndex111 }
  }};
  // csPowerIndices[i][j] is the index of A_ij
  static constexpr std::array<std::array<uint32_t, 4u>, 4u> csPowerIndices = {{
    { csPowerIndex00, csPowerIndex01, csPowerIndex02, csPowerIndex03 },
    { csPowerIndex10, csPowerIndex11, csPowerIndex12, cInvalid      },
    { csPowerIndex20, csPowerIndex21, cInvalid,       cInvalid      },
    { csPowerIndex30, cInvalid,       cInvalid,       cInvalid      }
  }};

  for(auto &coefficient : mPowerCoefficients) {
    coefficient = Vector::Zero();
  }
  // Each term is 3! / (i! j! k!) * P_ijk * u^i * v^j * (1 - u - v)^k, and the last factor is expanded by the trinomial theorem.
  for(auto const &exponents : csControlExponents) {
    auto const i = exponents[0u];
    auto const j = exponents[1u];
    auto const k = exponents[2u];
    double const bernstein = csFactorials[3u] / (csFactorials[i] * csFactorials[j] * csFactorials[k]);
    for(uint32_t a = 0u; a <= k; ++a) {
      for(uint32_t b = 0u; a + b <= k; ++b) {
        double const sign = ((a + b) % 2u == 0u ? 1.0 : -1.0);
        double const trinomial = csFactorials[k] / (csFactorials[a] * csFactorials[b] * csFactorials[k - a - b]);
        mPowerCoefficients[csPowerIndices[i + a][j + b]] += mControlPoints[exponents[3u]] * (bernstein * trinomial * sign);
      }
    }
  }
}

template<typename tCoefficient>
void BezierTriangle::evaluatePower(std::array<tCoefficient, csPowerCoefficientsSize> const &aCoefficients, double const aU, double const aV,
                                   tCoefficient &aDistance, tCoefficient &aDerivativeU, tCoefficient &aDerivativeV) {
  auto const u2 = aU * aU;
  auto const v2 = aV * aV;
  auto const uv = aU * aV;
  auto const &a = aCoefficients;

  aDistance    = a[csPowerIndex00] + a[csPowerIndex10] * aU + a[csPowerIndex20] * u2 + a[csPowerIndex30] * (u2 * aU) +
                 a[csPowerIndex01] * aV + a[csPowerIndex11] * uv + a[csPowerIndex21] * (u2 * aV) +
                 a[csPowerIndex02] * v2 + a[csPowerIndex12] * (aU * v2) + a[csPowerIndex03] * (v2 * aV);
  aDerivativeU = a[csPowerIndex10] + a[csPowerIndex20] * (2.0 * aU) + a[csPowerIndex30] * (3.0 * u2) +
                 a[csPowerIndex11] * aV + a[csPowerIndex21] * (2.0 * uv) + a[csPowerIndex12] * v2;
  aDerivativeV = a[csPowerIndex01] + a[csPowerIndex11] * aU + a[csPowerIndex21] * u2 +
                 a[csPowerIndex02] * (2.0 * aV) + a[csPowerIndex12] * (2.0 * uv) + a[csPowerIndex03] * (3.0 * v2);
}

void BezierTriangle::setMissingFields4(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious) {
  mNeighbourDividerPlanes[1u] = Plane::createFrom1vector2points(mUnderlyingPlane.mNormal + aTriangleNext.mUnderlyingPlane.mNormal,
                                                                mControlPoints[csControlIndexOriginalVertex1],
                                                                mControlPoints[csControlIndexAboveOriginalCentroid]);
  mNeighbourDividerPlanes[2u] = Plane::createFrom1vector2points(mUnderlyingPlane.mNormal + aTrianglePrevious.mUnderlyingPlane.mNormal,
                                                                mControlPoints[csControlIndexOriginalVertex0],
                                                                mControlPoints[csControlIndexAboveOriginalCentroid]);
  mNeighbourDividerPlanes[1u].makeDistancePositive(mControlPoints[csControlIndexMiddle]);
  mNeighbourDividerPlanes[2u].makeDistancePositive(mControlPoints[csControlIndexMiddle]);
}

Vertex BezierTriangle::interpolateLinear(double const aBary0, double const aBary1, double const aBary2) const {
  return mControlPoints[csControlIndexOriginalVertex0]       * aBary0 +
         mControlPoints[csControlIndexOriginalVertex1]       * aBary1 +
         mControlPoints[csControlIndexAboveOriginalCentroid] * aBary2;
}

Vertex BezierTriangle::interpolate(double const aBary0, double const aBary1, double const aBary2) const {
  auto const bary0_2 = aBary0 * aBary0;
  auto const bary1_2 = aBary1 * aBary1;
  auto const bary2_2 = aBary2 * aBary2;

  return mControlPoints[csControlIndexOriginalVertex0]             * aBary0 * bary0_2 +
         mControlPoints[csControlIndexOriginalVertex1]             * aBary1 * bary1_2 +
         mControlPoints[csControlIndexAboveOriginalCentroid]       * aBary2 * bary2_2 +
         3.0 *
        (mControlPoints[csControlIndexOnOriginalSide0]             * aBary1 * bary0_2 +
         mControlPoints[csControlIndexOnOriginalSide1]             * aBary0 * bary1_2 +
         mControlPoints[csControlIndexOnSideToOriginalCentroid0]   * aBary2 * bary1_2 +
         mControlPoints[csControlIndexOnSideToOriginalCentroid1]   * aBary1 * bary2_2 +
         mControlPoints[csControlIndexOnSideFromOriginalCentroid0] * aBary0 * bary2_2 +
         mControlPoints[csControlIndexOnSideFromOriginalCentroid1] * aBary2 * bary0_2) +
         mControlPoints[csControlIndexMiddle]                      * aBary0 * aBary1 * aBary2 * 6.0;
}

BezierIntersection BezierTriangle::intersectByBezierEvaluations(Ray const &aRay, LimitPlaneIntersection const aShouldLimitPlaneIntersection) const {
  auto inPlane = mUnderlyingPlane.intersect(aRay);
  BezierIntersection result;
  if(inPlane.mValid && ::abs(inPlane.mDistance) > -mHeightInside && ::abs(inPlane.mDistance) > mHeightOutside) {  // Make sure we don't intersect the same triangle again
    Vector barycentric = mBarycentricInverse * inPlane.mPoint;
    if(aShouldLimitPlaneIntersection == LimitPlaneIntersection::cNone ||
       (barycentric(0) >= 0.0 && barycentric(0) <= 1.0 &&
        barycentric(1) >= 0.0 && barycentric(1) <= 1.0 &&
        barycentric(2) >= 0.0 && barycentric(2) <= 1.0)) {
      auto distanceInside = mHeightInside / inPlane.mCosIncidence;
      auto distanceOutside = mHeightOutside / inPlane.mCosIncidence;
      double closer = inPlane.mDistance + (inPlane.mCosIncidence > 0.0 ? distanceInside : distanceOutside);
      double further = inPlane.mDistance + (inPlane.mCosIncidence > 0.0 ? distanceOutside : distanceInside);

      Vertex pointOnRay = aRay.mStart + aRay.mDirection * closer;
      Vertex barycentric = mBarycentricInverse * mUnderlyingPlane.project(pointOnRay);
      auto diffCloser = ::abs(mUnderlyingPlane.distance(pointOnRay)) - ::abs(mUnderlyingPlane.distance(interpolate(barycentric)));

      pointOnRay = aRay.mStart + aRay.mDirection * further;
      barycentric = mBarycentricInverse * mUnderlyingPlane.project(pointOnRay);
      auto diffFurther = ::abs(mUnderlyingPlane.distance(pointOnRay)) - ::abs(mUnderlyingPlane.distance(interpolate(barycentric)));

      double middle;
      auto denom = diffCloser - diffFurther;
      if(::abs(denom) < csIntersectionEstimationEpsilon) {
        middle = (closer + further) / 2.0;
      }
      else {
        middle = (diffCloser * further - diffFurther * closer) / denom;
      }
      Vector projectionDirection = mUnderlyingPlane.mNormal;

      for(uint32_t i = 0u; i < csRootSearchIterations4bezierEvaluations; ++i) {
        result.mIntersection.mDistance = middle;
        pointOnRay = aRay.mStart + aRay.mDirection * middle;
        auto intersection = mUnderlyingPlane.intersect(pointOnRay, projectionDirection);
        result.mBarycentric = mBarycentricInverse * intersection.mPoint;
        result.mNormal = getNormal(result.mBarycentric);
        result.mIntersection.mPoint = interpolate(result.mBarycentric);
        projectionDirection = (result.mIntersection.mPoint - intersection.mPoint).normalized();
        middle = ((result.mIntersection.mPoint - aRay.mStart).dot(result.mNormal) / aRay.mDirection.dot(result.mNormal));
      }
      if(aRay.getDistance(result.mIntersection.mPoint) > csMaxIntersectionDistanceFromRay || result.mIntersection.mDistance < (further - closer) * csMinimalRayDistance) {
        result.mWhat = BezierIntersection::What::cNone;
      }
      else {
        classifyIntersection(result, aRay);
      }
    }
    else {
      result.mWhat = BezierIntersection::What::cNone;
    }
  }
  else {
    result.mWhat = BezierIntersection::What::cNone;
  }
  return result;
}

BezierIntersection BezierTriangle::intersectKajiya(Ray const &aRay, LimitPlaneIntersection const aShouldLimitPlaneIntersection) const {
  auto inPlane = mUnderlyingPlane.intersect(aRay);
  BezierIntersection result;
  if(inPlane.mValid && ::abs(inPlane.mDistance) > -mHeightInside && ::abs(inPlane.mDistance) > mHeightOutside) {  // Make sure we don't intersect the same triangle again
    Vector barycentric = mBarycentricInverse * inPlane.mPoint;
    if(aShouldLimitPlaneIntersection == LimitPlaneIntersection::cNone ||
       (barycentric(0) >= 0.0 && barycentric(0) <= 1.0 &&
        barycentric(1) >= 0.0 && barycentric(1) <= 1.0 &&
        barycentric(2) >= 0.0 && barycentric(2) <= 1.0)) {
      auto distanceInside = mHeightInside / inPlane.mCosIncidence;
      auto distanceOutside = mHeightOutside / inPlane.mCosIncidence;
      double closer = inPlane.mDistance + (inPlane.mCosIncidence > 0.0 ? distanceInside : distanceOutside);
      double further = inPlane.mDistance + (inPlane.mCosIncidence > 0.0 ? distanceOutside : distanceInside);

      // The ray is the intersection of these two perpendicular planes, both containing it.
      Vector const normal1 = util::getAperpendicular(aRay.mDirection);
      Vector const normal2 = aRay.mDirection.cross(normal1);
      std::array<double, csPowerCoefficientsSize> projected1;
      std::array<double, csPowerCoefficientsSize> projected2;
      for(uint32_t i = 0u; i < csPowerCoefficientsSize; ++i) {
        projected1[i] = mPowerCoefficients[i].dot(normal1);
        projected2[i] = mPowerCoefficients[i].dot(normal2);
      }
      projected1[csPowerIndex00] -= normal1.dot(aRay.mStart);    // Now projected* describe the signed distance of S(u, v) from the planes.
      projected2[csPowerIndex00] -= normal2.dot(aRay.mStart);

      double u = barycentric(0);
      double v = barycentric(1);
      bool converged = false;
      for(uint32_t i = 0u; i < csRootSearchIterations4kajiya; ++i) {
        double distance1, derivative1u, derivative1v;
        double distance2, derivative2u, derivative2v;
        evaluatePower(projected1, u, v, distance1, derivative1u, derivative1v);
        evaluatePower(projected2, u, v, distance2, derivative2u, derivative2v);
        auto const determinant = derivative1u * derivative2v - derivative1v * derivative2u;
        if(::abs(determinant) < csKajiyaDeterminantEpsilon) {    // Ray parallel to the surface, no reliable step.
          break;
        }
        else { // Nothing to do
        }
        auto const deltaU = (distance1 * derivative2v - distance2 * derivative1v) / determinant;
        auto const deltaV = (derivative1u * distance2 - derivative2u * distance1) / determinant;
        u -= deltaU;
        v -= deltaV;
        if(::abs(deltaU) + ::abs(deltaV) < csKajiyaBarycentricEpsilon) {
          converged = true;
          break;
        }
        else { // Nothing to do
        }
      }

      Vector derivativeU, derivativeV;
      evaluatePower(mPowerCoefficients, u, v, result.mIntersection.mPoint, derivativeU, derivativeV);
      result.mBarycentric = Vertex{u, v, 1.0 - u - v};
      result.mNormal = derivativeU.cross(derivativeV).normalized();
      result.mIntersection.mDistance = aRay.mDirection.dot(result.mIntersection.mPoint - aRay.mStart);
      if(!converged || result.mIntersection.mDistance < (further - closer) * csMinimalRayDistance) {
        result.mWhat = BezierIntersection::What::cNone;
      }
      else {
        classifyIntersection(result, aRay);
      }
    }
    else {
      result.mWhat = BezierIntersection::What::cNone;
    }
  }
  else {
    result.mWhat = BezierIntersection::What::cNone;
  }
  return result;
}

void BezierTriangle::classifyIntersection(BezierIntersection &aResult, Ray const &aRay) const {
  uint32_t outside = (mNeighbourDividerPlanes[0].distance(aResult.mIntersection.mPoint) < 0.0 ? 1u : 0u);
  outside |= (mNeighbourDividerPlanes[1].distance(aResult.mIntersection.mPoint) < 0.0 ? 2u : 0u);
  outside |= (mNeighbourDividerPlanes[2].distance(aResult.mIntersection.mPoint) < 0.0 ? 4u : 0u);
  if(outside == 1u) {
    aResult.mWhat = BezierIntersection::What::cFollowSide0;
  }
  else if(outside == 2u) {
    aResult.mWhat = BezierIntersection::What::cFollowSide1;
  }
  else if(outside == 4u) {
    aResult.mWhat = BezierIntersection::What::cFollowSide2;
  }
  else { // Most probably not possible for 2 sides at the same time. If yes, that rare case is not interesting
    aResult.mWhat = BezierIntersection::What::cIntersect;
    aResult.mIntersection.mCosIncidence = aRay.mDirection.dot(aResult.mNormal);
  }
}

// The point where the line through aPoint in mLiftDirection meets aPlane, or the orthogonal projection for a zero mLiftDirection.
Vertex BezierTriangle::liftOnto(Plane const &aPlane, Vertex const &aPoint) const {
  Vertex result;
  if(mLiftDirection.isZero()) {
    result = aPlane.project(aPoint);
  }
  else {
    auto const intersection = aPlane.intersect(aPoint, mLiftDirection);
    result = (intersection.mValid ? intersection.mPoint : aPlane.project(aPoint));  // Parallel only in degenerate cases.
  }
  return result;
}

Vector BezierTriangle::getNormal(Vector const &aBarycentric) const {
  double bary0_2 = aBarycentric(0) * aBarycentric(0);
  double bary1_2 = aBarycentric(1) * aBarycentric(1);
  double bary2_2 = aBarycentric(2) * aBarycentric(2);

  Vector component0 = mControlPoints[csControlIndex300] * bary0_2 +
                      mControlPoints[csControlIndex102] * bary2_2 +
                      mControlPoints[csControlIndex120] * bary1_2 +
                      2.0 *
                     (mControlPoints[csControlIndex201] * aBarycentric(0) * aBarycentric(2) +
                      mControlPoints[csControlIndex210] * aBarycentric(0) * aBarycentric(1) +
                      mControlPoints[csControlIndex111] * aBarycentric(2) * aBarycentric(1));

  Vector component1 = mControlPoints[csControlIndex030] * bary1_2 +
                      mControlPoints[csControlIndex012] * bary2_2 +
                      mControlPoints[csControlIndex210] * bary0_2 +
                      2.0 *
                     (mControlPoints[csControlIndex111] * aBarycentric(0) * aBarycentric(2) +
                      mControlPoints[csControlIndex120] * aBarycentric(0) * aBarycentric(1) +
                      mControlPoints[csControlIndex021] * aBarycentric(1) * aBarycentric(2));

  Vector component2 = mControlPoints[csControlIndex003] * bary2_2 +
                      mControlPoints[csControlIndex201] * bary0_2 +
                      mControlPoints[csControlIndex021] * bary1_2 +
                      2.0 *
                     (mControlPoints[csControlIndex102] * aBarycentric(0) * aBarycentric(2) +
                      mControlPoints[csControlIndex012] * aBarycentric(1) * aBarycentric(2) +
                      mControlPoints[csControlIndex111] * aBarycentric(0) * aBarycentric(1));

  Vector componentA = mBezierDerivativeDirectionVectorA(0) * component0 +
                      mBezierDerivativeDirectionVectorA(1) * component1 +
                      mBezierDerivativeDirectionVectorA(2) * component2;
  Vector componentB = mBezierDerivativeDirectionVectorB(0) * component0 +
                      mBezierDerivativeDirectionVectorB(1) * component1 +
                      mBezierDerivativeDirectionVectorB(2) * component2;
  return componentA.cross(componentB).normalized();
}
