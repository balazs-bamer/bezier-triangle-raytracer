#include "bezierLens.h"


std::pair<Ray, RefractionResult> BezierLens::refract(Ray const &aRay, RefractionResult const aExpected) const {
  Ray result;
  RefractionResult statusLater;
  auto const intersect = mMesh.intersect(aRay);
  if(intersect.mWhat == BezierIntersection::What::cIntersect) {
    statusLater = (intersect.mIntersection.mCosIncidence < 0.0 ? RefractionResult::cInside : RefractionResult::cOutside);
    result.mStart = intersect.mIntersection.mPoint;
    auto const effectiveRefractiveFactor = (statusLater == RefractionResult::cInside ? 1.0 / mRefractiveIndex : mRefractiveIndex);
    auto const sin2refraction = effectiveRefractiveFactor * effectiveRefractiveFactor * (1.0 - intersect.mIntersection.mCosIncidence * intersect.mIntersection.mCosIncidence);
    if(sin2refraction < csMaxSin2refraction) {
      if(sin2refraction > csMinSin2refraction) {
        auto const directionFactor = (statusLater == RefractionResult::cInside ? 1.0 : -1.0);
        auto const normal = intersect.mNormal * directionFactor;
        auto cos1 = ::abs(intersect.mIntersection.mCosIncidence);
        auto cos2 = ::sqrt(1.0 - sin2refraction);
        result.mDirection = (aRay.mDirection * effectiveRefractiveFactor + normal * (effectiveRefractiveFactor * cos1 - cos2)).normalized();
      }
      else {
        result.mDirection = aRay.mDirection;        // Too little incidence, continue in the same direction.
      }
    }
    else {
      statusLater = RefractionResult::cNone;        // Total reflection, or close to it. Ignored.
    }
  }
  else {
    statusLater = RefractionResult::cNone;
  }
  statusLater = (statusLater == aExpected ? statusLater : RefractionResult::cNone);
  return std::make_pair(result, statusLater);
}

