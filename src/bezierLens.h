#ifndef BEZIER_TRIANGLE_RAYTRACER_BEZIERLENS
#define BEZIER_TRIANGLE_RAYTRACER_BEZIERLENS

#include "bezierMesh.h"


enum class RefractionResult : uint32_t {
  cNone    = 0u,
  cInside  = 1u,
  cOutside = 2u
};

// This class first will be used only to simulate discrete lenses. Compounds will come later, if any.
// The lens is bounded either by one closed surface, or by two open ones in the order of the light path. In both cases the surface normals
// must point out of the lens. The ray must enter through the first surface and leave through the last one (the same for a closed lens).
// A ray not hitting the very next surface first is terminated, as if it would hit the black tube holding the lenses.
class BezierLens final {
private:
  static constexpr double csMaxSin2refraction = 0.99;  // Squared, approximately 82 degrees.
  static constexpr double csMinSin2refraction = 1e-12; // Don't calculate refraction below this, let the ray pass as it came.

  double                  mRefractiveIndex;
  std::vector<BezierMesh> mSurfaces;

public:
  BezierLens(double const aRi, BezierMesh const & aMesh) : mRefractiveIndex(aRi), mSurfaces({aMesh}) {}
  BezierLens(double const aRi, BezierMesh&& aMesh)      : mRefractiveIndex(aRi) { mSurfaces.push_back(std::move(aMesh)); }
  BezierLens(double const aRi, std::vector<BezierMesh> const &aSurfaces) : mRefractiveIndex(aRi), mSurfaces(aSurfaces) { checkSurfaceCount(); }
  BezierLens(double const aRi, std::vector<BezierMesh>&& aSurfaces)     : mRefractiveIndex(aRi), mSurfaces(std::move(aSurfaces)) { checkSurfaceCount(); }

  // The ray either is refracted, or gets lost (if misses the very next surface or suffers total reflection).
  // aExpected is cInside when the ray should enter the lens, cOutside when it should leave it.
  std::pair<Ray, RefractionResult> refract(Ray const &aRay, RefractionResult const aExpected) const;

private:
  void checkSurfaceCount() const {
    if(mSurfaces.empty() || mSurfaces.size() > 2u) {
      throw "A lens needs one closed or two open surfaces.";
    }
    else { // Nothing to do
    }
  }

  std::pair<BezierIntersection, uint32_t> getClosestIntersection(Ray const &aRay) const;  // Among all the surfaces, with the index of the surface.
};

#endif // BEZIERLENS_H
