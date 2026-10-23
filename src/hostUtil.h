#ifndef BEZIER_TRIANGLE_RAYTRACER_HOSTUTIL
#define BEZIER_TRIANGLE_RAYTRACER_HOSTUTIL

#include "3dGeomUtil.h"
#include <vector>
#include <random>

class UniformHemisphere final {
private:
  double const                             cmBeltWidth;
  std::vector<std::pair<double, uint32_t>>  mBelts;
  std::ranlux24_base                       mRandomGenerator;
  std::uniform_real_distribution<double>    mUniform1;
  std::uniform_real_distribution<double>    mUniform2pi;
  uint32_t                                 mPatchCount;

public:
  UniformHemisphere(uint32_t const aBelts);

  uint32_t getPatchCount() const { return mPatchCount; }

  // Returns unit vector and the patch index.
  std::pair<Vector, uint32_t> getRandom();
};

#endif
