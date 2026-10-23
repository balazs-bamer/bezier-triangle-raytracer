#include "hostUtil.h"

UniformHemisphere::UniformHemisphere(uint32_t const aBelts)
  : cmBeltWidth(cgPi / 2.0 / aBelts)
  , mUniform1(0.0, 1.0)
  , mUniform2pi(0.0, cgPi * 2.0) {
  mBelts.reserve(aBelts);
  mPatchCount= 0u;
  for(uint32_t i = 0u; i < aBelts; ++i) {
    uint32_t now = ::ceil(4.0 * aBelts * ::sin((2.0 * i + 1.0) / (4.0 * aBelts) * cgPi));
    mBelts.emplace_back(cgPi * 2.0 / now, mPatchCount);
    mPatchCount += now;
  }
}

std::pair<Vector, uint32_t> UniformHemisphere::getRandom() {
  Vector direction;
  uint32_t index;
  double incidence = ::acos(mUniform1(mRandomGenerator));
  double beltRadius = ::sin(incidence);
  auto turn = mUniform2pi(mRandomGenerator);
  direction(0) = ::cos(incidence);
  direction(1) = beltRadius * ::cos(turn);
  direction(2) = beltRadius * ::sin(turn);
  direction.normalize();
  auto [patchWidth, soFar] = mBelts[static_cast<uint32_t>(incidence / cmBeltWidth)];
  index = soFar + static_cast<uint32_t>(turn / patchWidth);
  return std::pair(direction, index);
}
