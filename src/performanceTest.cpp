#include "bezierMesh.h"

#include<random>
#include<vector>
#include<chrono>
#include<string>
#include<cstdlib>
#include<iomanip>
#include<iostream>


// Performance comparison of BezierTriangle::intersectByBezierEvaluations and BezierTriangle::intersectKajiya.
// Usage: performance-test [calls per measurement]

constexpr uint32_t cgDefaultCallsPerMeasurement = 1000000u;
constexpr uint32_t cgRaysPerMeasurement         = 20000u;   // Distinct (triangle, ray) pairs, repeated until the call count is reached.
constexpr uint32_t cgRandomSeed                 = 20261023u;

struct TestMesh final {
  std::string mName;
  BezierMesh  mBezier;
  double      mRadius;    // Bounding radius around the origin.
};

struct RayCase final {
  uint32_t mTriangle;
  Ray      mRay;
};

enum class RayKind : uint32_t {
  cPerpendicular = 0u,  // From outside towards a surface point, close to the triangle normal.
  cOblique       = 1u,  // From outside towards a surface point, 30-70 degrees from the normal.
  cGrazing       = 2u,  // From outside towards a surface point, 75-88 degrees from the normal.
  cFromInside    = 3u,  // From the origin (inside the solid) towards a surface point.
  cNearMiss      = 4u,  // Towards a point just outside the triangle, plane limit on: mostly early-out.
  cNeighbour     = 5u,  // Same as cNearMiss, but without plane limit, as BezierMesh does when following a neighbour.
  cRandom        = 6u   // Random lines around the object: mostly early-out.
};

char const * const cgRayKindNames[] = { "perpendicular", "oblique", "grazing", "from inside", "near miss", "neighbour", "random" };

TestMesh makeTestMesh(std::string const &aName, int32_t const aSectors, int32_t const aBelts, std::function<double(double)> aEnvelope, Vector const &aSize) {
  Mesh mesh;
  mesh.makeSolidOfRevolution(aSectors, aBelts, aEnvelope, aSize);
  mesh.standardizeVertices();
  mesh.standardizeNormals();
  double radius = 0.0;
  for(auto const &triangle : mesh) {
    for(auto const &vertex : triangle) {
      radius = std::max(radius, vertex.norm());
    }
  }
  return TestMesh{aName, BezierMesh(mesh), radius};
}

// Direction making aTilt angle with aNormal, at a random azimuth.
Vector getTilted(Vector const &aNormal, double const aTilt, double const aAzimuth) {
  Vector const perpendicular0 = util::getAperpendicular(aNormal);
  Vector const perpendicular1 = aNormal.cross(perpendicular0);
  return aNormal * ::cos(aTilt) + (perpendicular0 * ::cos(aAzimuth) + perpendicular1 * ::sin(aAzimuth)) * ::sin(aTilt);
}

std::vector<RayCase> makeRays(TestMesh const &aMesh, RayKind const aKind, std::mt19937 &aGenerator) {
  std::uniform_int_distribution<uint32_t> triangleDistribution(0u, aMesh.mBezier.size() - 1u);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::uniform_real_distribution<double> azimuthDistribution(0.0, 2.0 * cgPi);
  double const startDistance = aMesh.mRadius * 3.0;

  std::vector<RayCase> result;
  result.reserve(cgRaysPerMeasurement);
  for(uint32_t i = 0u; i < cgRaysPerMeasurement; ++i) {
    auto const index = triangleDistribution(aGenerator);
    auto const &triangle = aMesh.mBezier[index];
    Vector const normal = util::getNormal(triangle.getControlPoint(0u), triangle.getControlPoint(1u), triangle.getControlPoint(2u)).normalized();

    // Uniform barycentric coordinates inside the triangle.
    double bary0 = unit(aGenerator);
    double bary1 = unit(aGenerator);
    if(bary0 + bary1 > 1.0) {
      bary0 = 1.0 - bary0;
      bary1 = 1.0 - bary1;
    }
    else { // Nothing to do
    }
    double bary2 = 1.0 - bary0 - bary1;
    if(aKind == RayKind::cNearMiss || aKind == RayKind::cNeighbour) {  // Just beyond the edge opposite to vertex 2.
      bary2 = -0.01 - 0.3 * unit(aGenerator);
      bary0 = (1.0 - bary2) * unit(aGenerator);
      bary1 = 1.0 - bary2 - bary0;
    }
    else { // Nothing to do
    }
    Vertex const target = triangle.interpolate(bary0, bary1, bary2);

    Ray ray;
    switch(aKind) {
    case RayKind::cPerpendicular:
      ray = Ray(Vertex::Zero(), -getTilted(normal, 0.1 * unit(aGenerator), azimuthDistribution(aGenerator)));
      ray.mStart = target - ray.mDirection * startDistance;
      break;
    case RayKind::cOblique:
      ray = Ray(Vertex::Zero(), -getTilted(normal, (30.0 + 40.0 * unit(aGenerator)) / 180.0 * cgPi, azimuthDistribution(aGenerator)));
      ray.mStart = target - ray.mDirection * startDistance;
      break;
    case RayKind::cGrazing:
      ray = Ray(Vertex::Zero(), -getTilted(normal, (75.0 + 13.0 * unit(aGenerator)) / 180.0 * cgPi, azimuthDistribution(aGenerator)));
      ray.mStart = target - ray.mDirection * startDistance;
      break;
    case RayKind::cFromInside:
      ray = Ray(Vertex::Zero(), target);
      break;
    case RayKind::cNearMiss:
    case RayKind::cNeighbour:
      ray = Ray(Vertex::Zero(), -getTilted(normal, 0.6 * unit(aGenerator), azimuthDistribution(aGenerator)));
      ray.mStart = target - ray.mDirection * startDistance;
      break;
    case RayKind::cRandom:
    default: {
        Vector const towards{unit(aGenerator) - 0.5, unit(aGenerator) - 0.5, unit(aGenerator) - 0.5};
        Vector const from = getTilted(Vector{0.0, 0.0, 1.0}, ::acos(1.0 - 2.0 * unit(aGenerator)), azimuthDistribution(aGenerator));
        ray = Ray(from * startDistance, towards * aMesh.mRadius - from * startDistance);
      }
      break;
    }
    result.push_back(RayCase{index, ray});
  }
  return result;
}

struct Measurement final {
  double   mNanosecondsPerCall;
  uint32_t mHits;         // Over one pass through the ray cases.
  double   mChecksum;     // Prevents the compiler from optimizing the calls away.
};

template<typename tMethod>
Measurement measure(TestMesh const &aMesh, std::vector<RayCase> const &aRays, BezierTriangle::LimitPlaneIntersection const aLimit, uint32_t const aCalls, tMethod &&aMethod) {
  Measurement result{0.0, 0u, 0.0};
  for(auto const &rayCase : aRays) {     // Warm-up pass, also counts the hits.
    auto const intersection = aMethod(aMesh.mBezier[rayCase.mTriangle], rayCase.mRay, aLimit);
    if(intersection.mWhat == BezierIntersection::What::cIntersect) {
      ++result.mHits;
    }
    else { // Nothing to do
    }
  }

  uint32_t const passes = (aCalls + aRays.size() - 1u) / aRays.size();
  auto const start = std::chrono::high_resolution_clock::now();
  for(uint32_t pass = 0u; pass < passes; ++pass) {
    for(auto const &rayCase : aRays) {
      auto const intersection = aMethod(aMesh.mBezier[rayCase.mTriangle], rayCase.mRay, aLimit);
      result.mChecksum += static_cast<double>(intersection.mWhat) + (intersection.mWhat == BezierIntersection::What::cIntersect ? intersection.mIntersection.mDistance : 0.0);
    }
  }
  auto const end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::nano> const diff = end - start;
  result.mNanosecondsPerCall = diff.count() / (static_cast<double>(passes) * aRays.size());
  return result;
}

int main(int argc, char **argv) {
  uint32_t const calls = (argc > 1 ? static_cast<uint32_t>(std::strtoul(argv[1], nullptr, 10)) : cgDefaultCallsPerMeasurement);

  auto ellipse = [](double const aX){ return ::sqrt(1.0 - aX * aX); };
  auto bulging = [](double const aX){    // Has concave parts, as in test.cpp
    auto x2 = aX * aX;
    return ::sqrt(1.0 - x2) + 0.7 * (::exp(-4.0) - ::exp(-4.0 * x2));
  };
  std::vector<TestMesh> meshes;
  meshes.push_back(makeTestMesh("sphere 7x3 coarse",    7,  3, ellipse, Vector{1.0, 1.0, 1.0}));
  meshes.push_back(makeTestMesh("sphere 21x15 fine",   21, 15, ellipse, Vector{1.0, 1.0, 1.0}));
  meshes.push_back(makeTestMesh("ellipsoid 7x3 coarse", 7,  3, ellipse, Vector{1.0, 4.0, 2.0}));
  meshes.push_back(makeTestMesh("ellipsoid 15x5",      15,  5, ellipse, Vector{1.0, 4.0, 2.0}));
  meshes.push_back(makeTestMesh("bulging 21x15",       21, 15, bulging, Vector{1.0, 4.0, 2.0}));

  auto byBezierEvaluations = [](BezierTriangle const &aTriangle, Ray const &aRay, BezierTriangle::LimitPlaneIntersection const aLimit) {
    return aTriangle.intersectByBezierEvaluations(aRay, aLimit);
  };
  auto kajiya = [](BezierTriangle const &aTriangle, Ray const &aRay, BezierTriangle::LimitPlaneIntersection const aLimit) {
    return aTriangle.intersectKajiya(aRay, aLimit);
  };

  std::cout << "Calls per measurement: " << calls << ", distinct rays per measurement: " << cgRaysPerMeasurement << "\n\n";
  std::cout << std::left << std::setw(22) << "mesh" << std::setw(15) << "rays"
            << std::right << std::setw(12) << "old ns" << std::setw(12) << "Kajiya ns" << std::setw(10) << "speedup"
            << std::setw(12) << "old hits" << std::setw(12) << "Kajiya hits" << '\n';

  std::mt19937 generator(cgRandomSeed);
  double checksum = 0.0;
  double totalOld = 0.0;
  double totalKajiya = 0.0;
  for(auto const &mesh : meshes) {
    for(uint32_t kind = static_cast<uint32_t>(RayKind::cPerpendicular); kind <= static_cast<uint32_t>(RayKind::cRandom); ++kind) {
      auto const rayKind = static_cast<RayKind>(kind);
      auto const rays = makeRays(mesh, rayKind, generator);
      auto const limit = (rayKind == RayKind::cNeighbour ? BezierTriangle::LimitPlaneIntersection::cNone : BezierTriangle::LimitPlaneIntersection::cThis);
      auto const old = measure(mesh, rays, limit, calls, byBezierEvaluations);
      auto const fast = measure(mesh, rays, limit, calls, kajiya);
      checksum += old.mChecksum + fast.mChecksum;
      totalOld += old.mNanosecondsPerCall;
      totalKajiya += fast.mNanosecondsPerCall;
      std::cout << std::left << std::setw(22) << mesh.mName << std::setw(15) << cgRayKindNames[kind]
                << std::right << std::fixed << std::setprecision(1)
                << std::setw(12) << old.mNanosecondsPerCall << std::setw(12) << fast.mNanosecondsPerCall
                << std::setw(9) << old.mNanosecondsPerCall / fast.mNanosecondsPerCall << 'x'
                << std::setw(12) << old.mHits << std::setw(12) << fast.mHits << '\n';
    }
  }
  std::cout << "\nSum of per-call times: old " << totalOld << " ns, Kajiya " << totalKajiya << " ns, speedup " << totalOld / totalKajiya << "x\n";
  std::cout << "Checksum (ignore): " << std::setprecision(3) << checksum << '\n';
  return 0;
}
