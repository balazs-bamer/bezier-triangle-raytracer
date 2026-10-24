#include "bezierMesh.h"


BezierMesh::BezierMesh(Mesh const &aMesh) : BezierMesh(aMesh, std::vector<bool>(aMesh.size(), true)) {}

BezierMesh::BezierMesh(Mesh const &aMesh, std::vector<bool> const &aMask)
  : mOriginalTriangles(aMesh.getMesh())
  , mOriginalNeighbours(aMesh.getFace2neighbours())
  , mFace2bezierBase(aMesh.size(), Mesh::csNoNeighbour) {
  for(uint32_t indexFace = 0u; indexFace < aMesh.size(); ++indexFace) {
    if(aMask[indexFace] && aMesh.isCoreFace(indexFace)) {          // A core triangle always has all its neighbours.
      mFace2bezierBase[indexFace] = mBuiltFaces.size() * 3u;
      mBuiltFaces.push_back(indexFace);
    }
    else { // Nothing to do
    }
  }
  mMesh.clear();
  mMesh.reserve(mBuiltFaces.size() * 3u);
  std::vector<std::pair<uint32_t, Plane>> rimPlanes;               // Bezier index and plane between the original neighbours
  auto const &vertex2averageNormals = aMesh.getVertex2averageNormals();
  for(auto const indexFace : mBuiltFaces) {
    auto const &neigh = mOriginalNeighbours[indexFace];
    auto const &originalTriangle = aMesh[indexFace];
    auto const originalCentroid = (originalTriangle[0u] + originalTriangle[1u] + originalTriangle[2u]) / 3.0;
    auto normal = util::getNormal(originalTriangle).normalized();
    for(uint32_t indexVertex = 0u; indexVertex < 3u; ++indexVertex) { // Appending a regular triangle means doing Clough-Tocher split on it and appending each small one.
      auto const &originalCommonVertex0 = originalTriangle[indexVertex];
      auto const &originalCommonVertex1 = originalTriangle[(indexVertex + 1u) % 3u];
      auto const &averageNormal0 = vertex2averageNormals.at(originalCommonVertex0);
      auto const &averageNormal1 = vertex2averageNormals.at(originalCommonVertex1);
                                                                                 // Neighbour of edge (index, index + 1)
      Plane planeBetweenOriginalNeighbours = Plane::createFrom1vector2points(normal + util::getNormal(aMesh[neigh.mFellowTriangles[indexVertex]]).normalized(),
                                                                             originalCommonVertex0, originalCommonVertex1);
      auto currentBase = mFace2bezierBase[indexFace];
      auto const neighbourBase = mFace2bezierBase[neigh.mFellowTriangles[indexVertex]];
      std::array<uint32_t, 3u> newNeighbourIndices({ (neighbourBase == Mesh::csNoNeighbour ? Mesh::csNoNeighbour : neighbourBase + neigh.mFellowCommonSideStarts[indexVertex]),
                                                     currentBase + (indexVertex + 1u) % 3u,
                                                     currentBase + (indexVertex + 2u) % 3u });
      if(neighbourBase == Mesh::csNoNeighbour) {
        Plane rimPlane = planeBetweenOriginalNeighbours;
        rimPlane.makeDistancePositive(originalCentroid);
        rimPlanes.emplace_back(std::make_pair(mMesh.size(), rimPlane));
      }
      else { // Nothing to do
      }
      mMesh.emplace_back(BezierTriangle(originalCommonVertex0, originalCommonVertex1, originalCentroid,
                                               averageNormal0, averageNormal1, planeBetweenOriginalNeighbours,
                                               newNeighbourIndices));
    }
  }
  setMissingFields(aMesh, &BezierTriangle::setMissingFields1);
  setMissingFields(aMesh, &BezierTriangle::setMissingFields2);
  setMissingFields(aMesh, &BezierTriangle::setMissingFields3);

  mBezier2rim.assign(mMesh.size(), Mesh::csNoNeighbour);
  for(auto const &[indexBezier, rimPlane] : rimPlanes) {           // Now the control point above the original centroid is known.
    auto const &bezier = mMesh[indexBezier];
    Vertex const vertex0 = bezier.getControlPoint(0u);
    Vertex const vertex1 = bezier.getControlPoint(1u);
    Vertex const vertex2 = bezier.interpolateAboveOriginalCentroid();
    mBezier2rim[indexBezier] = mRims.size();
    mRims.push_back(Rim{ Plane::createFrom3points(vertex0, vertex1, vertex2), util::getBarycentricInverse(vertex0, vertex1, vertex2), rimPlane });
  }
}

void BezierMesh::setMissingFields(Mesh const &aMesh, MissingFieldsMethod aMissingFieldsMethod) {
  Vertex originalCentroid;
  for(uint32_t i = 0u; i < mMesh.size(); ++i) {
    auto subTriangleIndex = i % 3u;
    auto indexBase = i - subTriangleIndex;
    auto const &triangleNext = mMesh[indexBase + (subTriangleIndex + 1u) % 3u];
    auto const &trianglePrevious = mMesh[indexBase + (subTriangleIndex + 2u) % 3u];
    if(subTriangleIndex % 3u == 0u) {
      auto const &originalTriangle = aMesh[mBuiltFaces[indexBase / 3u]];
      originalCentroid = (originalTriangle[0u] + originalTriangle[1u] + originalTriangle[2u]) / 3.0;
    }
    else { // Nothing to do
    }
    aMissingFieldsMethod(mMesh[i], originalCentroid, triangleNext, trianglePrevious);
  }
}

Mesh BezierMesh::interpolate(int32_t const aDivisor) const {
  Mesh result;
  Triangle barycentric{Vertex{{1.0, 0.0, 0.0}}, Vertex{{0.0, 1.0, 0.0}}, Vertex{{0.0, 0.0, 1.0}}};
  util::divide(barycentric, aDivisor, [&result, this](Triangle && aNewBary){
    for(auto const &bezier : mMesh) {
      result.push_back({ bezier.interpolate(aNewBary[0u]),
                         bezier.interpolate(aNewBary[1u]),
                         bezier.interpolate(aNewBary[2u]) });
    }
  } );
  return result;
}

std::vector<Vertex> BezierMesh::dumpControlPoints() const {
  std::vector<Vertex> result;
  result.reserve(mMesh.size() * BezierTriangle::csControlPointsSize);
  for(auto const &bezier : mMesh) {
    for(uint32_t i = 0u; i < BezierTriangle::csControlPointsSize; ++i) {
      result.push_back(bezier.getControlPoint(i));
    }
  }
  return result;
}

Mesh BezierMesh::splitThickBezierTriangles() const {
  static constexpr uint8_t csSplitMask[3u]  = { 1u, 2u, 4u };
  static constexpr uint8_t csSplitAll       = 7u;
  static constexpr uint8_t csSplitCount[8u] = { 1u, 2u, 2u, 3u, 2u, 3u, 3u, 4u };
  std::vector<uint8_t> splitSides(mOriginalNeighbours.size(), 0u);

  for(uint32_t indexOriginal = 0u; indexOriginal < mOriginalNeighbours.size(); ++indexOriginal) {
    auto indexSplit = mFace2bezierBase[indexOriginal];
    if(indexSplit != Mesh::csNoNeighbour) {                         // Only built triangles can be thick.
      Triangle original { mMesh[indexSplit].getControlPoint(0u),
                          mMesh[indexSplit + 1u].getControlPoint(0u),
                          mMesh[indexSplit + 2u].getControlPoint(0u) };
      auto plane = Plane::createFromTriangle(original);
      double max = ::abs(plane.distance(mMesh[indexSplit].interpolateAboveOriginalCentroid()));
      for(uint32_t i = 0; i < 3u; ++i) {
        for(auto const ratio : csSampleRatiosOriginalSide) {
          max = std::max(max, ::abs(plane.distance(mMesh[indexSplit + i].interpolate(ratio, 1.0 - ratio, 0.0))));
        }
      }
      if(max / util::getPerimeter(original) > csBezierHeightPerPerimeterLimit) {
        splitSides[indexOriginal] = csSplitAll;
        auto const &neigh = mOriginalNeighbours[indexOriginal];
        for(uint32_t side = 0u; side < 3u; ++side) {
          if(neigh.mFellowTriangles[side] != Mesh::csNoNeighbour) {
            splitSides[neigh.mFellowTriangles[side]] |= csSplitMask[neigh.mFellowCommonSideStarts[side]];
          }
          else { // Nothing to do
          }
        }
      }
      else { // Nothing to do
      }
    }
    else { // Nothing to do
    }
  }
  // TODO go over original edges to see if that edge is not split but other edges of the 2 containing
  // original triangles are. In this case this edge should be split as well to prevent it become a depressed edge.
  // Or perhaps inspect each remaining original edge if the convexity changes there.
  uint32_t newCount = 0u;
  newCount = std::accumulate(splitSides.cbegin(), splitSides.cend(), newCount, [](uint32_t aSofar, uint8_t aNew){ return aSofar + csSplitCount[aNew]; });
  Mesh result;
  result.reserve(newCount);
  for(uint32_t indexOriginal = 0u; indexOriginal < mOriginalNeighbours.size(); ++indexOriginal) {
    auto const &original = mOriginalTriangles[indexOriginal];
    uint8_t split = splitSides[indexOriginal];
    newCount = csSplitCount[split];
    if(newCount == 1u) {
      result.push_back(original);
    }
    else {
      std::array<Vertex, 3u> splitVertices;                         // splitVertices[i] is on side (i, i + 1) if that one is split
      for(uint32_t side = 0u; side < 3u; ++side) {
        if((split & csSplitMask[side]) != 0u) {
          splitVertices[side] = getSplitVertex(indexOriginal, side);
        }
        else { // Nothing to do
        }
      }
      if(newCount == 2u) {
        append2split(result, original, split, splitVertices);
      }
      else if(newCount == 3u) {
        append3split(result, original, split, splitVertices);
      }
      else { // newCount == 4u
        append4split(result, original, splitVertices);
      }
    }
  }
  return result;
}

// The split vertex of a side not built is taken from the built neighbour, because the Bezier curve of the common side is the same.
Vertex BezierMesh::getSplitVertex(uint32_t const aIndexOriginal, uint32_t const aSide) const {
  auto const indexBase = mFace2bezierBase[aIndexOriginal];
  if(indexBase != Mesh::csNoNeighbour) {
    return interpolate(indexBase + aSide, 0.5, 0.5, 0.0);
  }
  else {
    auto const &neigh = mOriginalNeighbours[aIndexOriginal];
    return interpolate(mFace2bezierBase[neigh.mFellowTriangles[aSide]] + neigh.mFellowCommonSideStarts[aSide], 0.5, 0.5, 0.0);
  }
}

//       1
//      /|\
//     / | \
//    /  |  \
//   /   |   \
//  /____|____\
// 0           2   index2 == 2
//
void BezierMesh::append2split(Mesh &aResult, Triangle const &aOriginalTriangle, uint8_t const aSplit, std::array<Vertex, 3u> const &aSplitVertices) {
  static constexpr uint8_t csIndexFor2onSide[8u] = { 3u, 0u, 1u, 3u, 2u, 3u, 3u, 3u };
  uint32_t index2 = csIndexFor2onSide[aSplit];
  auto const &splitVertex = aSplitVertices[index2];
  uint32_t indexAfter2 = (index2 + 1u) % 3u;
  uint32_t indexBefore2 = (index2 + 2u) % 3u;
  aResult.push_back(Triangle{ aOriginalTriangle[indexAfter2], aOriginalTriangle[indexBefore2], splitVertex });
  aResult.push_back(Triangle{ aOriginalTriangle[indexBefore2], aOriginalTriangle[index2], splitVertex });
}

//       1
//      /|\
//     / | \SA1
//    /  | /\
//   /   |/  \     or not 1-SB1 but 0-SA1, whichever is shorter
//  /____/____\
// 0    SB1    2   index1 == 0
//
void BezierMesh::append3split(Mesh &aResult, Triangle const &aOriginalTriangle, uint8_t const aSplit, std::array<Vertex, 3u> const &aSplitVertices) {
  static constexpr uint8_t csIndexFor1onSide[8u] = { 3u, 3u, 3u, 2u, 3u, 1u, 0u, 3u };
  uint32_t index1 = csIndexFor1onSide[aSplit];
  uint32_t indexAfter1 = (index1 + 1u) % 3u;
  uint32_t indexBefore1 = (index1 + 2u) % 3u;
  auto const &splitVertexBefore1 = aSplitVertices[indexBefore1];
  auto const &splitVertexAfter1  = aSplitVertices[indexAfter1];
  aResult.push_back(Triangle{ aOriginalTriangle[indexBefore1], splitVertexBefore1, splitVertexAfter1 });
  if((aOriginalTriangle[indexAfter1] - splitVertexBefore1).norm() < (aOriginalTriangle[index1] - splitVertexAfter1).norm()) {
    aResult.push_back(Triangle{ aOriginalTriangle[indexAfter1], splitVertexAfter1, splitVertexBefore1 });
    aResult.push_back(Triangle{ aOriginalTriangle[index1], aOriginalTriangle[indexAfter1], splitVertexBefore1 });
  }
  else {
    aResult.push_back(Triangle{ aOriginalTriangle[indexAfter1], splitVertexAfter1, aOriginalTriangle[index1] });
    aResult.push_back(Triangle{ aOriginalTriangle[index1], splitVertexAfter1, splitVertexBefore1 });
  }
}

//        1
//       /\
//      /  \
//   S0/____\S1
//    /\    /\
//   /  \  /  \
//  /____\/____\
// 0     S2     2
//
void BezierMesh::append4split(Mesh &aResult, Triangle const &aOriginalTriangle, std::array<Vertex, 3u> const &aSplitVertices) {
  Triangle const &middle = aSplitVertices;
  aResult.push_back(middle);
  for(uint32_t i = 0u; i < 3u; ++i) {
    aResult.push_back(Triangle{ aOriginalTriangle[i], middle[i], middle[(i + 2u) % 3u] });
  }
}

Vertex BezierMesh::interpolate(uint32_t const aIndex, double const aBary0, double const aBary1, double const aBary2) const {
  auto const &triangle = mMesh[aIndex];
  return triangle.interpolate(aBary0, aBary1, aBary2) * csSplitBezierInterpolateFactor +
         triangle.interpolateLinear(aBary0, aBary1, aBary2) * (1.0 - csSplitBezierInterpolateFactor);
}

BezierIntersection BezierMesh::intersect(Ray const &aRay) const {
  // Brute-force for now
  BezierIntersection result;
  result.mIntersection.mDistance = std::numeric_limits<double>::max();
  result.mWhat = BezierIntersection::What::cNone;
  for(uint32_t indexBezier = 0u; indexBezier < mMesh.size(); ++indexBezier) {
    auto const &bezier = mMesh[indexBezier];
    auto candidate = bezier.intersectKajiya(aRay, BezierTriangle::LimitPlaneIntersection::cThis);
    if(candidate.mWhat == BezierIntersection::What::cNone && mBezier2rim[indexBezier] != Mesh::csNoNeighbour) {
      candidate = intersectBeyondRim(mRims[mBezier2rim[indexBezier]], bezier, aRay);
    }
    else { // Nothing to do
    }
    if(candidate.mWhat == BezierIntersection::What::cFollowSide0 ||
       candidate.mWhat == BezierIntersection::What::cFollowSide1 ||
       candidate.mWhat == BezierIntersection::What::cFollowSide2) {
      auto const indexNeighbour = bezier.getNeighbours()[static_cast<uint32_t>(candidate.mWhat)];
      if(indexNeighbour == Mesh::csNoNeighbour) {                   // Beyond the rim of an open mesh.
        candidate.mWhat = BezierIntersection::What::cNone;
      }
      else {
        candidate = mMesh[indexNeighbour].intersectKajiya(aRay, BezierTriangle::LimitPlaneIntersection::cNone);
      }
    }
    else { // Nothing to do
    }
    if(candidate.mWhat == BezierIntersection::What::cIntersect && candidate.mBarycentric.minCoeff() < csMinBarycentricOfHit) {
      candidate.mWhat = BezierIntersection::What::cNone;            // On the polynomial continuation of the patch, see csMinBarycentricOfHit.
    }
    else { // Nothing to do
    }
    if(candidate.mWhat == BezierIntersection::What::cIntersect && candidate.mIntersection.mDistance < result.mIntersection.mDistance) {
      result = candidate;
    }
    else { // Nothing to do
    }
  }
  return result;
}

// A ray hitting the Bezier surface close to the rim may cross the underlying plane just beyond the rim edge, where no built neighbour
// would catch it. Such rays get a second chance without the planar limit, but only when the planar intersection is near beyond
// the rim edge, to let Newton start close to the hit.
BezierIntersection BezierMesh::intersectBeyondRim(Rim const &aRim, BezierTriangle const &aBezier, Ray const &aRay) const {
  BezierIntersection result;
  result.mWhat = BezierIntersection::What::cNone;
  auto const inPlane = aRim.mUnderlyingPlane.intersect(aRay);
  if(inPlane.mValid) {
    Vector const barycentric = aRim.mBarycentricInverse * inPlane.mPoint;
    if(barycentric(0) >= 0.0 && barycentric(1) >= 0.0 && barycentric(2) < 0.0 && barycentric(2) > -csRimBarycentricBand) {
      result = aBezier.intersectKajiya(aRay, BezierTriangle::LimitPlaneIntersection::cNone);
      if(result.mWhat == BezierIntersection::What::cIntersect && aRim.mRimPlane.distance(result.mIntersection.mPoint) < 0.0) {
        result.mWhat = BezierIntersection::What::cNone;            // Outside two divider planes, see BezierTriangle::classifyIntersection
      }
      else { // Nothing to do
      }
    }
    else { // Nothing to do
    }
  }
  else { // Nothing to do
  }
  return result;
}
