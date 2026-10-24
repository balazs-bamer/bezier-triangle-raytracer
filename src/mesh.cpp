#include "mesh.h"


double Mesh::getSmallestSide() const {
  double smallestSide = std::numeric_limits<double>::max();
  for(auto &triangle : mMesh) {
    for(uint32_t i = 0u; i < 3u; ++i) {
      smallestSide = std::min(smallestSide, (triangle[i] - triangle[(i + 1u) % 3u]).norm());
    }
  }
  return smallestSide;
}

void Mesh::projectVertices(ProjectedIndices &aCurrentProjectedSorted, int32_t const aDimension) const {
  for(uint32_t indexFace = 0u; indexFace < mMesh.size(); ++indexFace) {
    auto &face = mMesh[indexFace];
    for(uint32_t indexVertex = 0u; indexVertex < 3u; ++indexVertex) {
      auto &vertex = face[indexVertex];
      aCurrentProjectedSorted.emplace(std::make_pair(vertex[aDimension], std::make_pair(indexFace, indexVertex)));
    }
  }
}

double Mesh::makeProximityIntervals(ProjectedIndices const &aCurrentProjectedSorted, Intervals &aCurrentIntervals, double const aEpsilon) const {
  bool was = false;
  double startValue;
  uint32_t max = 0u, counter;
  Iterator startIterator = aCurrentProjectedSorted.cend();;
  for(Iterator i = aCurrentProjectedSorted.cbegin(); i != aCurrentProjectedSorted.cend(); ++i) {
    if(was) {
      auto valueNow = i->first;
      if(valueNow - startValue >= aEpsilon) {
        aCurrentIntervals.emplace_back(std::make_pair(startIterator, i));
        startValue = valueNow;
        startIterator = i;
        max = std::max(max, counter);
        counter = 1u;
      }
      else {
        ++counter;
      }
    }
    else {
      startIterator = i;
      startValue = i->first;
      was = true;
      counter = 1u;
    }
  }
  Iterator last = aCurrentProjectedSorted.cend();
  aCurrentIntervals.emplace_back(std::make_pair(startIterator, last));
  max = std::max(max, counter);
  return max;
}

void Mesh::standardizeInIntervals(Intervals const &aIntervals, double const aEpsilonSquared) {
  for(auto const& interval : aIntervals) {                               // We need to check pairwise proximity only inside each interval.
    for(auto i = interval.first; i != interval.second; ++i) {
      auto &v1 = mMesh[i->second.first][i->second.second];
      for(auto j = interval.first; j != interval.second; ++j) {
        auto &v2 = mMesh[j->second.first][j->second.second];
        if((v1-v2).squaredNorm() < aEpsilonSquared && (v1[0] < v2[0] || (v1[0] == v2[0] && v1[1] < v2[1]) || (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] < v2[2]))) {
          v1 = v2;                                                        // Use lexicographic sorting to choose an unique one among the points too close to each other.
        }
        else { // nothing to do
        }
      }
    }
  }
}

void Mesh::standardizeVertices() {
  auto epsilon = getSmallestSide() * csStandardizeVerticesEpsilonFactor; // Vertices closer to each other than this will be forced to one point.

  std::array<uint32_t, 3> maximums;
  std::array<ProjectedIndices, 3u> projectedSorted;                           // Holds the projected indices sorted by a coordinate, separate for x, y and z.
  std::array<Intervals, 3u> proximityIntervals;                               // Holds iterator pairs to the previous stuff in between the projections are closer to each other than epsilon.

  for(uint32_t dimension = 0u; dimension < 3u; ++dimension) {
    auto& currentProjectedSorted = projectedSorted[dimension];
    auto& currentIntervals = proximityIntervals[dimension];

    projectVertices(currentProjectedSorted, dimension);
    maximums[dimension] = makeProximityIntervals(currentProjectedSorted, currentIntervals, epsilon);
  }
  auto minIt = std::min_element(maximums.cbegin(), maximums.cend());          // Best coordinate is where the maximum interval length is minimal.
  auto whichDim = minIt - maximums.cbegin();
  auto const& bestIntervals = proximityIntervals[whichDim];

  standardizeInIntervals(bestIntervals, epsilon * epsilon);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typename Mesh::Vertices Mesh::getVertices() const {
  Vertices result;
  for(auto const &triangle : mMesh) {
    for(auto const &vertex : triangle) {
      result.insert(vertex);
    }
  }
  return result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

uint32_t Mesh::getIndependentFrom(Triangle const &aFaceTarget, Triangle const &aFaceOther) {
  for(uint32_t i = 0u; i < 3u; ++i) {
    if(std::find(aFaceOther.cbegin(), aFaceOther.cend(), aFaceTarget[i]) == aFaceOther.cend()) {
      return i;
    }
    else { // Nothing to do
    }
  }
  return 3u;                                                          // Won't get here, but avoid compiler complaint.
}

std::tuple<typename Mesh::Edge2face, typename Mesh::Face2vertex, typename Mesh::Vertex2index> Mesh::createEdge2faceFace2vertex() const {
  Edge2face edge2face;
  Face2vertex face2vertex;
  face2vertex.reserve(mMesh.size());
  std::unordered_map<Vertex, uint32_t, VertexHash> was;               // For presence testing, maps each vertex to its index in the deque below.
  for(uint32_t indexFace = 0u; indexFace < mMesh.size(); ++indexFace) {
    auto const &face = mMesh[indexFace];
    std::array<uint32_t, 3u> faceVertexIndices;                       // We collect here the vertex indices for the actual face.
    for(uint32_t indexInFace = 0u; indexInFace < 3u; ++indexInFace) {
      auto const vertex = face[indexInFace];
      uint32_t indexInVertices;
      auto found = was.find(vertex);
      if(found == was.end()) {
        indexInVertices = was.size();
        faceVertexIndices[indexInFace] = indexInVertices;
        was.emplace(std::make_pair(vertex, indexInVertices));
      }
      else {
        indexInVertices = found->second;
        faceVertexIndices[indexInFace] = indexInVertices;
      }
    }
    face2vertex.push_back(faceVertexIndices);
    for(uint32_t i = 0u; i < 3u; ++i) {
      auto low = faceVertexIndices[i];
      auto high = faceVertexIndices[(i + 1u) % 3u];
      if(low > high) {
        std::swap(low, high);
      }
      else { // Nothing to do
      }
      edge2face.emplace(std::make_pair(std::pair(low, high), indexFace));
    }
  }
  return std::make_tuple(edge2face, face2vertex, was);
}

std::pair<double, std::unordered_set<uint32_t>> Mesh::getSmallestXstuff(Face2vertex const &aFace2vertex, Vertex2index const &aVertex2index) const {
  double smallestX = std::numeric_limits<double>::max();
  uint32_t smallestXverticeIndex;
  for(uint32_t indexFace = 0u; indexFace < mMesh.size(); ++indexFace) {
    auto const &face = mMesh[indexFace];
    for(uint32_t indexInFace = 0u; indexInFace < 3u; ++indexInFace) {
      auto const vertex = face[indexInFace];
      if(vertex[0] < smallestX) {
        smallestX = vertex[0];
        smallestXverticeIndex = aVertex2index.at(vertex);
      }
      else { // Nothing to do
      }
    }
  }
  std::unordered_set<uint32_t> facesAtSmallestX;
  for(uint32_t indexFace = 0u; indexFace < mMesh.size(); ++indexFace) {
    auto const &face = aFace2vertex[indexFace];
    Neighbours neighbours;
    for(uint32_t indexInFace = 0u; indexInFace < 3u; ++indexInFace) {
      if(face[indexInFace] == smallestXverticeIndex) {
        facesAtSmallestX.insert(indexFace);
      }
      else { // Nothing to do
      }
    }
  }
  return std::make_pair(smallestX, facesAtSmallestX);
}

void Mesh::createFace2neighbour(Edge2face const &aEdge2face, Face2vertex const &aFace2vertex, bool const aAllowOpen) {
  mFace2neighbours.clear();
  mFace2neighbours.reserve(mMesh.size());
  mBoundaryVertices.clear();
  for(uint32_t indexFace = 0u; indexFace < mMesh.size(); ++indexFace) {
    auto const &face = aFace2vertex[indexFace];
    Neighbours neighbours;
    for(uint32_t indexInFace = 0u; indexInFace < 3u; ++indexInFace) {
      auto vertexIndex0 = face[indexInFace];
      auto vertexIndex1 = face[(indexInFace + 1u) % 3u];
      auto edge = std::make_pair(vertexIndex0, vertexIndex1);
      if(edge.first > edge.second) {
        std::swap(edge.first, edge.second);
      }
      else { // Nothing to do
      }
      auto[begin, end] = aEdge2face.equal_range(edge);
      auto const faceCount = std::distance(begin, end);
      if(faceCount > 2) {
        throw "Edge belonging to more than two triangles detected.";
      }
      else if(faceCount == 1) {                                       // Only this face, so the edge is on the rim.
        if(aAllowOpen) {
          neighbours.mFellowTriangles[indexInFace] = csNoNeighbour;
          neighbours.mFellowCommonSideStarts[indexInFace] = csNoCommonSideStart;
          mBoundaryVertices.insert(mMesh[indexFace][indexInFace]);
          mBoundaryVertices.insert(mMesh[indexFace][(indexInFace + 1u) % 3u]);
        }
        else {
          throw "Vertex on edge detected.";
        }
      }
      else {
        if(begin->second == indexFace) {                              // We need the other, neighbouring face.
          ++begin;
        }
        else { // Nothing to do
        }
        auto otherFaceIndex = begin->second;
        neighbours.mFellowTriangles[indexInFace] = otherFaceIndex;    // Neighbour of edge (index, index + 1)
        auto const &otherFace = aFace2vertex[otherFaceIndex];
        auto otherFaceRawIndex0 = std::find(otherFace.cbegin(), otherFace.cend(), vertexIndex0) - otherFace.cbegin();
        auto otherFaceRawIndex1 = std::find(otherFace.cbegin(), otherFace.cend(), vertexIndex1) - otherFace.cbegin();
        uint32_t const resolve[3u][3u] = {{3u, 0u, 2u}, {0u, 3u, 1u}, {2u, 1u, 3u}};

        neighbours.mFellowCommonSideStarts[indexInFace] = resolve[otherFaceRawIndex0][otherFaceRawIndex1];
      }
    }
    mFace2neighbours.push_back(neighbours);
  }
}

// In an open mesh a vertex lying on an other triangle's edge would look like a rim, so we need to check it explicitly.
void Mesh::checkVerticesOnBoundaryEdges() const {
  std::vector<std::pair<Vertex, Vertex>> boundaryEdges;
  for(uint32_t indexFace = 0u; indexFace < mMesh.size(); ++indexFace) {
    for(uint32_t indexInFace = 0u; indexInFace < 3u; ++indexInFace) {
      if(mFace2neighbours[indexFace].mFellowTriangles[indexInFace] == csNoNeighbour) {
        boundaryEdges.emplace_back(std::make_pair(mMesh[indexFace][indexInFace], mMesh[indexFace][(indexInFace + 1u) % 3u]));
      }
      else { // Nothing to do
      }
    }
  }
  if(!boundaryEdges.empty()) {
    double const epsilon = getSmallestSide() * csBoundaryVertexOnEdgeEpsilonFactor;
    for(auto const &vertex : mBoundaryVertices) {
      for(auto const &[start, end] : boundaryEdges) {
        if(vertex != start && vertex != end) {
          Vector const side = end - start;
          double const proportion = (vertex - start).dot(side) / side.squaredNorm();
          if(proportion > 0.0 && proportion < 1.0 && (start + side * proportion - vertex).norm() < epsilon) {
            throw "Vertex on edge detected.";
          }
          else { // Nothing to do
          }
        }
        else { // Nothing to do
        }
      }
    }
  }
  else { // Nothing to do
  }
}

bool Mesh::isCoreFace(uint32_t const aIndex) const {
  auto const &face = mMesh[aIndex];
  return std::none_of(face.cbegin(), face.cend(), [this](Vertex const &aVertex){ return mBoundaryVertices.count(aVertex) > 0u; });
}

uint32_t Mesh::getInitialFaceIndex(std::unordered_set<uint32_t> const &aFacesAtSmallestX, Vertex const &aDesiredVector) const {
  double maxAbsoluteDotProduct = -std::numeric_limits<double>::max();
  uint32_t initialFaceIndex;                                          // First determine the initial face for which we want (f[1]-f[0])x(f[2]-f[0]) point outwards.
  for(auto const indexFace : aFacesAtSmallestX) {
    auto const &face = mMesh[indexFace];
    auto normal = util::getNormal(face).normalized();
    auto absoluteDotProduct = std::abs(aDesiredVector.dot(normal));
    if(absoluteDotProduct > maxAbsoluteDotProduct) {
      maxAbsoluteDotProduct = absoluteDotProduct;
      initialFaceIndex = indexFace;
    }
    else { // Nothing to do
    }
  }
  return initialFaceIndex;
}

void Mesh::normalize(Triangle &aFace, Vertex const &aDesiredVector) {
  auto normal = util::getNormal(aFace);
  if(aDesiredVector.dot(normal) < 0.0) {
    std::swap(aFace[0u], aFace[1u]);
  }
  else { // Nothing to do
  }
}

void Mesh::normalize(Triangle const &aFaceKnown, Triangle &aFaceUnknown) {  // Calculates everything twice for each triangle, but won't cache now.
  uint32_t independentIndexKnown = getIndependentFrom(aFaceKnown, aFaceUnknown);
  uint32_t commonIndex1known = (independentIndexKnown + 1u) % 3u;
  uint32_t commonIndex2known = (independentIndexKnown + 2u) % 3u;
  uint32_t independentIndexUnknown = getIndependentFrom(aFaceUnknown, aFaceKnown);
  uint32_t commonIndex1unknown = (independentIndexUnknown + 1u) % 3u;
  uint32_t commonIndex2unknown = (independentIndexUnknown + 2u) % 3u;

  auto altitudeKnown = util::getAltitude(aFaceKnown[commonIndex1known], aFaceKnown[commonIndex2known], aFaceKnown[independentIndexKnown]);
  Vertex altitudeUnknown = util::getAltitude(aFaceUnknown[commonIndex1unknown], aFaceUnknown[commonIndex2unknown], aFaceUnknown[independentIndexUnknown]);
  auto dotFaceAltitudes = altitudeKnown.dot(altitudeUnknown);

  auto normalKnown = util::getNormal(aFaceKnown);
  auto normalUnknown = util::getNormal(aFaceUnknown);
  auto knownDotUnknown = normalKnown.dot(normalUnknown);
  if(std::abs(knownDotUnknown / (normalKnown.norm() * normalUnknown.norm())) < csStandardizeNormalsEpsilon) {
    auto newIndependentUnknown = aFaceUnknown[independentIndexUnknown] + csStandardizeNormalsIndependentMoveFactor * (aFaceKnown[independentIndexKnown] - (aFaceKnown[commonIndex1known] + aFaceKnown[commonIndex2known]) / 2.0);
    auto newFaceUnknown = aFaceUnknown;
    newFaceUnknown[independentIndexUnknown] = newIndependentUnknown;

    altitudeUnknown = util::getAltitude(aFaceUnknown[commonIndex1unknown], aFaceUnknown[commonIndex2unknown], newIndependentUnknown);
    dotFaceAltitudes = altitudeKnown.dot(altitudeUnknown);
    normalUnknown = util::getNormal(newFaceUnknown);
    knownDotUnknown = normalKnown.dot(normalUnknown);
  }
  else { // Nothing to do
  }
  if(dotFaceAltitudes * knownDotUnknown > 0.0) {
    std::swap(aFaceUnknown[commonIndex1unknown], aFaceUnknown[commonIndex2unknown]);
  }
  else { // Nothing to do
  }
}

void Mesh::calculateNormalAverages4vertices() {
  std::unordered_multimap<Vertex, uint32_t, VertexHash> vertex2triangleIndex;
  for(uint32_t indexInMesh = 0u; indexInMesh < mMesh.size(); ++indexInMesh) {
    auto const &triangle = mMesh[indexInMesh];
    for(uint32_t i = 0u; i < 3u; ++i) {
      vertex2triangleIndex.emplace(std::make_pair(triangle[i], indexInMesh));
    }
  }
  for (auto iter=vertex2triangleIndex.cbegin(); iter != vertex2triangleIndex.cend(); ) {
    auto const &vertex = iter->first;
    auto const range = vertex2triangleIndex.equal_range(vertex);
    Vector sum = Vector::Zero();
    for (auto inRange = range.first; inRange != range.second; ++inRange) {
      auto const &triangle = mMesh[inRange->second];
      uint32_t const whichVertex = std::find(triangle.cbegin(), triangle.cend(), vertex) - triangle.cbegin();
      Vector sideA = triangle[(whichVertex + 1u) % 3u] - triangle[whichVertex];         // We use angle here
      Vector sideB = triangle[(whichVertex + 2u) % 3u] - triangle[whichVertex];         // Because small triangles are equally important
      double cosAngle = sideA.dot(sideB) / (sideA.norm() * sideB.norm());                // as large ones since they probably help the surface bend.
      sum += util::getNormal(mMesh[inRange->second]).normalized() * acos(cosAngle);
    }
    sum.normalize();
    mVertex2averageNormals.emplace(std::make_pair(vertex, sum));
    iter = range.second;
  }
}

void Mesh::standardizeNormals() {
  auto [edge2face, face2vertex, vertex2index] = createEdge2faceFace2vertex();
  auto [smallestX, facesAtSmallestX] = getSmallestXstuff(face2vertex, vertex2index);

  Vertex desiredVector;                                                 // A known unit vector pointing outwards for a definite face.
  desiredVector << -1.0, 0.0, 0.0;
  orientFaces(edge2face, face2vertex, facesAtSmallestX, desiredVector, false);
}

void Mesh::standardizeNormals(Vector const &aFrontDirection) {
  auto [edge2face, face2vertex, vertex2index] = createEdge2faceFace2vertex();
  std::unordered_set<uint32_t> allFaces;                                // Any face can be the initial one.
  for(uint32_t indexFace = 0u; indexFace < mMesh.size(); ++indexFace) {
    allFaces.insert(indexFace);
  }
  orientFaces(edge2face, face2vertex, allFaces, aFrontDirection, true);
}

void Mesh::orientFaces(Edge2face const &aEdge2face, Face2vertex const &aFace2vertex, std::unordered_set<uint32_t> const &aSeedCandidates,
                       Vector const &aDesiredVector, bool const aAllowOpen) {
  {
    createFace2neighbour(aEdge2face, aFace2vertex, aAllowOpen);
    checkVerticesOnBoundaryEdges();

    uint32_t initialFaceIndex = getInitialFaceIndex(aSeedCandidates, aDesiredVector);

    struct KnownUnknownFacePair {
      uint32_t mKnownIndex;
      uint32_t mUnknownIndex;
    };
    std::list<KnownUnknownFacePair> queue;                              // Faces with normals yet to be normalized.
    normalize(mMesh[initialFaceIndex], aDesiredVector);
    std::vector<bool> remaining;
    remaining.reserve(mFace2neighbours.size());
    std::fill_n(std::back_inserter(remaining), mFace2neighbours.size(), true);
    for(auto const indexFace : mFace2neighbours[initialFaceIndex].mFellowTriangles) {
      if(indexFace != csNoNeighbour) {
        queue.emplace_back(KnownUnknownFacePair{initialFaceIndex, indexFace});
      }
      else { // Nothing to do
      }
    }
    remaining[initialFaceIndex] = false;
    while(!queue.empty()) {
      auto actualPair = queue.back();
      queue.pop_back();
      if(remaining[actualPair.mUnknownIndex]) {
        normalize(mMesh[actualPair.mKnownIndex], mMesh[actualPair.mUnknownIndex]);
      }
      else { // Nothing to do
      }
      remaining[actualPair.mUnknownIndex] = false;
      for(auto const indexFace : mFace2neighbours[actualPair.mUnknownIndex].mFellowTriangles) {
        if(indexFace != csNoNeighbour && remaining[indexFace] && actualPair.mUnknownIndex != indexFace) {
          queue.emplace_back(KnownUnknownFacePair{actualPair.mUnknownIndex, indexFace});
        }
        else { // Nothing to do
        }
      }
    }
    if(std::find(remaining.cbegin(), remaining.cend(), true) != remaining.cend()) {
      throw "Mesh consisting of more than one piece detected.";
    }
    else { // Nothing to do
    }
  }
  {
    auto [edge2face, face2vertex, vertex2index] = createEdge2faceFace2vertex();
    createFace2neighbour(edge2face, face2vertex, aAllowOpen);           // Need to call it once more, because we swapped vertices of triangles
  }
  calculateNormalAverages4vertices();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Mesh::transform(Transform const &aTransform, Vertex const aDisplacement) {
  for(auto &triangle : mMesh) {
    for(auto &vertex : triangle) {
      vertex = aTransform * vertex + aDisplacement;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Mesh::divideTriangle(TheMesh &result, Triangle const &aTriangle, int32_t const aDivisor) {
  util::divide(aTriangle, aDivisor, [&result](Triangle && aNew){ result.push_back(aNew); } );
}

void Mesh::splitTriangles(double const aMaxTriangleSide) {
  TheMesh result;
  for(auto const &triangle : mMesh) {
    double maxSide = (triangle[0] - triangle[1]).norm();
    maxSide = std::max(maxSide, (triangle[0] - triangle[2]).norm());
    maxSide = std::max(maxSide, (triangle[1] - triangle[2]).norm());
    int32_t divisor = static_cast<int32_t>(std::ceil(maxSide / aMaxTriangleSide));
    divideTriangle(result, triangle, divisor);
  }
  mMesh = result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Mesh::splitTriangles(int32_t const aDivisor) {
  TheMesh result;
  for(auto const &triangle : mMesh) {
    divideTriangle(result, triangle, aDivisor);
  }
  mMesh = result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Mesh::readMesh(std::string const &aFilename) {
  mMesh.clear();
  mFace2neighbours.clear();
  stl_reader::StlMesh<double, int32_t> mesh(aFilename);
  for(int32_t indexTriangle = 0; indexTriangle < mesh.num_tris(); ++indexTriangle) {
    Triangle triangle;
    for(int32_t indexCorner = 0; indexCorner < 3; ++indexCorner) {
      double const * const coords = mesh.tri_corner_coords(indexTriangle, indexCorner);
      Vertex in;
      for(int32_t i = 0; i < 3; ++i) {
        in(i) = coords[i];
      }
      triangle[indexCorner] = in;
    }
    mMesh.push_back(triangle);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Mesh::writeMesh(std::string const& aFilename) const {
  std::ofstream out(aFilename);
  out << "solid Exported from Blender-2.82 (sub 7)\n";
  for(auto const & triangle : mMesh) {
    out << "facet normal 0.000000 0.000000 0.000000\nouter loop\n";
    for(auto const & vertex : triangle) {
      out << "vertex " << vertex(0) << ' ' << vertex(1) << ' ' << vertex(2) << '\n';
    }
    out << "endloop\nendfacet\n";
  }
  out << "endsolid Exported from Blender-2.82 (sub 7)\n";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Mesh::makeSolidOfRevolution(int32_t const aSectors, int32_t const aBelts, std::function<double(double)> aEnvelope, Vector const &aSize) {
  mMesh.clear();
  mMesh.reserve(static_cast<uint32_t>(2 * aSectors * aBelts));
  mFace2neighbours.clear();
  double sectorAngleHalf = cgPi / aSectors;
  double sectorAngleFull = sectorAngleHalf * 2.0;
  double beltAngle       = cgPi / (aBelts + 1.0);
  double bias = 0.0;
  double beltAngleUp = 0.0;
  double beltAngleMiddle = beltAngle;
  double beltAngleDown = 2.0 * beltAngle;
  double beltRadiusUp = 0.0;
  double beltRadiusMiddle = aSize(0) * aEnvelope(std::cos(beltAngleMiddle));
  double beltRadiusDown = aSize(0) * aEnvelope(std::cos(beltAngleDown));
  double beltZup = aSize(2);
  double beltZmiddle = aSize(2) * std::cos(beltAngleMiddle);
  double beltZdown = aSize(2) * std::cos(beltAngleDown);
  for(int32_t belt = 0; belt < aBelts; ++belt) {
    double sectorAngleUpDown = bias + sectorAngleHalf;
    double sectorAngleMiddle1 = bias + 0.0;
    double sectorAngleMiddle2 = bias + sectorAngleFull;
    for(int32_t sector = 0; sector < aSectors; ++sector) {
      Vertex corner1(beltRadiusUp * std::sin(sectorAngleUpDown), aSize(1) * beltRadiusUp * std::cos(sectorAngleUpDown), beltZup);
      Vertex corner2(beltRadiusMiddle * std::sin(sectorAngleMiddle1), aSize(1) * beltRadiusMiddle * std::cos(sectorAngleMiddle1), beltZmiddle);
      Vertex corner3(beltRadiusMiddle * std::sin(sectorAngleMiddle2), aSize(1) * beltRadiusMiddle * std::cos(sectorAngleMiddle2), beltZmiddle);
      mMesh.push_back({corner1, corner2, corner3});
      corner1 = {aSize(0) * beltRadiusDown * std::sin(sectorAngleUpDown), aSize(1) * beltRadiusDown * std::cos(sectorAngleUpDown), beltZdown};
      mMesh.push_back({corner2, corner3, corner1});
      sectorAngleUpDown += sectorAngleFull;
      sectorAngleMiddle1 = sectorAngleMiddle2;
      sectorAngleMiddle2 += sectorAngleFull;
    }
    beltAngleUp = beltAngleMiddle;
    beltAngleMiddle = beltAngleDown;
    beltAngleDown += beltAngle;
    beltRadiusUp = beltRadiusMiddle;
    beltRadiusMiddle = beltRadiusDown;
    beltRadiusDown = aSize(0) * aEnvelope(std::cos(beltAngleDown));
    beltZup = beltZmiddle;
    beltZmiddle = beltZdown;
    beltZdown = aSize(2) * std::cos(beltAngleDown);
    bias += sectorAngleHalf;
  }
}
