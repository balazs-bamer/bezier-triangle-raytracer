#ifndef BEZIER_TRIANGLE_RAYTRACER_MESH
#define BEZIER_TRIANGLE_RAYTRACER_MESH

#include "3dGeomUtil.h"
#include "stl_reader.h"
#include <map>
#include <list>
#include <array>
#include <deque>
#include <vector>
#include <string>
#include <fstream>
#include <limits>
#include <functional>
#include <unordered_map>
#include <unordered_set>


class Mesh final {
private:
  static constexpr double csStandardizeVerticesEpsilonFactor        = 0.2;
  static constexpr double csStandardizeNormalsEpsilon               = 0.01;
  static constexpr double csStandardizeNormalsIndependentMoveFactor = 0.2;
  static constexpr double csBoundaryVertexOnEdgeEpsilonFactor       = 0.01;

public:
  using TheMesh         = std::vector<Triangle>;
  using value_type      = Triangle;

  static constexpr uint32_t csNoNeighbour       = std::numeric_limits<uint32_t>::max(); // In Neighbours::mFellowTriangles for an edge on the rim of an open mesh.
  static constexpr uint8_t  csNoCommonSideStart = 3u;                                   // In Neighbours::mFellowCommonSideStarts for the same.

  struct Neighbours final {
    std::array<uint32_t, 3u> mFellowTriangles;        // Triangle indices: neighbour of edge (index, index + 1), or csNoNeighbour
    std::array<uint8_t, 3u>  mFellowCommonSideStarts; // Vertice index in each fellow triangle where the common side starts,
                                                      // such that the side in the neighbouring triangle is (i, i+1), or csNoCommonSideStart
  };

  struct VertexHash {
    std::size_t operator()(Vertex const &aVertex) const {
      return std::hash<double>{}(aVertex[0]) ^ (std::hash<double>{}(aVertex[1]) << 1u) ^ (std::hash<double>{}(aVertex[2]) << 2u);
    }
  };

  using Face2neighbours       = std::vector<Neighbours>;
  using Vertex2averageNormals = std::unordered_map<Vertex, Vector, VertexHash>;
  using Vertices              = std::unordered_set<Vertex, VertexHash>;

private:
  TheMesh               mMesh;
  Face2neighbours       mFace2neighbours;
  Vertex2averageNormals mVertex2averageNormals;
  Vertices              mBoundaryVertices;      // Vertices on the rim of an open mesh, empty for a closed one.

public:
  Mesh() = default;
  Mesh(Mesh &&) = default;
  Mesh(Mesh const&) = default;
  Mesh& operator=(Mesh &&) = default;
  Mesh& operator=(Mesh const&) = default;

  auto size() const { return mMesh.size(); }
  auto begin() const { return mMesh.begin(); }
  auto end() const { return mMesh.end(); }
  auto cbegin() const { return mMesh.cbegin(); }
  auto cend() const { return mMesh.cend(); }
  void reserve(uint32_t const aSize) { mMesh.reserve(aSize); }
  void push_back(Triangle const &aTriangle) { mMesh.push_back(aTriangle); }
  void clear() { mMesh.clear(); }
  auto& operator[](uint32_t const aI) { return mMesh[aI]; }
  auto const &operator[](uint32_t const aI) const { return mMesh[aI]; }

  TheMesh const&               getMesh() const                  { return mMesh; }
  Face2neighbours const&       getFace2neighbours() const       { return mFace2neighbours; }
  Vertex2averageNormals const& getVertex2averageNormals() const { return mVertex2averageNormals; } // TODO perhaps implement a way to import these values
                                                                                                   // along with the triangles when the original surface is known
                                                                                                   // and dervatives are present.
  void standardizeVertices();
  Vertices getVertices() const;

/// Makes all normalvectors (V[1]-V[0])x(V[2]-V[0]) point outwards.
/// Should run after standardizeVertices.
/// But not after auto divideLargeTriangles(double const aMaxTriangleSide), because it will probably put vertices on edges of other triangles.
/// Only this method populates the member variable mFace2neighbours and mVertex2averageNormals
/// For closed meshes only, throws if an edge belongs to only one triangle.
  void standardizeNormals();

/// The same for open meshes (and closed ones too), where "outwards" has no meaning. Instead the triangle facing
/// most in aFrontDirection (or opposite to it) will have its normal pointing to aFrontDirection, and all others
/// accordingly. Edges on the rim get csNoNeighbour, and their vertices form getBoundaryVertices().
  void standardizeNormals(Vector const &aFrontDirection);

/// These are valid after standardizeNormals.
  bool            isClosed() const            { return mBoundaryVertices.empty(); }
  Vertices const& getBoundaryVertices() const { return mBoundaryVertices; }
/// True if no vertex of the triangle is on the rim, so all its vertex normals are averages of complete fans.
  bool            isCoreFace(uint32_t const aIndex) const;

  void transform(Transform const &aTransform, Vertex const aDisplacement);

  Mesh& operator+=(Vector const aDisplacement)  { transform(Transform::Identity(),           aDisplacement);  return *this; }
  Mesh& operator*=(Transform const &aTransform) { transform(aTransform,                      Vertex::Zero()); return *this; }
  Mesh& operator*=(double const &aFactor)        { transform(Transform::Identity() * aFactor, Vertex::Zero()); return *this; }

/// WARNING! Likely to put new vertices on edges!
  void splitTriangles(double const aMaxTriangleSide);

/// Splits every triangle such that each will become (aDivisor + 1)(aDivisor + 2)/2 new small ones.
  void splitTriangles(int32_t const aDivisor);

  void readMesh(std::string const &aFilename);

  void writeMesh(std::string const &aFilename) const;

  /// aEnvelope is a C2 function in [-1, 1], with prefereably vertical tangents in -1 and 1 and values of 0 in -1 and 1.
  void makeSolidOfRevolution(int32_t const aSectors, int32_t const aBelts, std::function<double(double)> aEnvelope, Vector const &aSize);
  void makeEllipsoid(int32_t const aSectors, int32_t const aBelts, Vector const &aSize) { makeSolidOfRevolution(aSectors, aBelts, [](double const aX){ return ::sqrt(1 - aX * aX); }, aSize); }
  void makeUnitSphere(int32_t const aSectors, int32_t const aBelts) { makeEllipsoid(aSectors, aBelts, Vector(1.0, 1.0, 1.0)); }

private:
  double getSmallestSide() const;

  using ProjectedIndices = std::multimap<double, std::pair<uint32_t, uint32_t>>;
  using Iterator = typename ProjectedIndices::const_iterator;
  using Intervals = std::deque<std::pair<Iterator, Iterator>>;

  void projectVertices(ProjectedIndices &aCurrentProjectedSorted, int32_t const aDimension) const;
  double makeProximityIntervals(ProjectedIndices const &aCurrentProjectedSorted, Intervals &aCurrentIntervals, double const aEpsilon) const;
  void standardizeInIntervals(Intervals const &aIntervals, double const aEpsilonSquared);
  static uint32_t getIndependentFrom(Triangle const &aFaceTarget, Triangle const &aFaceOther);

  struct PairHash {
    template<typename t1, typename t2>
    std::size_t operator()(std::pair<t1, t2> const &aPair) const {
      return std::hash<t1>{}(aPair.first) ^ (std::hash<t2>{}(aPair.second) << 1u);
    }
  };

  using Edge2face = std::unordered_multimap<std::pair<uint32_t, uint32_t>, uint32_t, PairHash>;
  using Face2vertex = std::vector<std::array<uint32_t, 3u>>;
  using Vertex2index = std::unordered_map<Vertex, uint32_t, VertexHash>;

  std::tuple<Edge2face, Face2vertex, Vertex2index> createEdge2faceFace2vertex() const;
  std::pair<double, std::unordered_set<uint32_t>> getSmallestXstuff(Face2vertex const &aFace2vertex, Vertex2index const &) const;
  void createFace2neighbour(Edge2face const &aEdge2face, Face2vertex const &aFace2vertex, bool const aAllowOpen);
  void checkVerticesOnBoundaryEdges() const;
  void orientFaces(Edge2face const &aEdge2face, Face2vertex const &aFace2vertex, std::unordered_set<uint32_t> const &aSeedCandidates,
                   Vector const &aDesiredVector, bool const aAllowOpen);
  uint32_t getInitialFaceIndex(std::unordered_set<uint32_t> const &aFacesAtSmallestX, Vertex const &aDesiredVector) const;
  static void normalize(Triangle &aFace, Vertex const &aDesiredVector);
  static void normalize(Triangle const &aFaceKnown, Triangle &aFaceUnknown);
  void calculateNormalAverages4vertices();
  static void divideTriangle(TheMesh &aMesh, Triangle const &aTriangle, int32_t const aDivisor);
};

#endif
