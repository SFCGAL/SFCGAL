// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/algorithm/meshToPolyhedralSurface.h"

#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace SFCGAL::detail::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

template <typename MeshType>
using VertexIndex = typename boost::graph_traits<MeshType>::vertex_descriptor;

/**
 * Checks whether two faces are coplanar within given tolerances.
 *
 *
 * @param mesh       The input surface mesh.
 * @param face1      The first face to compare.
 * @param face2      The second face to compare.
 * @param normal1    normal vector of the first face.
 * @param normal2    normal vector of the second face.
 * @param epsAngle   Maximum allowed angle (in degrees) between the two normals.
 * @param epsDist    Maximum allowed distance between vertices of a face and the
 * plane of another face to consider them coplanar. of face1.
 *
 * @return true if the two faces are considered coplanar, false otherwise.
 */
template <typename MeshType>
static auto
areCoplanarFaces(const MeshType &mesh, FaceIndex<MeshType> face1,
                 FaceIndex<MeshType> face2, const Vector_3 &normal1,
                 const Vector_3 &normal2, const Kernel::FT &epsAngle,
                 const Kernel::FT &epsDist) -> bool
{
  // A degenerate face has a null normal, hence no plane: both tests below
  // would then compare 0 to 0 and swallow it into any group.
  if (normal1 == CGAL::NULL_VECTOR || normal2 == CGAL::NULL_VECTOR) {
    return false;
  }

  const Kernel::FT len1Sq = normal1.squared_length();
  const Kernel::FT len2Sq = normal2.squared_length();

  // Test if the normals are parallels:
  //   |n1.n2| >= cos(epsAngle) * |n1| * |n2|
  // squared on both sides to stay in exact arithmetic, since Kernel::FT has no
  // exact square root
  const Kernel::FT deg2rad = CGAL_PI / Kernel::FT(180.0);
  const Kernel::FT cosEps  = std::cos(CGAL::to_double(epsAngle * deg2rad));
  const Kernel::FT dot     = normal1 * normal2;
  if (dot * dot < cosEps * cosEps * len1Sq * len2Sq) {
    return false;
  }

  // Test if every vertex of face2 lies in the plane of face1:
  //   (n1.(pt - pt0))^2 <= epsDist^2 * |n1|^2
  auto           vertexPointMap = get(CGAL::vertex_point, mesh);
  const Point_3 &pt0 = get(vertexPointMap, source(halfedge(face1, mesh), mesh));
  const Kernel::FT maxDistSq = epsDist * epsDist * len1Sq;
  auto halfedges = CGAL::halfedges_around_face(halfedge(face2, mesh), mesh);
  return std::all_of(
      halfedges.begin(), halfedges.end(), [&](auto halfedge) -> auto {
        const Point_3   &pt      = get(vertexPointMap, source(halfedge, mesh));
        const Kernel::FT dotProj = normal1 * (pt - pt0);
        return dotProj * dotProj <= maxDistSq;
      });
}

/// @} end of private section

template <typename MeshType>
auto
groupCoplanarFaces(const MeshType &mesh, const Kernel::FT &epsAngle,
                   const Kernel::FT &epsDist)
    -> std::vector<std::vector<FaceIndex<MeshType>>>
{
  // A Polyhedron_3 face carries no index, so it cannot key a vector; the
  // dynamic property map falls back to a hash map only for those meshes.
  auto visited = get(CGAL::dynamic_face_property_t<bool>(), mesh);
  auto normals = get(CGAL::dynamic_face_property_t<Vector_3>(), mesh);

  for (FaceIndex<MeshType> face : faces(mesh)) {
    put(normals, face, PMP::compute_face_normal(face, mesh));
  }

  std::vector<std::vector<FaceIndex<MeshType>>> result;

  // Iterate over all faces to look for coplanar components
  for (const auto &seedFace : faces(mesh)) {
    if (get(visited, seedFace)) {
      continue;
    }

    std::queue<FaceIndex<MeshType>> queue;
    queue.push(seedFace);
    put(visited, seedFace, true);
    std::vector<FaceIndex<MeshType>> coplanarGroup;

    // look for connected coplanar faces
    while (!queue.empty()) {
      FaceIndex<MeshType> currentFace = queue.front();
      queue.pop();
      coplanarGroup.emplace_back(currentFace);

      // Iterate over all neighboring faces
      for (const auto &halfedge :
           CGAL::halfedges_around_face(halfedge(currentFace, mesh), mesh)) {

        FaceIndex<MeshType> oppositeFace = face(opposite(halfedge, mesh), mesh);

        if (oppositeFace == boost::graph_traits<MeshType>::null_face() ||
            get(visited, oppositeFace)) {
          continue;
        }

        if (areCoplanarFaces(mesh, seedFace, oppositeFace,
                             get(normals, seedFace), get(normals, oppositeFace),
                             epsAngle, epsDist)) {
          queue.push(oppositeFace);
          put(visited, oppositeFace, true);
        }
      }
    }

    result.push_back(std::move(coplanarGroup));
  }

  return result;
}

template <typename MeshType>
auto
createRing(const MeshType                            &mesh,
           const CGAL::Face_filtered_graph<MeshType> &filteredMesh,
           const HalfedgeIndex<MeshType>             &ringIdx,
           const Kernel::Vector_3 &faceNormal, LineString &ring) -> Kernel::FT
{
  for (const VertexIndex<MeshType> vertex :
       vertices_around_face(ringIdx, filteredMesh)) {
    ring.addPoint(Point(get(CGAL::vertex_point, mesh, vertex)));
  }

  Vector_3          sum        = CGAL::NULL_VECTOR;
  const std::size_t nrPts      = ring.numPoints();
  const Point_3     firstPoint = ring.startPoint().toPoint_3();
  for (std::size_t i = 0; i < nrPts; ++i) {
    sum += CGAL::cross_product(ring.pointN(i).toPoint_3() - firstPoint,
                               ring.pointN((i + 1) % nrPts).toPoint_3() -
                                   firstPoint);
  }
  ring.closes();

  return sum * faceNormal;
}

// Explicit instantiations
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template std::vector<
    std::vector<typename boost::graph_traits<Polyhedron_3>::face_descriptor>>
groupCoplanarFaces<Polyhedron_3>(const Polyhedron_3 &, const Kernel::FT &,
                                 const Kernel::FT &);

template std::vector<std::vector<
    typename boost::graph_traits<detail::MarkedPolyhedron>::face_descriptor>>
groupCoplanarFaces<detail::MarkedPolyhedron>(const detail::MarkedPolyhedron &,
                                             const Kernel::FT &,
                                             const Kernel::FT &);

template std::vector<
    std::vector<typename boost::graph_traits<Surface_mesh_3>::face_descriptor>>
groupCoplanarFaces<Surface_mesh_3>(const Surface_mesh_3 &, const Kernel::FT &,
                                   const Kernel::FT &);

template Kernel::FT
createRing<Polyhedron_3>(const Polyhedron_3 &,
                         const CGAL::Face_filtered_graph<Polyhedron_3> &,
                         const HalfedgeIndex<Polyhedron_3> &,
                         const Kernel::Vector_3 &, LineString &);

template Kernel::FT
createRing<detail::MarkedPolyhedron>(
    const detail::MarkedPolyhedron &,
    const CGAL::Face_filtered_graph<detail::MarkedPolyhedron> &,
    const HalfedgeIndex<detail::MarkedPolyhedron> &, const Kernel::Vector_3 &,
    LineString &);

template Kernel::FT
createRing<Surface_mesh_3>(const Surface_mesh_3 &,
                           const CGAL::Face_filtered_graph<Surface_mesh_3> &,
                           const HalfedgeIndex<Surface_mesh_3> &,
                           const Kernel::Vector_3 &, LineString &);
#endif

} // namespace SFCGAL::detail::algorithm
