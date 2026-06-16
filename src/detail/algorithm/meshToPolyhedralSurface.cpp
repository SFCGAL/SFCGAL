// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/algorithm/meshToPolyhedralSurface.h"

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
  // Test if the normals are parallels
  const Kernel::FT deg2rad = CGAL_PI / Kernel::FT(180.0);
  const Kernel::FT cosEps  = std::cos(CGAL::to_double(epsAngle * deg2rad));
  const Kernel::FT dot =
      std::clamp(normal1 * normal2, Kernel::FT(-1.0), Kernel::FT(1.0));
  if (CGAL::abs(dot) < cosEps) {
    return false;
  }

  // Test if a vertex of face2 lies in the plane of face1
  auto           vertexPointMap = get(CGAL::vertex_point, mesh);
  const Point_3 &pt0 = get(vertexPointMap, source(halfedge(face1, mesh), mesh));
  const CGAL::Plane_3<Kernel> plane(pt0, normal1);
  auto halfedges = CGAL::halfedges_around_face(halfedge(face2, mesh), mesh);
  return std::all_of(
      halfedges.begin(), halfedges.end(), [&](auto halfedge) -> auto {
        const Point_3   &pt = get(vertexPointMap, source(halfedge, mesh));
        const Kernel::FT signedDist = plane.a() * pt.x() + plane.b() * pt.y() +
                                      plane.c() * pt.z() + plane.d();
        return signedDist * signedDist <= epsDist * epsDist;
      });
}

/// @} end of private section

template <typename MeshType>
auto
groupCoplanarFaces(const MeshType &mesh, const Kernel::FT &epsAngle,
                   const Kernel::FT &epsDist)
    -> std::vector<std::vector<FaceIndex<MeshType>>>
{
  std::unordered_map<FaceIndex<MeshType>, bool>     visited;
  std::unordered_map<FaceIndex<MeshType>, Vector_3> normals;

  for (FaceIndex<MeshType> face : faces(mesh)) {
    visited[face] = false;
    normals[face] = PMP::compute_face_normal(face, mesh);
  }

  std::vector<std::vector<FaceIndex<MeshType>>> result;

  // Iterate over all faces to look for coplanar components
  for (const auto &seedFace : faces(mesh)) {
    if (visited[seedFace]) {
      continue;
    }

    std::queue<FaceIndex<MeshType>> queue;
    queue.push(seedFace);
    visited[seedFace] = true;
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
            visited[oppositeFace]) {
          continue;
        }

        if (areCoplanarFaces(mesh, seedFace, oppositeFace, normals[seedFace],
                             normals[oppositeFace], epsAngle, epsDist)) {
          queue.push(oppositeFace);
          visited[oppositeFace] = true;
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
    std::vector<typename boost::graph_traits<Surface_mesh_3>::face_descriptor>>
groupCoplanarFaces<Surface_mesh_3>(const Surface_mesh_3 &, const Kernel::FT &,
                                   const Kernel::FT &);

template Kernel::FT
createRing<Surface_mesh_3>(const Surface_mesh_3 &,
                           const CGAL::Face_filtered_graph<Surface_mesh_3> &,
                           const HalfedgeIndex<Surface_mesh_3> &,
                           const Kernel::Vector_3 &, LineString &);
#endif

} // namespace SFCGAL::detail::algorithm
