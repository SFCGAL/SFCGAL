// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/meshToPolyhedralSurface.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/detail/algorithm/meshToPolyhedralSurface.h"
#include "SFCGAL/version.h"

#if SFCGAL_CGAL_VERSION_NUM >= SFCGAL_CGAL_MAKE_VERSION(6, 2, 0)
  #include <CGAL/boost/graph/border.h>
#else
  #include <CGAL/Polygon_mesh_processing/border.h>
#endif
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

using FaceIndex     = Surface_mesh_3::Face_index;
using HalfedgeIndex = Surface_mesh_3::Halfedge_index;

/// @} end of private section

// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
meshToPolyhedralSurface(const Surface_mesh_3 &mesh, const Kernel::FT &epsAngle,
                        const Kernel::FT &epsDist)
    -> std::unique_ptr<PolyhedralSurface>
{
  auto result = std::make_unique<PolyhedralSurface>();

  // extract planar surfaces from the mesh
  std::vector<std::vector<FaceIndex>> coplanarGroups =
      detail::algorithm::groupCoplanarFaces(mesh, epsAngle, epsDist);

  // Extract the vertices of a face, constructs the corresponding ring,
  // and computes its signed area projected onto the face normal.
  // the are will be used to orient the ring consistently

  // for each surface, extract the boundary cycles
  // if there is only one cycle: no hole.
  // if there is more than one cycle, it means there are
  // holes
  for (const auto &planarIdx : coplanarGroups) {
    Polygon patch;

    // extract mesh cycles
    std::vector<HalfedgeIndex>                cyclesIdx;
    CGAL::Face_filtered_graph<Surface_mesh_3> filteredMesh(mesh, planarIdx);
#if SFCGAL_CGAL_VERSION_NUM >= SFCGAL_CGAL_MAKE_VERSION(6, 2, 0)
    CGAL::extract_boundary_cycles(filteredMesh, std::back_inserter(cyclesIdx));
#else
    PMP::extract_boundary_cycles(filteredMesh, std::back_inserter(cyclesIdx));
#endif

    // extract normal of current group
    Kernel::Vector_3 faceNormal =
        PMP::compute_face_normal(planarIdx.front(), mesh);

    // no holes
    // can directly set exteriorRing
    if (cyclesIdx.size() == 1) {
      LineString       exteriorRing;
      const Kernel::FT signedArea = detail::algorithm::createRing(
          mesh, filteredMesh, cyclesIdx[0], faceNormal, exteriorRing);
      if (signedArea < Kernel::FT(0)) {
        exteriorRing.reverse();
      }
      patch.setExteriorRing(exteriorRing);
    }

    // holes
    // extract all rings and sort by area
    // the largest one is the exterior ring
    else {
      std::vector<std::pair<LineString, Kernel::FT>> rings;
      for (const auto &cycleHalfedge : cyclesIdx) {
        LineString       ring;
        const Kernel::FT signedArea = detail::algorithm::createRing(
            mesh, filteredMesh, cycleHalfedge, faceNormal, ring);
        rings.emplace_back(ring, signedArea);
      }

      std::sort(rings.begin(), rings.end(),
                [](const auto &a, const auto &b) -> auto {
                  return CGAL::abs(a.second) > CGAL::abs(b.second);
                });

      for (std::size_t i = 0; i < rings.size(); ++i) {
        LineString       ring       = rings[i].first;
        const Kernel::FT signedArea = rings[i].second;
        if (i == 0) {
          if (signedArea < Kernel::FT(0)) {
            ring.reverse();
          }
          patch.setExteriorRing(ring);
        } else {
          if (signedArea > Kernel::FT(0)) {
            ring.reverse();
          }
          patch.addInteriorRing(ring);
        }
      }
    }

    result->addPatch(patch);
  }

  return result;
}
// NOLINTEND(readability-function-cognitive-complexity)

} // namespace SFCGAL::algorithm
