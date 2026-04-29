// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_DETAIL_MERGE_COPLANAR_FACES_H_
#define SFCGAL_ALGORITHM_DETAIL_MERGE_COPLANAR_FACES_H_

#include "SFCGAL/Kernel.h"

namespace SFCGAL::algorithm::detail {

using faceDescriptor     = Surface_mesh_3::Face_index;
using halfedgeDescriptor = Surface_mesh_3::Halfedge_index;

/**
 * Checks whether two faces are coplanar within given tolerances.
 *
 * @param mesh: input surface mesh
 * @param face1 first face to compare
 * @param face2 second face to compare
 * @param epsAngle angular tolerance. Two faces whose normals
 *                 differ by less than this angle are considered parallel.
 * @param epsDist geometric tolerance for distance to plane
 * @return true if the two faces can be considered coplanar, false otherwise.
 */
auto
areCoplanarFaces(const Surface_mesh_3 &mesh, faceDescriptor face1,
                 faceDescriptor face2, const Kernel::FT &epsAngle,
                 const Kernel::FT &epsDist) -> bool;

/**
 * Checks the preconditions to join faces on edge.
 *
 * join face requires:
 *   - the edge is not a border edge
 *   - both halfedges are incident to valid and distinct faces
 *   - the two vertices of the edge have degree >= 3
 *     (otherwise the mesh topology would become invalid)
 *
 * @param mesh input surface mesh
 * @param halfedge edge to inspect
 * @return true if the two faces can be considered coplanar, false otherwise.
 */
auto
canJoinFace(const Surface_mesh_3 &mesh, halfedgeDescriptor halfedge) -> bool;

} // namespace SFCGAL::algorithm::detail

#endif // SFCGAL_ALGORITHM_DETAIL_MERGE_COPLANAR_FACES_H_
