// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_MERGE_COPLANAR_FACES_H_
#define SFCGAL_ALGORITHM_MERGE_COPLANAR_FACES_H_

#include <memory>

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/export.h"

namespace SFCGAL::algorithm {

/**
 * Merges adjacent coplanar faces
 *
 * @param mesh The mesh to modify
 * @param epsAngle angular tolerance. Two faces whose normals
 *                 differ by less than this angle are considered parallel.
 * @param epsDist distance to plane tolerance
 *
 * @since 2.3
 */
SFCGAL_API
void
mergeCoplanarFaces(Surface_mesh_3   &mesh,
                   const Kernel::FT &epsAngle = Kernel::FT(0.5),
                   const Kernel::FT &epsDist  = Kernel::FT(1e-8));

/**
 * Merges adjacent coplanar faces
 *
 * @param geometry The geometry to modify
 * @param epsAngle angular tolerance in degrees. Two faces whose normals
 *                 differ by less than this angle are considered parallel.
 * @param epsDist distance to plane tolerance
 *
 * @return A unique pointer to the simplified geometry
 *
 * @since 2.3
 */
SFCGAL_API
auto
mergeCoplanarFaces(const Geometry   &geometry,
                   const Kernel::FT &epsAngle = Kernel::FT(0.5),
                   const Kernel::FT &epsDist  = Kernel::FT(1e-8))
    -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_MERGE_COPLANAR_FACES_H_
