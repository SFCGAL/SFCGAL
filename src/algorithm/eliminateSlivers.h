// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_ELIMINATE_SLIVERS_H_
#define SFCGAL_ALGORITHM_ELIMINATE_SLIVERS_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL::algorithm {

/**
 * @brief Collapses the sliver faces of a volumic geometry.
 *
 * A boolean operation is computed exactly, but on meshes that are not: a
 * geometry built in floating point only touches its neighbour approximately,
 * so the intersection points land a rounding error away from the vertices they
 * should coincide with. The needle faces between them are valid on their own,
 * and fail later, when the shell is triangulated again or written with a
 * precision coarser than their own size.
 *
 * The boolean operations already clean up after themselves at the scale of
 * that rounding error. This entry point exists for the coarser case: a
 * geometry whose coordinates were rounded before the operation carries slivers
 * as large as that rounding, and only the caller knows how large it was.
 *
 * A geometry already written with a precision coarser than its slivers is
 * repaired too: both ends of a needle then share the same coordinates, which
 * folds its ring onto itself and leaves a patch without a plane. Those patches
 * are dropped first, they carried no area.
 *
 * Solids, polyhedral surfaces and triangulated surfaces are repaired, in place
 * for every shell; collections are walked through. Any other geometry is
 * returned unchanged. The result is triangulated.
 *
 * @param geometry      The geometry to repair.
 * @param maxEdgeLength Longest edge a collapse may consume, in the unit of the
 * coordinates. Zero or less asks for the default used by the boolean
 * operations, EPSILON_RELATIVE times the diagonal of the bounding box.
 *
 * @return A repaired copy of @p geometry.
 *
 * @warning A collapse moves a vertex, which can push it through a neighbouring
 * face. When that happens the shell is left as it was and a warning is
 * emitted: an intersecting shell is worse than the slivers it replaces.
 *
 * @since 3.0
 */
SFCGAL_API auto
eliminateSlivers(const Geometry &geometry, double maxEdgeLength = 0.0)
    -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_ELIMINATE_SLIVERS_H_
