// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_ALGORITHM_COLLAPSE_SLIVERS_H_
#define SFCGAL_DETAIL_ALGORITHM_COLLAPSE_SLIVERS_H_

#include "SFCGAL/detail/TypeForDimension.h"

namespace SFCGAL::detail::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

/**
 * Default bound of collapseSlivers(), relative to the diagonal of the bounding
 * box. It has to cover the rounding error of a mesh built in double precision,
 * about 1e-16 relative, and to stay far below anything a user models: on a
 * 100 m object it allows 10 nm.
 */
constexpr double SLIVER_RELATIVE_LENGTH = 1e-10;

/**
 * Collapses the sliver faces a boolean operation leaves behind.
 *
 * The corefinement is exact, but the meshes it is given rarely are: a cutter
 * built from floating point trigonometry only touches its target
 * approximately, so the intersection points land a rounding error away from
 * the vertices they should coincide with. The result is correct, and carries
 * needle faces whose shortest edge measures a fraction of a nanometre. They
 * survive until something downstream triangulates the shell again, or until
 * the geometry is written with a coarser precision than their own size, and
 * fail there instead of here.
 *
 * Only faces whose short edge is negligible against the extent of the mesh are
 * touched, so the collapse stays several orders of magnitude below any
 * modelled feature. A mesh free of such faces is left untouched.
 *
 * @param mesh The triangulated mesh to clean up, modified in place.
 * @param maxEdgeLength Longest edge a collapse may consume. A value below or
 * equal to zero asks for the default, SLIVER_RELATIVE_LENGTH times the
 * diagonal of the bounding box.
 */
void
collapseSlivers(MarkedPolyhedron &mesh, double maxEdgeLength = 0.0);

/// @} end of private section

} // namespace SFCGAL::detail::algorithm

#endif // SFCGAL_DETAIL_ALGORITHM_COLLAPSE_SLIVERS_H_
