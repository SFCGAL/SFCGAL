// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_TESSELLATE_H_
#define SFCGAL_ALGORITHM_TESSELLATE_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"

namespace SFCGAL::algorithm {
struct NoValidityCheck;

/**
 * Tessellate a geometry: this will triangulate surfaces (including polyhedral
 * and solid's surfaces) and keep untouched points, lines, etc.
 * @return tessellated geometry
 * @pre g is a valid geometry
 */
SFCGAL_API auto
tessellate(const Geometry &) -> std::unique_ptr<SFCGAL::Geometry>;

/**
 * Tessellate a geometry: this will triangulate surfaces (including polyhedral
 * and solid's surfaces) and keep untouched points, lines, etc.
 * @return tessellated geometry
 * @pre g is a valid geometry
 * @warning No actual validity check is done.
 */
SFCGAL_API auto
tessellate(const Geometry &, NoValidityCheck)
    -> std::unique_ptr<SFCGAL::Geometry>;

} // namespace SFCGAL::algorithm

#endif
