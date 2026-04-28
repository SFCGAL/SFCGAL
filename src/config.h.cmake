/**
 *   SFCGAL
 *
 *   Copyright (C) 2012-2013 Oslandia <infos@oslandia.com>
 *   Copyright (C) 2012-2013 IGN (http://www.ign.fr)
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Library General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Library General Public License for more details.

 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */
#ifndef SFCGAL_CONFIG_H_
#define SFCGAL_CONFIG_H_

#ifndef CGAL_USE_GMPXX
#define CMAKE_OVERRIDDEN_DEFAULT_ENT_BACKEND 0 // GMP
#endif

#include "SFCGAL/export.h"

/**
 * Maximum recursion depth for processing nested geometry collections
 */
#define SFCGAL_MAX_RECURSION_DEPTH @SFCGAL_MAX_RECURSION_DEPTH@
/**
 * Maximum number of GMP limbs allowed during deserialization
 */
#define SFCGAL_MAX_GMP_LIMBS @SFCGAL_MAX_GMP_LIMBS@
/**
 * Maximum total number of coordinates in WKB to prevent memory exhaustion
 */
#define SFCGAL_MAX_WKB_TOTAL_COORDINATES @SFCGAL_MAX_WKB_TOTAL_COORDINATES@
/**
 * Maximum total number of geometry elements in WKB to prevent memory exhaustion
 */
#define SFCGAL_MAX_WKB_TOTAL_ELEMENTS @SFCGAL_MAX_WKB_TOTAL_ELEMENTS@
/**
 * Maximum number of vertices allowed in OBJ file
 */
#define SFCGAL_MAX_OBJ_VERTICES @SFCGAL_MAX_OBJ_VERTICES@
/**
 * Maximum number of faces allowed in OBJ file
 */
#define SFCGAL_MAX_OBJ_FACES @SFCGAL_MAX_OBJ_FACES@
/**
 * Maximum number of lines allowed in OBJ file
 */
#define SFCGAL_MAX_OBJ_LINES @SFCGAL_MAX_OBJ_LINES@
/**
 * Maximum number of points allowed in OBJ file
 */
#define SFCGAL_MAX_OBJ_POINTS @SFCGAL_MAX_OBJ_POINTS@


/**
 * indicates if OpenSceneGraph dependency is activated
 */
#cmakedefine SFCGAL_WITH_OSG

#endif
