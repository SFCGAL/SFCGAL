// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_ROOFGENERATION_H_
#define SFCGAL_ALGORITHM_ROOFGENERATION_H_

#include <memory>
#include <vector>

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/config.h"

namespace SFCGAL {
class Point;
}

namespace SFCGAL::algorithm {

struct NoValidityCheck;

/**
 * @brief Enumeration of roof types supported by the roof generation algorithms.
 */
enum class RoofType : std::uint8_t {
  FLAT,     ///< Flat roof (simple extrusion)
  HIPPED,   ///< Hipped roof (straight skeleton extrusion)
  SKILLION, ///< Skillion roof (alias for pitched roof)
  GABLE     ///< Gable roof (dual symmetric slopes)
};

/**
 * @brief Parameters for roof generation algorithms.
 */
struct RoofParameters {
  RoofType type       = RoofType::GABLE; ///< Roof Type (cf. RoofType)
  double   slopeAngle = 30.0;            ///< Slope angle in degrees (0-90) for
                                         ///< gable/skillion roofs
  double roofHeight = 3.0; ///< Maximum roof height for flat/hipped roofs (>= 0)
  size_t primaryEdgeIndex = 0; ///< Index of the primary slope edge (skillion)
};

/**
 * @brief Generate a roof using unified parameters.
 * @param footprint The 2D polygon footprint of the building
 * @param params The parameters for roof generation
 * @return A unique_ptr to the generated roof geometry
 */
SFCGAL_API auto
generateRoof(const Polygon &footprint, const RoofParameters &params)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Generate a roof using unified parameters without validity check.
 * @param footprint The 2D polygon footprint of the building
 * @param params The parameters for roof generation
 * @param nvc A NoValidityCheck object
 * @return A unique_ptr to the generated roof geometry
 */
SFCGAL_API auto
generateRoof(const Polygon &footprint, const RoofParameters &params,
             NoValidityCheck &nvc) -> std::unique_ptr<Geometry>;

/**
 * @brief Gable roof generation (Surface only).
 * @param polygon The 2D polygon footprint
 * @param clippingHeight Optional maximum height for clipping (0.0 for auto)
 * @param slopeAngle The slope angle in degrees
 * @return A unique_ptr to the generated polyhedral surface
 */
SFCGAL_API auto
extrudeGableRoof(const Polygon &polygon, double clippingHeight,
                 double slopeAngle = 45.0)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Skillion roof generation (Surface only).
 * @param polygon The 2D polygon footprint
 * @param clippingHeight Optional maximum height for clipping (0.0 for auto)
 * @param primaryEdgeIndex Index of the edge where the slope starts
 * @param primaryAngle The slope angle in degrees
 * @param secondaryAngle The angle for other sides (usually 90.0 for skillion)
 * @return A unique_ptr to the generated polyhedral surface
 */
SFCGAL_API auto
extrudeSkillionRoof(const Polygon &polygon, double clippingHeight,
                    size_t primaryEdgeIndex = 0, double primaryAngle = 30.0,
                    double secondaryAngle = 90.0)
    -> std::unique_ptr<PolyhedralSurface>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_ROOFGENERATION_H_
