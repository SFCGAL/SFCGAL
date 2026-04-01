// Copyright (c) 2025-2026, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_SWEEP_H_
#define SFCGAL_ALGORITHM_SWEEP_H_

#include <SFCGAL/config.h>

#include <SFCGAL/Geometry.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>

#include <memory>
#include <optional>
#include <vector>

namespace SFCGAL::algorithm {

/**
 * @brief Options for controlling sweep operation behavior
 */
struct SweepOptions {
  /**
   * @brief Frame computation method for sweep operation
   */
  enum class FrameMethod : std::uint8_t {
    /**
     * @brief Rotation Minimizing Frame (RMF) / Double Reflection Method.
     *
     * Best for smooth, organic paths (pipes, cables, knots).
     * Computes frames that minimize twist around the tangent.
     * Based on: Wang, W., et al. (2008). "Computation of rotation minimizing
     * frames". ACM Transactions on Graphics, 27(1), 1-18.
     */
    ROTATION_MINIMIZING,

    /**
     * @brief Frenet-Serret Frame.
     *
     * Classic differential geometry frame defined by tangent and curvature.
     * Warning: Undefined on straight lines and flips at inflection points.
     * Use only if specifically required by mathematical definition.
     */
    FRENET,

    /**
     * @brief Segment Aligned (Discrete Miter).
     *
     * Best for architecture, mechanical parts, and CAD (beams, walls).
     * Does not use continuous frame propagation. Instead, computes a constant
     * orientation for each segment and projects vertices onto bisector planes
     * at corners to create perfect miter joins.
     */
    SEGMENT_ALIGNED
  };

  /**
   * @brief End cap style for sweep termination
   */
  enum class EndCapStyle : std::uint8_t {
    NONE, ///< Tube remains open at ends
    FLAT, ///< Tube is closed with flat faces
  };

  /// Method used to compute the reference frame along the path
  FrameMethod frame_method = FrameMethod::ROTATION_MINIMIZING;

  /// Style for the start of the sweep
  EndCapStyle start_cap = EndCapStyle::FLAT;

  /// Style for the end of the sweep
  EndCapStyle end_cap = EndCapStyle::FLAT;

  /// If true, connects the end of the path back to the start
  bool closed_path = false;

  /**
   * @brief Anchor point X coordinate in profile space (default: 0.0)
   * The anchor point of the profile is the point that will travel exactly
   * along the path.
   */
  double anchor_x = 0.0;

  /**
   * @brief Anchor point Y coordinate in profile space (default: 0.0)
   */
  double anchor_y = 0.0;

  /**
   * @brief Optional reference normal to orient the frame
   *
   * If provided, the frame's normal will be aligned with this vector
   * (projected onto the plane perpendicular to the path tangent).
   * This ensures consistent orientation for chamfer/fillet operations
   * on arbitrarily oriented solids.
   */
  std::optional<Kernel::Vector_3> reference_normal = std::nullopt;
};

/**
 * @brief Frame structure for sweep operations
 * Represents a local coordinate system (Coordinate Frame) at a point on the
 * path.
 */
struct Frame {
  Kernel::Vector_3 tangent;  ///< Z-axis of local frame (direction of motion)
  Kernel::Vector_3 normal;   ///< X-axis of local frame
  Kernel::Vector_3 binormal; ///< Y-axis of local frame
};

/**
 * @brief Sweep a 2D profile along a 3D path to create a 3D surface
 *
 * This function extrudes a 2D profile (polygon or linestring) along a 3D path
 * to create a polyhedral surface. The profile is defined in the XY plane where:
 * - X corresponds to the Normal direction of the path frame
 * - Y corresponds to the Binormal direction of the path frame
 *
 * @param path 3D LineString defining the sweep path
 * @param profile 2D Geometry (LineString or Polygon) defining the cross-section
 * @param options SweepOptions controlling frame method, end caps, etc.
 * @return std::unique_ptr<PolyhedralSurface> The resulting 3D surface
 * @throws std::invalid_argument if inputs are invalid
 *
 * @example
 * // Sweep with default options (RMF, flat caps, origin on path)
 * auto surface1 = sweep(helix, circle);
 *
 * // Sweep with anchor at rectangle center (1, 0.5)
 * SweepOptions opts;
 * opts.anchor_x = 1.0;
 * opts.anchor_y = 0.5;
 * auto surface2 = sweep(helix, rect, opts);
 */
SFCGAL_API auto
sweep(const LineString &path, const Geometry &profile,
      const SweepOptions &options = {}) -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Create a circular profile for sweep operations
 * Generates a Polygon approximating a circle centered at (0,0).
 * @param radius Radius of the circle
 * @param segments Number of segments to approximate the circle
 * @return std::unique_ptr<LineString> (Closed LineString)
 */
SFCGAL_API auto
create_circular_profile(double radius, int segments = 32)
    -> std::unique_ptr<LineString>;

/**
 * @brief Create a rectangular profile for sweep operations
 * Generates a rectangular Polygon centered at (0,0).
 * @param width Width (X dimension)
 * @param height Height (Y dimension)
 * @return std::unique_ptr<Polygon>
 */
SFCGAL_API auto
create_rectangular_profile(double width, double height)
    -> std::unique_ptr<Polygon>;

/**
 * @brief Create a triangular profile for chamfer operations
 *
 * Generates a right-angled triangle in the 3rd quadrant, intended for
 * subtracting material (chamfer) from a 90-degree corner.
 *
 * The triangle vertices are:
 * (0, 0)       - The corner
 * (0, -radius_y)
 * (-radius_x, 0)
 *
 * @param radius_x Length of the chamfer along X axis
 * @param radius_y Length of the chamfer along Y axis (if < 0, uses radius_x)
 * @return std::unique_ptr<Polygon>
 */
SFCGAL_API auto
create_chamfer_profile(double radius_x, double radius_y = -1.0)
    -> std::unique_ptr<Polygon>;

/**
 * @brief Create a fillet (rounded) profile for rounding operations
 *
 * Generates a profile representing the material to be removed to create a
 * rounded corner (fillet) of a specific radius. The shape is located in the
 * 3rd quadrant.
 *
 * The shape connects:
 * (0, 0) -> (0, -radius) -> arc to (-radius, 0) -> (0, 0)
 *
 * @param radius Radius of the fillet
 * @param segments Number of segments to approximate the 90-degree arc
 * @return std::unique_ptr<Polygon>
 */
SFCGAL_API auto
create_fillet_profile(double radius, int segments = 4)
    -> std::unique_ptr<Polygon>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_SWEEP_H_
