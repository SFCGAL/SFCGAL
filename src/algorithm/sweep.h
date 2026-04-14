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
   * @brief Optional reference vector to orient the profile consistently.
   *
   * If provided, the frame's normal direction is aligned with this vector
   * (projected onto the plane perpendicular to the path tangent).
   * Used internally by chamfer to align the cutting profile with face normals.
   * Only affects SEGMENT_ALIGNED frame method.
   */
  std::optional<Kernel::Vector_3> reference_normal = std::nullopt;
};

/**
 * @brief Right-handed orthonormal basis at a point on the sweep path.
 *
 * Defines how the 2D profile is oriented in 3D space at each path vertex.
 * The profile's X-axis maps to @c normal and its Y-axis maps to @c binormal.
 * The three vectors satisfy: binormal = tangent × normal.
 */
struct Frame {
  Kernel::Vector_3 tangent;  ///< Direction of motion along the path
  Kernel::Vector_3 normal;   ///< Profile X-axis direction
  Kernel::Vector_3 binormal; ///< Profile Y-axis direction (tangent × normal)
};

/**
 * @brief Sweep a 2D profile along a 3D path to create a polyhedral surface.
 *
 * Extrudes a 2D cross-section (the @p profile) along a 3D curve (the @p path).
 * At each path vertex, the profile is placed in a local coordinate frame:
 * - Profile X → frame Normal direction
 * - Profile Y → frame Binormal direction
 *
 * The anchor point (anchor_x, anchor_y) controls which point of the profile
 * rides exactly along the path (default: origin).
 *
 * @par Limitations
 * - Profile must be a Polygon or LineString; only exterior ring is used.
 * - Only X and Y coordinates of the profile are used; Z is discarded.
 * - FRENET method is undefined on straight lines; prefer ROTATION_MINIMIZING
 *   or SEGMENT_ALIGNED.
 * - End caps are only added on open (non-closed) paths.
 *
 * @param path 3D LineString (>= 2 points).
 * @param profile 2D Polygon or LineString cross-section.
 * @param options Frame method, end caps, anchor, closure.
 * @return The resulting 3D polyhedral surface.
 *
 * @throws std::invalid_argument If path has < 2 points.
 * @throws std::invalid_argument If profile is not Polygon or LineString.
 * @throws std::invalid_argument If profile has < 2 distinct points.
 */
SFCGAL_API auto
sweep(const LineString &path, const Geometry &profile,
      const SweepOptions &options = {}) -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Create a circular profile for sweep operations
 * Generates a Polygon approximating a circle centered at (0,0).
 * @param radius Radius of the circle
 * @param segments Number of segments to approximate the circle
 * @return std::unique_ptr<Polygon>
 * @throws std::invalid_argument If radius <= 0 or segments < 3.
 */
SFCGAL_API auto
create_circular_profile(double radius, int segments = 32)
    -> std::unique_ptr<Polygon>;

/**
 * @brief Create a rectangular profile for sweep operations
 * Generates a rectangular Polygon centered at (0,0).
 * @param width Width (X dimension)
 * @param height Height (Y dimension)
 * @return std::unique_ptr<Polygon>
 * @throws std::invalid_argument If width <= 0 or height <= 0.
 */
SFCGAL_API auto
create_rectangular_profile(double width, double height)
    -> std::unique_ptr<Polygon>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_SWEEP_H_
