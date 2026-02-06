// Copyright (c) 2025-2026, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_SWEEP_H_
#define SFCGAL_ALGORITHM_SWEEP_H_

#include <SFCGAL/config.h>

#include <SFCGAL/Geometry.h>
#include <SFCGAL/Kernel.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>

#include <cstdint>
#include <memory>

namespace SFCGAL::algorithm {

/**
 * @brief Options for controlling sweep operation behavior
 */
struct SweepOptions {
  /**
   * @brief Frame computation method for sweep operation
   */
  enum class FrameMethod : std::uint8_t {
    ROTATION_MINIMIZING, ///< Rotation Minimizing Frames (RMF) - minimal twist
    FRENET,              ///< Frenet-Serret frames - natural curvature-based
    FIXED_UP,            ///< Fixed up vector - constant reference direction
    SEGMENT_ALIGNED ///< Segment-aligned frames with miter joins (like buffer3D)
  };

  /**
   * @brief End cap style for sweep termination
   */
  enum class EndCapStyle : std::uint8_t {
    NONE, ///< No end caps (open tube)
    FLAT, ///< Flat planar caps perpendicular to path
    ROUND ///< Rounded hemispherical caps (future)
  };

  /// Frame computation method (default: RMF for minimal twist)
  FrameMethod frame_method = FrameMethod::ROTATION_MINIMIZING;

  /// Start cap style
  EndCapStyle start_cap = EndCapStyle::FLAT;

  /// End cap style
  EndCapStyle end_cap = EndCapStyle::FLAT;

  /// Whether path is closed (first and last points coincide)
  /// When true, applies holonomy correction for seamless closure
  bool closed_path = false;

  /// Fixed up vector for FIXED_UP frame method (default: Z-axis)
  Kernel::Vector_3 fixed_up_vector = Kernel::Vector_3(0, 0, 1);

  /// Tolerance for closed path detection (squared distance)
  double closed_path_tolerance = 1e-10;

  /// Anchor point X coordinate in profile space (default: origin)
  /// The anchor point is positioned on the path during sweep
  double anchor_x = 0.0;

  /// Anchor point Y coordinate in profile space (default: origin)
  /// The anchor point is positioned on the path during sweep
  double anchor_y = 0.0;
};

/**
 * @brief Frame structure for sweep operations
 *
 * Represents an orthonormal frame at each point along the path:
 * - tangent: direction along the path (T)
 * - normal: first perpendicular direction (U) - maps to profile X-axis
 * - binormal: second perpendicular direction (V) - maps to profile Y-axis
 *
 * Forms a right-handed coordinate system: binormal = tangent × normal
 */
struct Frame {
  Kernel::Vector_3 tangent;  ///< T: direction along path
  Kernel::Vector_3 normal;   ///< U: first perpendicular (profile X)
  Kernel::Vector_3 binormal; ///< V: second perpendicular (profile Y)
};

/**
 * @brief Sweep a 2D profile along a 3D path to create a 3D surface
 *
 * This function extrudes a 2D profile (polygon or linestring) along a 3D path
 * to create a polyhedral surface. The profile is defined in the XY plane where:
 * - X-axis maps to the frame normal direction (U)
 * - Y-axis maps to the frame binormal direction (V)
 * - The anchor point (specified in options) is positioned at the path center
 *
 * @param path 3D LineString defining the sweep path
 * @param profile 2D Polygon or LineString defining the cross-section
 *                Must be 2D (Z coordinates ignored if present)
 * @param options SweepOptions controlling frame method, end caps, etc.
 *
 * @return PolyhedralSurface representing the swept geometry
 *
 * @throws std::invalid_argument if path has < 2 points
 * @throws std::invalid_argument if profile has < 2 points
 * @throws std::invalid_argument if profile is not 2D
 *
 * @note Profile orientation:
 *       - For closed profiles (Polygon): exterior ring swept as tube wall,
 *         interior rings create holes in the wall
 *       - For open profiles (LineString): creates a ribbon/blade surface
 *
 * @note Anchor point:
 *       The anchor point (anchor_x, anchor_y) specifies which point in the
 *       profile coordinate system is positioned on the path. By default
 *       (anchor_x=0, anchor_y=0), the profile origin follows the path.
 *       For example, to center a rectangle on the path, set the anchor
 *       to the rectangle's center coordinates.
 *
 * @example
 * @code
 * // Example 1: Circular profile along helical path
 * std::vector<Point> circle_pts;
 * for (int i = 0; i < 16; ++i) {
 *   double angle = 2.0 * M_PI * i / 16.0;
 *   circle_pts.emplace_back(cos(angle), sin(angle), 0);
 * }
 * circle_pts.emplace_back(0, 0, 0)
 * LineString circle(circle_pts);
 *
 * std::vector<Point> helix_pts;
 * for (int i = 0; i <= 20; ++i) {
 *   double t = i / 20.0;
 *   double angle = 4.0 * M_PI * t;
 *   helix_pts.emplace_back(2*cos(angle), 2*sin(angle), 10*t);
 * }
 * LineString helix(helix_pts);
 *
 * // Sweep with default options (RMF, flat caps, origin on path)
 * auto surface1 = sweep(helix, circle);
 *
 * // Example 2: Rectangular profile with custom anchor point
 * std::vector<Point> rect_pts;
 * rect_pts.emplace_back(0, 0, 0);
 * rect_pts.emplace_back(2, 0, 0);
 * rect_pts.emplace_back(2, 1, 0);
 * rect_pts.emplace_back(0, 1, 0);
 * rect_pts.emplace_back(0, 0, 0);
 * LineString rect(rect_pts);
 *
 * // Sweep with anchor at rectangle center (1, 0.5)
 * SweepOptions opts;
 * opts.anchor_x = 1.0;
 * opts.anchor_y = 0.5;
 * auto surface2 = sweep(helix, rect, opts);
 * @endcode
 */
SFCGAL_API auto
sweep(const LineString &path, const Geometry &profile,
      const SweepOptions &options = {}) -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Create a circular profile for sweep operations
 *
 * Helper function to generate a circular profile centered at origin.
 * Useful for creating pipe/tube geometries.
 *
 * @param radius Circle radius
 * @param segments Number of segments (vertices) in the circle
 * @return LineString representing the circular profile in XY plane
 */
SFCGAL_API auto
create_circular_profile(double radius, int segments)
    -> std::unique_ptr<LineString>;

/**
 * @brief Create a rectangular profile for sweep operations
 *
 * Helper function to generate a rectangular profile centered at origin.
 * Useful for creating beam/extrusion geometries.
 *
 * @param width Width of rectangle (X direction)
 * @param height Height of rectangle (Y direction)
 * @return Polygon representing the rectangular profile in XY plane
 */
SFCGAL_API auto
create_rectangular_profile(double width, double height)
    -> std::unique_ptr<Polygon>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_SWEEP_H_
