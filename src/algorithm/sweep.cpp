// Copyright (c) 2025-2026, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/algorithm/sweep.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/algorithm/isClosed.h>
#include <SFCGAL/detail/tools/Log.h>
#include <SFCGAL/numeric.h>

#include <CGAL/Surface_mesh.h>

#include <cmath>
#include <stdexcept>

namespace SFCGAL::algorithm {

using Surface_mesh_3 = CGAL::Surface_mesh<Kernel::Point_3>;

// Forward declarations of internal helpers
namespace {

// Numerical thresholds used throughout the sweep implementation.
// Mirrors the named constants in Chamfer.cpp for consistency.
constexpr double ZERO_SQ_LEN       = 1e-20; // squared-length below which a vector is zero
constexpr double NEAR_ZERO_SQ_LEN  = 1e-12; // squared-length for degenerate checks
constexpr double NEAR_PARALLEL_COS = 0.99;  // parallelism guard for helper axis
constexpr double TOLERANCE_ALIGNMENT = 1e-6; // threshold for normal vector alignment
constexpr double TOLERANCE_PERPENDICULAR = 1e-10; // threshold for fallback perpendicular vector

/**
 * @brief Extract 2D profile points from geometry
 * @return Vector of 2D points (X, Y) in profile coordinate system
 */
auto
extract_profile_points(const Geometry &profile) -> std::vector<Kernel::Point_2>
{
  std::vector<Kernel::Point_2> points;

  if (profile.geometryTypeId() == TYPE_LINESTRING) {
    const auto &linestring = profile.as<LineString>();
    size_t num_points = linestring.numPoints();
    if (num_points > 1 && linestring.pointN(0) == linestring.pointN(num_points - 1)) {
      num_points--;
    }
    for (size_t i = 0; i < num_points; ++i) {
      const Point &point = linestring.pointN(i);
      points.emplace_back(point.x(), point.y());
    }
  } else if (profile.geometryTypeId() == TYPE_POLYGON) {
    const auto &polygon = profile.as<Polygon>();
    const auto &ring    = polygon.exteriorRing();
    size_t      num_points = ring.numPoints();
    if (num_points > 0) {
      num_points--;
    }
    for (size_t i = 0; i < num_points; ++i) {
      const Point &point = ring.pointN(i);
      points.emplace_back(point.x(), point.y());
    }
  } else {
    throw std::invalid_argument("Profile cannot be of this Geometry Type. Must be LineString or Polygon.");
  }

  if (points.size() < 2) {
    throw std::invalid_argument("Profile must have at least 2 points, got: " +
                                std::to_string(points.size()));
  }

  return points;
}

/**
 * @brief Compute initial frame from tangent direction
 * Used as a starting point for RMF propagation.
 */
auto
compute_initial_frame(const Kernel::Vector_3 &tangent) -> Frame
{
  Frame frame;
  frame.tangent = normalizeVector(tangent);

  // Choose helper vector perpendicular to tangent
  // Try X-axis first, if too parallel try Y-axis
  Kernel::Vector_3 helper(1, 0, 0);
  if (std::abs(CGAL::to_double(frame.tangent * helper)) > NEAR_PARALLEL_COS) {
    helper = Kernel::Vector_3(0, 1, 0);
  }

  // Gram-Schmidt orthogonalization
  frame.normal   = normalizeVector(CGAL::cross_product(frame.tangent, helper));
  frame.binormal = CGAL::cross_product(frame.tangent, frame.normal);

  return frame;
}

/**
 * @brief Propagate frame using Double Reflection Method (RMF)
 *
 * Implements the minimal rotation frame propagation described in:
 * Wang, W., Jüttler, B., Zheng, D., & Liu, Y. (2008).
 * "Computation of rotation minimizing frames."
 * ACM Transactions on Graphics (TOG), 27(1), 1-18.
 *
 * This method reflects the frame twice:
 * 1. Across the plane bisecting the chord (P_i, P_{i+1})
 * 2. Across the plane bisecting the tangents (T_i, T_{i+1})
 */
auto
propagate_frame_rmf(const Frame &prev_frame, const Kernel::Point_3 &prev_point,
                    const Kernel::Point_3  &curr_point,
                    const Kernel::Vector_3 &curr_tangent) -> Frame
{
  Frame frame;
  frame.tangent = normalizeVector(curr_tangent);

  // Vector from previous to current point (Chord)
  const Kernel::Vector_3 chord_vec    = curr_point - prev_point;
  const Kernel::FT       chord_len_sq = chord_vec * chord_vec;

  if (chord_len_sq < Kernel::FT(ZERO_SQ_LEN)) {
    // Points are coincident, keep same frame
    frame.normal   = prev_frame.normal;
    frame.binormal = prev_frame.binormal;
    return frame;
  }

  // First reflection: reflect normal across plane perpendicular to chord_vec
  // Formula 11 from Wang et al. (R_1)
  // R = P - 2(P.n)n where n is unit normal to reflection plane.
  // Here reflection plane normal is chord_vec/|chord_vec|.
  // R_1(n) = n - 2(n.chord_vec/|chord_vec|^2)chord_vec
  const Kernel::FT       reflection_factor_1 = Kernel::FT(2) / chord_len_sq;
  const Kernel::Vector_3 reflected_normal =
      prev_frame.normal - reflection_factor_1 * (chord_vec * prev_frame.normal) * chord_vec;
  const Kernel::Vector_3 reflected_tangent =
      prev_frame.tangent - reflection_factor_1 * (chord_vec * prev_frame.tangent) * chord_vec;

  // Second reflection to align with new tangent
  // Reflection plane normal is bisector of (reflected_tangent, new_tangent)
  const Kernel::Vector_3 tangent_diff    = frame.tangent - reflected_tangent;
  const Kernel::FT       tangent_diff_sq = tangent_diff * tangent_diff;

  if (tangent_diff_sq < Kernel::FT(ZERO_SQ_LEN)) {
    // Tangents are parallel, use first reflection result
    frame.normal = reflected_normal;
  } else {
    // Formula 11 again (R_2)
    const Kernel::FT reflection_factor_2 = Kernel::FT(2) / tangent_diff_sq;
    frame.normal        = reflected_normal - reflection_factor_2 * (tangent_diff * reflected_normal) * tangent_diff;
  }

  frame.normal   = normalizeVector(frame.normal);
  frame.binormal = CGAL::cross_product(frame.tangent, frame.normal);

  return frame;
}

/**
 * @brief Correct holonomy for closed paths
 * Distributes the closing error angle smoothly across all frames.
 */
void
correct_holonomy(std::vector<Frame> &frames)
{
  if (frames.size() < 2) {
    return;
  }

  const Frame &first_frame = frames.front();
  const Frame &last_frame  = frames.back();

  // Compute angle between first and last normal vectors
  const Kernel::FT dot_product = first_frame.normal * last_frame.normal;
  const double     cos_angle   = CGAL::to_double(dot_product);

  // If frames are already aligned, no correction needed
  if (SFCGAL::almostEqual(cos_angle, 1.0, TOLERANCE_ALIGNMENT)) {
    return;
  }

  // Compute total angle to distribute
  double total_angle = std::acos(std::clamp(cos_angle, -1.0, 1.0));

  // Determine rotation direction
  const Kernel::Vector_3 cross =
      CGAL::cross_product(last_frame.normal, first_frame.normal);
  if (cross * first_frame.tangent < Kernel::FT(0)) {
    total_angle = -total_angle;
  }

  // Apply incremental rotation to each frame
  const size_t n_frames = frames.size();
  for (size_t i = 1; i < n_frames; ++i) {
    const double fraction =
        static_cast<double>(i) / static_cast<double>(n_frames - 1);
    const double angle = total_angle * fraction;

    const double cos_a = std::cos(angle);
    const double sin_a = std::sin(angle);

    const Kernel::Vector_3 original_normal   = frames[i].normal;
    const Kernel::Vector_3 original_binormal = frames[i].binormal;

    frames[i].normal   = cos_a * original_normal + sin_a * original_binormal;
    frames[i].binormal = -sin_a * original_normal + cos_a * original_binormal;

    frames[i].normal   = normalizeVector(frames[i].normal);
    frames[i].binormal = normalizeVector(frames[i].binormal);
  }
}

/**
 * @brief Compute RMF frames along path (Discrete/Continuous hybrid)
 */
auto
compute_rmf_frames(const std::vector<Kernel::Point_3> &points, bool closed)
    -> std::vector<Frame>
{
  std::vector<Frame> frames;
  if (points.size() < 2) {
    return frames;
  }

  // Initialize first frame
  Kernel::Vector_3 first_tangent = points[1] - points[0];
  frames.push_back(compute_initial_frame(first_tangent));

  // Propagate frames along the path
  for (size_t i = 1; i < points.size(); ++i) {
    Kernel::Vector_3 tangent;
    if (i < points.size() - 1) {
      tangent = points[i + 1] - points[i];
    } else {
      tangent = points[i] - points[i - 1];
    }

    frames.push_back(
        propagate_frame_rmf(frames.back(), points[i - 1], points[i], tangent));
  }

  // Apply holonomy correction for closed paths
  if (closed) {
    correct_holonomy(frames);
  }

  return frames;
}


/**
 * @brief Compute bisector plane at a corner (miter join)
 */
auto
compute_bisector_plane(const Kernel::Point_3 &prev_point, const Kernel::Point_3 &curr_point,
                       const Kernel::Point_3 &next_point) -> Kernel::Plane_3
{
  const Kernel::Vector_3 vec_in   = normalizeVector(curr_point - prev_point);
  const Kernel::Vector_3 vec_out  = normalizeVector(next_point - curr_point);
  const Kernel::Vector_3 bisector = vec_in + vec_out; // Bisector direction

  // If segments are collinear (180 deg) or backtracking (0 deg), handle
  // robustly
  if (bisector * bisector < Kernel::FT(NEAR_ZERO_SQ_LEN)) {
    // Degenerate case: plane passes through curr_point perpendicular to vec_in
    return {curr_point, vec_in};
  }

  return {curr_point, bisector};
}

/**
 * @brief Project point onto plane along a direction (ray-plane intersection)
 */
auto
project_onto_bisector_plane(const Kernel::Point_3 &ray_origin,
                            const Kernel::Point_3 &ray_target,
                            const Kernel::Plane_3 &plane) -> Kernel::Point_3
{
  const Kernel::Vector_3 ray_dir = ray_target - ray_origin;

  const Kernel::FT numerator = -(plane.a() * ray_origin.x() + plane.b() * ray_origin.y() +
                                 plane.c() * ray_origin.z() + plane.d());
  const Kernel::FT denominator =
      plane.a() * ray_dir.x() + plane.b() * ray_dir.y() + plane.c() * ray_dir.z();

  if (CGAL::abs(denominator) < Kernel::FT(NEAR_ZERO_SQ_LEN)) {
    return ray_target; // Ray parallel to plane — keep original position
  }

  const Kernel::FT intersection_param = numerator / denominator;
  return ray_origin + intersection_param * ray_dir;
}

/**
 * @brief Compute frame perpendicular to segment axis (for Segment Aligned)
 */
auto
compute_segment_frame(
    const Kernel::Vector_3                &axis,
    const std::optional<Kernel::Vector_3> &ref_normal = std::nullopt) -> Frame
{
  Frame                  frame;
  const Kernel::Vector_3 normalized_axis = normalizeVector(axis);

  frame.tangent = normalized_axis;

  Kernel::Vector_3 perpendicular;

  if (ref_normal.has_value()) {
    // Use reference normal: project it onto plane perpendicular to tangent
    perpendicular = *ref_normal - (*ref_normal * frame.tangent) * frame.tangent;
  }

  // Fallback if no reference or degenerate reference
  if (!ref_normal.has_value() ||
      perpendicular * perpendicular < Kernel::FT(NEAR_ZERO_SQ_LEN)) {
    // Try Z-axis first, if parallel try Y-axis (original method)
    perpendicular =
        CGAL::cross_product(normalized_axis, Kernel::Vector_3(0, 0, 1));
    if (perpendicular * perpendicular < Kernel::FT(TOLERANCE_PERPENDICULAR)) {
      perpendicular =
          CGAL::cross_product(normalized_axis, Kernel::Vector_3(0, 1, 0));
    }
  }

  frame.normal = normalizeVector(perpendicular);
  frame.binormal =
      normalizeVector(CGAL::cross_product(frame.tangent, frame.normal));

  return frame;
}

/**
 * @brief Transform 2D profile point to 3D using frame
 */
auto
transform_profile_point(const Kernel::Point_2 &profile_pt,
                        const Kernel::Point_3 &path_pt, const Frame &frame,
                        double anchor_x, double anchor_y) -> Kernel::Point_3
{
  const double           profile_x = CGAL::to_double(profile_pt.x()) - anchor_x;
  const double           profile_y = CGAL::to_double(profile_pt.y()) - anchor_y;
  const Kernel::Vector_3 offset    = profile_x * frame.normal + profile_y * frame.binormal;
  return path_pt + offset;
}

/**
 * @brief Connect consecutive profile rings with quad faces.
 * @pre All entries in @p profiles must have the same number of points.
 */
auto
build_sweep_mesh(Surface_mesh_3                                  &mesh,
                 const std::vector<std::vector<Kernel::Point_3>> &profiles)
    -> std::vector<std::vector<Surface_mesh_3::Vertex_index>>
{
  std::vector<std::vector<Surface_mesh_3::Vertex_index>> vertex_rings;

  if (profiles.size() < 2) {
    return vertex_rings;
  }

  size_t n_path_points    = profiles.size();
  size_t n_profile_points = profiles[0].size();

  vertex_rings.reserve(n_path_points);

  for (const auto &profile : profiles) {
    std::vector<Surface_mesh_3::Vertex_index> ring;
    ring.reserve(profile.size());
    for (const auto &point : profile) {
      ring.push_back(mesh.add_vertex(point));
    }
    vertex_rings.push_back(std::move(ring));
  }

  for (size_t i = 0; i < n_path_points - 1; ++i) {
    const auto &ring1 = vertex_rings[i];
    const auto &ring2 = vertex_rings[i + 1];

    for (size_t j = 0; j < n_profile_points; ++j) {
      size_t next_j = (j + 1) % n_profile_points;
      mesh.add_face(ring1[j], ring2[j], ring2[next_j], ring1[next_j]);
    }
  }

  return vertex_rings;
}

/**
 * @brief Add flat end caps to close the tube.
 *
 * Tries reversed winding first, falls back to direct order if rejected
 * by the mesh (winding depends on path direction vs profile orientation).
 */
void
add_flat_caps(
    Surface_mesh_3                                               &mesh,
    const std::vector<std::vector<Surface_mesh_3::Vertex_index>> &vertex_rings,
    bool add_start, bool add_end)
{
  if (vertex_rings.empty()) {
    return;
  }

  const auto &first_ring = vertex_rings.front();
  const auto &last_ring  = vertex_rings.back();

  // Caps close the tube ends. The winding must be opposite to the
  // lateral faces' edge orientation on the boundary.
  if (add_start) {
    // Try reversed first, fall back to direct order
    std::vector<Surface_mesh_3::Vertex_index> cap_rev(first_ring.rbegin(),
                                                      first_ring.rend());
    auto                                      face_idx = mesh.add_face(cap_rev);
    if (face_idx == Surface_mesh_3::null_face()) {
      face_idx = mesh.add_face(first_ring);
    }
    if (face_idx == Surface_mesh_3::null_face()) {
      SFCGAL_WARNING("Failed to add start cap face");
    }
  }

  if (add_end) {
    auto face_idx = mesh.add_face(last_ring);
    if (face_idx == Surface_mesh_3::null_face()) {
      std::vector<Surface_mesh_3::Vertex_index> cap_rev(last_ring.rbegin(),
                                                        last_ring.rend());
      face_idx = mesh.add_face(cap_rev);
    }
    if (face_idx == Surface_mesh_3::null_face()) {
      SFCGAL_WARNING("Failed to add end cap face");
    }
  }
}

// ----------------------------------------------------------------------------
// Discrete Sweep Implementation (Segment Aligned)
// ----------------------------------------------------------------------------

/**
 * @brief Discrete sweep with miter joins at corners.
 *
 * Each path segment gets a constant frame. At corners, profile vertices
 * are projected onto the bisector plane to create clean miter joins.
 * Best for rectilinear paths (beams, walls, CAD extrusions).
 */
// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
sweep_discrete(const std::vector<Kernel::Point_3> &path_points,
               const Geometry &profile, const SweepOptions &options)
    -> std::unique_ptr<PolyhedralSurface>
{
  Surface_mesh_3                                         mesh;
  std::vector<std::vector<Surface_mesh_3::Vertex_index>> cap_rings;
  std::vector<Kernel::Point_2> profile_points = extract_profile_points(profile);

  const bool closed = options.closed_path;

  // 1. Pre-compute bisector planes
  std::vector<Kernel::Plane_3> bisector_planes;
  if (path_points.size() >= 3) {
    bisector_planes.reserve(path_points.size() - 2);
    for (size_t i = 1; i < path_points.size() - 1; ++i) {
      bisector_planes.push_back(compute_bisector_plane(
          path_points[i - 1], path_points[i], path_points[i + 1]));
    }

    if (closed) {
      bisector_planes.push_back(compute_bisector_plane(
          path_points[path_points.size() - 2],
          path_points[path_points.size() - 1], path_points[1]));
    }
  }

  // Previous segment's end ring - reused as next start ring
  std::vector<Surface_mesh_3::Vertex_index> prev_end_ring;
  std::vector<Surface_mesh_3::Vertex_index> first_start_ring;

  // Process each segment
  for (size_t seg_idx = 0; seg_idx < path_points.size() - 1; ++seg_idx) {
    // Constant frame for this segment
    const Kernel::Vector_3 axis          = path_points[seg_idx + 1] - path_points[seg_idx];
    const Frame            segment_frame = compute_segment_frame(axis, options.reference_normal);

    // --- Start Ring ---
    std::vector<Surface_mesh_3::Vertex_index> start_ring;

    if (seg_idx == 0) {
      // Create first ring from scratch
      std::vector<Kernel::Point_3> start_profile;
      start_profile.reserve(profile_points.size());

      if (closed && !bisector_planes.empty()) {
        // For closed path: the first ring must lie on the closing bisector
        // (miter join between the last segment and the first segment).
        // Project natural last-segment positions at path[n-2] through natural
        // last-segment positions at path[0] onto the closing bisector plane.
        // This follows the same segment-direction ray convention used at
        // every other interior corner.
        Kernel::Vector_3 last_axis =
            path_points[0] - path_points[path_points.size() - 2];
        Frame           last_frame = compute_segment_frame(last_axis, options.reference_normal);
        const Kernel::Plane_3 &closing_bisector = bisector_planes.back();

        for (const auto &profile_pt : profile_points) {
          Kernel::Point_3 p_prev = transform_profile_point(
              profile_pt, path_points[path_points.size() - 2], last_frame,
              options.anchor_x, options.anchor_y);
          Kernel::Point_3 p_curr = transform_profile_point(
              profile_pt, path_points[0], last_frame, options.anchor_x,
              options.anchor_y);
          start_profile.push_back(
              project_onto_bisector_plane(p_prev, p_curr, closing_bisector));
        }
      } else {
        for (const auto &profile_pt : profile_points) {
          start_profile.push_back(transform_profile_point(
              profile_pt, path_points[seg_idx], segment_frame, options.anchor_x,
              options.anchor_y));
        }
      }

      start_ring.reserve(profile_points.size());
      for (const auto &pt : start_profile) {
        start_ring.push_back(mesh.add_vertex(pt));
      }
      first_start_ring = start_ring;
    } else {
      // Reuse previous
      start_ring = prev_end_ring;
    }

    // --- End Ring ---
    std::vector<Surface_mesh_3::Vertex_index> end_ring;

    if (closed && seg_idx == path_points.size() - 2) {
      // Closing segment: first_start_ring already lies on the closing bisector.
      // Simply loop back to close the tube topology.
      end_ring = first_start_ring;
    } else {
      // Normal segment
      std::vector<Kernel::Point_3> end_profile;
      end_profile.reserve(profile_points.size());

      for (const auto &profile_pt : profile_points) {
        end_profile.push_back(transform_profile_point(
            profile_pt, path_points[seg_idx + 1], segment_frame,
            options.anchor_x, options.anchor_y));
      }

      end_ring.reserve(profile_points.size());
      for (size_t pt_idx = 0; pt_idx < profile_points.size(); ++pt_idx) {
        Kernel::Point_3 end_point = end_profile[pt_idx];

        // If there is a corner ahead, project onto bisector plane along
        // the ray from start to end (preserves correct miter geometry)
        if (seg_idx < path_points.size() - 2) {
          const Kernel::Point_3 &start_point = mesh.point(start_ring[pt_idx]);
          end_point = project_onto_bisector_plane(start_point, end_point,
                                                  bisector_planes[seg_idx]);
        }
        end_ring.push_back(mesh.add_vertex(end_point));
      }
    }

    prev_end_ring = end_ring;

    // Connect faces between start and end rings.
    // Use triangles instead of quads to avoid non-planar face issues
    // at miter join corners (where bisector projection shifts vertices
    // by different amounts).
    const size_t n_pts = profile_points.size();
    for (size_t i = 0; i < n_pts - 1; ++i) {
      mesh.add_face(start_ring[i], start_ring[i + 1], end_ring[i + 1]);
      mesh.add_face(start_ring[i], end_ring[i + 1], end_ring[i]);
    }
    // Close the last profile edge (last vertex back to first)
    if (profile.geometryTypeId() == TYPE_POLYGON ||
        (profile.geometryTypeId() == TYPE_LINESTRING &&
         profile.as<LineString>().isClosed())) {
      mesh.add_face(start_ring[n_pts - 1], start_ring[0], end_ring[0]);
      mesh.add_face(start_ring[n_pts - 1], end_ring[0], end_ring[n_pts - 1]);
    }

    if (seg_idx == 0) {
      cap_rings.push_back(start_ring);
    }
    if (seg_idx == path_points.size() - 2) {
      cap_rings.push_back(end_ring);
    }
  }

  if (!closed && (options.start_cap == SweepOptions::EndCapStyle::FLAT ||
                  options.end_cap == SweepOptions::EndCapStyle::FLAT)) {
    const bool add_start = (options.start_cap == SweepOptions::EndCapStyle::FLAT);
    const bool add_end   = (options.end_cap == SweepOptions::EndCapStyle::FLAT);
    add_flat_caps(mesh, cap_rings, add_start, add_end);
  }

  return std::make_unique<PolyhedralSurface>(mesh);
}
// NOLINTEND(readability-function-cognitive-complexity)

// ----------------------------------------------------------------------------
// Continuous Sweep Implementation (Point Based)
// ----------------------------------------------------------------------------

/**
 * @brief Continuous sweep for smooth curves (RMF / Frenet)
 */
auto
sweep_continuous(const std::vector<Kernel::Point_3> &path_points,
                 const Geometry &profile, const SweepOptions &options)
    -> std::unique_ptr<PolyhedralSurface>
{
  std::vector<Kernel::Point_2> profile_points = extract_profile_points(profile);
  const bool                   closed         = options.closed_path;

  // 1. Compute Frames
  std::vector<Frame> frames;
  if (options.frame_method == SweepOptions::FrameMethod::ROTATION_MINIMIZING) {
    frames = compute_rmf_frames(path_points, closed);
  }

  // 2. Transform profile at each point
  std::vector<std::vector<Kernel::Point_3>> transformed_profiles;
  transformed_profiles.reserve(path_points.size());

  for (size_t i = 0; i < path_points.size(); ++i) {
    std::vector<Kernel::Point_3> profile_3d;
    profile_3d.reserve(profile_points.size());

    for (const auto &profile_pt : profile_points) {
      profile_3d.push_back(transform_profile_point(profile_pt, path_points[i],
                                                   frames[i], options.anchor_x,
                                                   options.anchor_y));
    }
    transformed_profiles.push_back(std::move(profile_3d));
  }

  // 3. Build Mesh
  Surface_mesh_3 mesh;
  auto           vertex_rings = build_sweep_mesh(mesh, transformed_profiles);

  if (!closed && (options.start_cap == SweepOptions::EndCapStyle::FLAT ||
                  options.end_cap == SweepOptions::EndCapStyle::FLAT)) {
    bool add_start = (options.start_cap == SweepOptions::EndCapStyle::FLAT);
    bool add_end   = (options.end_cap == SweepOptions::EndCapStyle::FLAT);
    add_flat_caps(mesh, vertex_rings, add_start, add_end);
  }

  return std::make_unique<PolyhedralSurface>(mesh);
}

} // anonymous namespace

// ----------------------------------------------------------------------------
// Public API Implementation
// ----------------------------------------------------------------------------

auto
sweep(const LineString &path, const Geometry &profile,
      const SweepOptions &options) -> std::unique_ptr<PolyhedralSurface>
{
  if (path.numPoints() < 2) {
    throw std::invalid_argument("Path must have at least 2 points");
  }

  std::vector<Kernel::Point_3> path_points;
  path_points.reserve(path.numPoints());
  for (size_t i = 0; i < path.numPoints(); ++i) {
    const Point &point = path.pointN(i);
    path_points.emplace_back(point.x(), point.y(), point.z());
  }

  // Determine if path is closed (explicit option or geometric closure)
  bool closed = options.closed_path || isClosed(path);

  // If closed but not geometrically closed, duplicate start point
  if (closed && !path_points.empty()) {
    if (path_points.front() != path_points.back()) {
      path_points.push_back(path_points.front());
    }
  }

  // Build effective options with resolved closed flag
  SweepOptions opts = options;
  opts.closed_path  = closed;

  if (opts.frame_method == SweepOptions::FrameMethod::SEGMENT_ALIGNED) {
    return sweep_discrete(path_points, profile, opts);
  }
  return sweep_continuous(path_points, profile, opts);
}

auto
create_circular_profile(double radius, int segments) -> std::unique_ptr<Polygon>
{
  if (radius <= 0.0) {
    throw std::invalid_argument("Radius must be positive");
  }
  if (segments < 3) {
    throw std::invalid_argument("Segments must be at least 3");
  }

  std::vector<Point> points;
  points.reserve(segments + 1);

  for (int i = 0; i < segments; ++i) {
    const double angle = 2.0 * M_PI * i / segments;
    const double x     = radius * std::cos(angle);
    const double y     = radius * std::sin(angle);
    points.emplace_back(x, y, 0);
  }
  points.emplace_back(points[0]); // Close

  LineString ring(points);
  return std::make_unique<Polygon>(ring);
}

auto
create_rectangular_profile(double width, double height)
    -> std::unique_ptr<Polygon>
{
  if (width <= 0.0 || height <= 0.0) {
    throw std::invalid_argument("Dimensions must be positive");
  }

  const double half_width  = width / 2.0;
  const double half_height = height / 2.0;

  std::vector<Point> points;
  points.reserve(5);
  points.emplace_back(-half_width, -half_height, 0);
  points.emplace_back(half_width, -half_height, 0);
  points.emplace_back(half_width, half_height, 0);
  points.emplace_back(-half_width, half_height, 0);
  points.emplace_back(-half_width, -half_height, 0);

  LineString ring(points);
  return std::make_unique<Polygon>(ring);
}

} // namespace SFCGAL::algorithm
