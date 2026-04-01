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
#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/detail/tools/Registry.h>

#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Surface_mesh.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace SFCGAL::algorithm {

using Surface_mesh_3 = CGAL::Surface_mesh<Kernel::Point_3>;

namespace PMP = CGAL::Polygon_mesh_processing;

// Forward declarations of internal helpers
namespace {

/**
 * @brief Normalize a vector to unit length
 */
auto
normalizeVector(const Kernel::Vector_3 &vector) -> Kernel::Vector_3
{
  Kernel::FT len_sq = vector * vector;
  if (len_sq < Kernel::FT(1e-20)) {
    return {0, 0, 0};
  }
  double len = std::sqrt(CGAL::to_double(len_sq));
  return vector / len;
}

/**
 * @brief Visitor to extract 2D profile points from geometry
 */
class ProfilePointExtractor : public ConstGeometryVisitor {
public:
  std::vector<Kernel::Point_2> points;

  void
  visit(const LineString &linestring) override
  {
    size_t num_points = linestring.numPoints();

    // Check if the linestring is closed (last point equals first point)
    // If so, skip the last point to avoid duplicates
    if (num_points > 1 &&
        linestring.pointN(0) == linestring.pointN(num_points - 1)) {
      num_points--;
    }

    for (size_t i = 0; i < num_points; ++i) {
      const Point &point = linestring.pointN(i);
      points.emplace_back(point.x(), point.y());
    }
  }

  void
  visit(const Polygon &polygon) override
  {
    const auto &ring       = polygon.exteriorRing();
    size_t      num_points = ring.numPoints();

    // Polygon rings are always closed (last point equals first)
    // Skip the last point to avoid duplicates
    if (num_points > 0) {
      num_points--;
    }

    for (size_t i = 0; i < num_points; ++i) {
      const Point &point = ring.pointN(i);
      points.emplace_back(point.x(), point.y());
    }
  }

  // Default implementations for unsupported types
  void
  visit(const Point & /*point*/) override
  {
    throw std::invalid_argument("Profile cannot be a Point");
  }
  void
  visit(const Triangle & /*triangle*/) override
  {
    throw std::invalid_argument("Profile cannot be a Triangle");
  }
  void
  visit(const Solid & /*solid*/) override
  {
    throw std::invalid_argument("Profile cannot be a Solid");
  }
  void
  visit(const MultiPoint & /*multipoint*/) override
  {
    throw std::invalid_argument("Profile cannot be a MultiPoint");
  }
  void
  visit(const MultiLineString & /*multilinestring*/) override
  {
    throw std::invalid_argument("Profile cannot be a MultiLineString");
  }
  void
  visit(const MultiPolygon & /*multipolygon*/) override
  {
    throw std::invalid_argument("Profile cannot be a MultiPolygon");
  }
  void
  visit(const MultiSolid & /*multisolid*/) override
  {
    throw std::invalid_argument("Profile cannot be a MultiSolid");
  }
  void
  visit(const GeometryCollection & /*collection*/) override
  {
    throw std::invalid_argument("Profile cannot be a GeometryCollection");
  }
  void
  visit(const PolyhedralSurface & /*surface*/) override
  {
    throw std::invalid_argument("Profile cannot be a PolyhedralSurface");
  }
  void
  visit(const TriangulatedSurface & /*surface*/) override
  {
    throw std::invalid_argument("Profile cannot be a TriangulatedSurface");
  }
};

/**
 * @brief Extract 2D profile points from geometry using visitor pattern
 * @return Vector of 2D points (X, Y) in profile coordinate system
 */
auto
extract_profile_points(const Geometry &profile) -> std::vector<Kernel::Point_2>
{
  ProfilePointExtractor extractor;
  profile.accept(extractor);

  if (extractor.points.size() < 2) {
    throw std::invalid_argument("Profile must have at least 2 points, got: " +
                                std::to_string(extractor.points.size()));
  }

  return extractor.points;
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
  if (std::abs(CGAL::to_double(frame.tangent * helper)) > 0.9) {
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
  Kernel::Vector_3 v1           = curr_point - prev_point;
  Kernel::FT       v1_length_sq = v1 * v1;

  if (v1_length_sq < Kernel::FT(1e-20)) {
    // Points are coincident, keep same frame
    frame.normal   = prev_frame.normal;
    frame.binormal = prev_frame.binormal;
    return frame;
  }

  // First reflection: reflect normal across plane perpendicular to v1
  // Formula 11 from Wang et al. (R_1)
  // R = P - 2(P.n)n where n is unit normal to reflection plane.
  // Here reflection plane normal is v1/|v1|.
  // R_1(n) = n - 2(n.v1/|v1|^2)v1
  Kernel::FT       c1 = Kernel::FT(2) / v1_length_sq;
  Kernel::Vector_3 reflected_normal =
      prev_frame.normal - c1 * (v1 * prev_frame.normal) * v1;
  Kernel::Vector_3 reflected_tangent =
      prev_frame.tangent - c1 * (v1 * prev_frame.tangent) * v1;

  // Second reflection to align with new tangent
  // Reflection plane normal is bisector of (reflected_tangent, new_tangent)
  Kernel::Vector_3 v2           = frame.tangent - reflected_tangent;
  Kernel::FT       v2_length_sq = v2 * v2;

  if (v2_length_sq < Kernel::FT(1e-20)) {
    // Tangents are parallel, use first reflection result
    frame.normal = reflected_normal;
  } else {
    // Formula 11 again (R_2)
    Kernel::FT c2 = Kernel::FT(2) / v2_length_sq;
    frame.normal  = reflected_normal - c2 * (v2 * reflected_normal) * v2;
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
correct_holonomy(std::vector<Frame> &frames,
                 const std::vector<Kernel::Point_3> & /*points*/)
{
  if (frames.size() < 2) {
    return;
  }

  const Frame &first_frame = frames.front();
  const Frame &last_frame  = frames.back();

  // Compute angle between first and last normal vectors
  Kernel::FT dot_product = first_frame.normal * last_frame.normal;
  double     cos_angle   = CGAL::to_double(dot_product);

  // If frames are already aligned, no correction needed
  if (std::abs(cos_angle - 1.0) < 1e-6) {
    return;
  }

  // Compute total angle to distribute
  double total_angle = std::acos(std::max(-1.0, std::min(1.0, cos_angle)));

  // Determine rotation direction
  Kernel::Vector_3 cross =
      CGAL::cross_product(last_frame.normal, first_frame.normal);
  if (cross * first_frame.tangent < Kernel::FT(0)) {
    total_angle = -total_angle;
  }

  // Apply incremental rotation to each frame
  size_t n_frames = frames.size();
  for (size_t i = 1; i < n_frames; ++i) {
    double fraction =
        static_cast<double>(i) / static_cast<double>(n_frames - 1);
    double angle = total_angle * fraction;

    double cos_a = std::cos(angle);
    double sin_a = std::sin(angle);

    Kernel::Vector_3 original_normal   = frames[i].normal;
    Kernel::Vector_3 original_binormal = frames[i].binormal;

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
    correct_holonomy(frames, points);
  }

  return frames;
}

/**
 * @brief Compute Frenet-Serret frames along path
 *
 * T(t) = P'(t) / |P'(t)|
 * B(t) = (P'(t) x P''(t)) / |P'(t) x P''(t)|
 * N(t) = B(t) x T(t)
 *
 * Warning: This method is unstable on straight lines (curvature = 0)
 * and at inflection points.
 */
auto
compute_frenet_frames(const std::vector<Kernel::Point_3> &points, bool closed)
    -> std::vector<Frame>
{
  (void)closed;
  std::vector<Frame> frames;
  if (points.size() < 2) {
    return frames;
  }

  for (size_t i = 0; i < points.size(); ++i) {
    Frame frame;

    // 1. Compute Tangent (Velocity)
    Kernel::Vector_3 tangent;
    if (i < points.size() - 1) {
      tangent = points[i + 1] - points[i];
    } else {
      tangent = points[i] - points[i - 1];
    }
    frame.tangent = normalizeVector(tangent);

    // 2. Compute Acceleration (Curvature direction)
    // We use central difference approximation for internal points
    Kernel::Vector_3 acceleration(0, 0, 0);

    if (i > 0 && i < points.size() - 1) {
      // Central difference: (P_{i+1} - 2P_i + P_{i-1})
      Kernel::Vector_3 v_in  = points[i] - points[i - 1];
      Kernel::Vector_3 v_out = points[i + 1] - points[i];
      acceleration           = v_out - v_in;
    } else if (i == 0 && points.size() > 2) {
      // Forward difference approximation
      // P'' ≈ P_2 - 2P_1 + P_0
      acceleration = (points[2] - points[1]) - (points[1] - points[0]);
    } else if (i == points.size() - 1 && points.size() > 2) {
      // Backward difference
      acceleration =
          (points[i] - points[i - 1]) - (points[i - 1] - points[i - 2]);
    }

    // 3. Compute Binormal = Tangent x Acceleration
    Kernel::Vector_3 binormal =
        CGAL::cross_product(frame.tangent, acceleration);
    Kernel::FT binormal_len_sq = binormal * binormal;

    if (binormal_len_sq > Kernel::FT(1e-10)) {
      // Normal curvature exists
      frame.binormal = normalizeVector(binormal);
      frame.normal   = CGAL::cross_product(frame.binormal, frame.tangent);
    } else {
      // Degenerate case (Straight line or inflection point)
      // Fallback:
      if (frames.empty()) {
        // If it's the very first point and it's straight, pick arbitrary
        frame = compute_initial_frame(frame.tangent);
      } else {
        // Propagate previous frame (Parallel transport approximation)
        Frame prev = frames.back();
        // Ensure orthogonality with new tangent
        // Project previous normal onto plane perpendicular to new tangent
        Kernel::Vector_3 proj_normal =
            prev.normal - (prev.normal * frame.tangent) * frame.tangent;

        if (proj_normal * proj_normal > Kernel::FT(1e-10)) {
          frame.normal   = normalizeVector(proj_normal);
          frame.binormal = CGAL::cross_product(frame.tangent, frame.normal);
        } else {
          // New tangent is parallel to old normal? Unlikely if curve is smooth.
          // Fallback to initial frame computation
          frame = compute_initial_frame(frame.tangent);
        }
      }
    }

    frames.push_back(frame);
  }

  return frames;
}

/**
 * @brief Compute bisector plane at a corner (miter join)
 */
auto
compute_bisector_plane(const Kernel::Point_3 &p1, const Kernel::Point_3 &p2,
                       const Kernel::Point_3 &p3) -> Kernel::Plane_3
{
  Kernel::Vector_3 v1       = normalizeVector(p2 - p1);
  Kernel::Vector_3 v2       = normalizeVector(p3 - p2);
  Kernel::Vector_3 bisector = v1 + v2; // Bisector direction

  // If segments are collinear (180 deg) or backtracking (0 deg), handle
  // robustly
  if (bisector * bisector < Kernel::FT(1e-10)) {
    // Degenerate case: plane passes through p2 perpendicular to v1
    return Kernel::Plane_3(p2, v1);
  }

  return Kernel::Plane_3(p2, bisector);
}

/**
 * @brief Project point onto plane along a direction (ray-plane intersection)
 */
auto
project_onto_bisector_plane(const Kernel::Point_3 &p_prev,
                            const Kernel::Point_3 &p,
                            const Kernel::Plane_3 &plane) -> Kernel::Point_3
{
  Kernel::Vector_3 v = p - p_prev;

  Kernel::FT numerator = -(plane.a() * p_prev.x() + plane.b() * p_prev.y() +
                           plane.c() * p_prev.z() + plane.d());
  Kernel::FT denominator =
      plane.a() * v.x() + plane.b() * v.y() + plane.c() * v.z();

  if (CGAL::abs(denominator) < Kernel::FT(1e-10)) {
    return p; // Parallel
  }

  Kernel::FT t = numerator / denominator;
  return p_prev + t * v;
}

/**
 * @brief Compute frame perpendicular to segment axis (for Segment Aligned)
 */
auto
compute_segment_frame(const Kernel::Vector_3 &axis,
                      const std::optional<Kernel::Vector_3> &ref_normal = std::nullopt) -> Frame
{
  Frame            frame;
  Kernel::Vector_3 normalized_axis = normalizeVector(axis);

  frame.tangent = normalized_axis;

  Kernel::Vector_3 perpendicular;

  if (ref_normal.has_value()) {
    // Use reference normal: project it onto plane perpendicular to tangent
    perpendicular = *ref_normal - (*ref_normal * frame.tangent) * frame.tangent;
  }

  // Fallback if no reference or degenerate reference
  if (!ref_normal.has_value() ||
      perpendicular * perpendicular < Kernel::FT(1e-10)) {
    // Try Z-axis first, if parallel try Y-axis (original method)
    perpendicular =
        CGAL::cross_product(normalized_axis, Kernel::Vector_3(0, 0, 1));
    if (perpendicular * perpendicular < Kernel::FT(1e-10)) {
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
  double           x      = CGAL::to_double(profile_pt.x()) - anchor_x;
  double           y      = CGAL::to_double(profile_pt.y()) - anchor_y;
  Kernel::Vector_3 offset = x * frame.normal + y * frame.binormal;
  return path_pt + offset;
}

/**
 * @brief Build sweep mesh by connecting consecutive profile instances
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
      // Two triangles instead of quad (avoids non-planar faces at corners)
      mesh.add_face(ring1[j], ring2[j], ring2[next_j]);
      mesh.add_face(ring1[j], ring2[next_j], ring1[next_j]);
    }
  }

  return vertex_rings;
}

/**
 * @brief Add flat end caps to mesh
 *
 * For closed profiles, adds the profile itself as a single polygonal face.
 * For open profiles, no caps are added (open-ended tube).
 *
 * @param add_start If true, adds the start cap
 * @param add_end If true, adds the end cap
 */
void
add_flat_caps(
    Surface_mesh_3                                               &mesh,
    [[maybe_unused]] const std::vector<Kernel::Point_3>          &path,
    const std::vector<std::vector<Surface_mesh_3::Vertex_index>> &vertex_rings,
    bool add_start, bool add_end)
{
  if (vertex_rings.empty()) {
    return;
  }

  const auto &first_ring = vertex_rings.front();
  const auto &last_ring  = vertex_rings.back();

  // Start cap: reversed orientation (inward-facing normal)
  if (add_start) {
    std::vector<Surface_mesh_3::Vertex_index> start_cap_vertices;
    start_cap_vertices.reserve(first_ring.size());
    for (auto it = first_ring.rbegin(); it != first_ring.rend(); ++it) {
      start_cap_vertices.push_back(*it);
    }

    auto start_face = mesh.add_face(start_cap_vertices);
    if (start_face == Surface_mesh_3::null_face()) {
      std::cerr << "Warning: Failed to add start cap face" << std::endl;
    }
  }

  // End cap: normal orientation (outward-facing normal)
  if (add_end) {
    auto end_face = mesh.add_face(last_ring);
    if (end_face == Surface_mesh_3::null_face()) {
      std::cerr << "Warning: Failed to add end cap face" << std::endl;
    }
  }
}

// ----------------------------------------------------------------------------
// Discrete Sweep Implementation (Segment Aligned)
// ----------------------------------------------------------------------------

/**
 * @brief Discrete sweep for Miter joins (Architecture/CAD style)
 *
 * This implementation treats each segment as a rigid tube and computes
 * intersections (miters) at corners. It avoids the "twisting" issues of
 * continuous frames at sharp corners.
 */
auto
sweep_discrete(const std::vector<Kernel::Point_3> &path_points,
               const Geometry &profile, const SweepOptions &options)
    -> std::unique_ptr<PolyhedralSurface>
{
  Surface_mesh_3                                         mesh;
  std::vector<std::vector<Surface_mesh_3::Vertex_index>> vertex_rings;
  std::vector<Kernel::Point_2> profile_points = extract_profile_points(profile);

  bool closed = options.closed_path;

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
    Kernel::Vector_3 axis = path_points[seg_idx + 1] - path_points[seg_idx];
    Frame            segment_frame = compute_segment_frame(axis, options.reference_normal);

    // --- Start Ring ---
    std::vector<Surface_mesh_3::Vertex_index> start_ring;

    if (seg_idx == 0) {
      // Create first ring from scratch
      std::vector<Kernel::Point_3> start_profile;
      start_profile.reserve(profile_points.size());

      for (const auto &profile_pt : profile_points) {
        start_profile.push_back(transform_profile_point(
            profile_pt, path_points[seg_idx], segment_frame, options.anchor_x,
            options.anchor_y));
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
      // Closing segment: Project start ring onto closing bisector
      std::vector<Kernel::Point_3> end_profile;
      end_profile.reserve(profile_points.size());

      for (const auto &profile_pt : profile_points) {
        end_profile.push_back(transform_profile_point(
            profile_pt, path_points[seg_idx + 1], segment_frame,
            options.anchor_x, options.anchor_y));
      }

      Kernel::Plane_3 closing_bisector = bisector_planes.back();

      for (size_t pt_idx = 0; pt_idx < first_start_ring.size(); ++pt_idx) {
        Kernel::Point_3 start_point  = mesh.point(start_ring[pt_idx]);
        Kernel::Point_3 target_point = end_profile[pt_idx];

        Kernel::Point_3 projected = project_onto_bisector_plane(
            start_point, target_point, closing_bisector);

        mesh.point(first_start_ring[pt_idx]) = projected;
      }
      end_ring = first_start_ring; // Loop back
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
          Kernel::Point_3 start_point = mesh.point(start_ring[pt_idx]);
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
    size_t n_pts = profile_points.size();
    for (size_t i = 0; i < n_pts - 1; ++i) {
      mesh.add_face(start_ring[i], start_ring[i + 1], end_ring[i + 1]);
      mesh.add_face(start_ring[i], end_ring[i + 1], end_ring[i]);
    }
    // Close tube
    if (profile.geometryTypeId() == TYPE_POLYGON ||
        (profile.geometryTypeId() == TYPE_LINESTRING &&
         profile.as<LineString>().isClosed())) {
      mesh.add_face(start_ring[n_pts - 1], start_ring[0], end_ring[0]);
      mesh.add_face(start_ring[n_pts - 1], end_ring[0], end_ring[n_pts - 1]);
    }

    if (seg_idx == 0) {
      vertex_rings.push_back(start_ring);
    }
    if (seg_idx == path_points.size() - 2) {
      vertex_rings.push_back(end_ring);
    }
  }

  if (!closed && (options.start_cap == SweepOptions::EndCapStyle::FLAT ||
                  options.end_cap == SweepOptions::EndCapStyle::FLAT)) {
    bool add_start = (options.start_cap == SweepOptions::EndCapStyle::FLAT);
    bool add_end   = (options.end_cap == SweepOptions::EndCapStyle::FLAT);
    add_flat_caps(mesh, path_points, vertex_rings, add_start, add_end);
  }

  return std::make_unique<PolyhedralSurface>(mesh);
}

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
  bool                         closed         = options.closed_path;

  // 1. Compute Frames
  std::vector<Frame> frames;
  if (options.frame_method == SweepOptions::FrameMethod::ROTATION_MINIMIZING) {
    frames = compute_rmf_frames(path_points, closed);
  } else if (options.frame_method == SweepOptions::FrameMethod::FRENET) {
    frames = compute_frenet_frames(path_points, closed);
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
    add_flat_caps(mesh, path_points, vertex_rings, add_start, add_end);
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

  // Force closed flag if geometry is closed
  SweepOptions opts = options;
  if (static_cast<bool>(isClosed(path))) {
    opts.closed_path = true;
  }

  // If closed is requested but path is not geometrically closed, close it by
  // duplicating start point
  if (opts.closed_path && !path_points.empty()) {
    if (path_points.front() != path_points.back()) {
      path_points.push_back(path_points.front());
    }
  }

  if (opts.frame_method == SweepOptions::FrameMethod::SEGMENT_ALIGNED) {
    return sweep_discrete(path_points, profile, opts);
  } else {
    return sweep_continuous(path_points, profile, opts);
  }
}

auto
create_circular_profile(double radius, int segments)
    -> std::unique_ptr<LineString>
{
  if (radius <= 0.0) {
    throw std::invalid_argument("Radius must be positive");
  }
  if (segments < 3) {
    throw std::invalid_argument("Segments must be at least 3");
  }

  std::vector<Point> points;
  points.reserve(segments);

  for (int i = 0; i < segments; ++i) {
    double angle = 2.0 * M_PI * i / segments;
    double x     = radius * std::cos(angle);
    double y     = radius * std::sin(angle);
    points.emplace_back(x, y, 0);
  }
  points.emplace_back(points[0]); // Close

  return std::make_unique<LineString>(points);
}

auto
create_rectangular_profile(double width, double height)
    -> std::unique_ptr<Polygon>
{
  if (width <= 0.0 || height <= 0.0) {
    throw std::invalid_argument("Dimensions must be positive");
  }

  double hw = width / 2.0;
  double hh = height / 2.0;

  std::vector<Point> points;
  points.reserve(5);
  points.emplace_back(-hw, -hh, 0);
  points.emplace_back(hw, -hh, 0);
  points.emplace_back(hw, hh, 0);
  points.emplace_back(-hw, hh, 0);
  points.emplace_back(-hw, -hh, 0);

  LineString ring(points);
  return std::make_unique<Polygon>(ring);
}

auto
create_chamfer_profile(double radius_x, double radius_y)
    -> std::unique_ptr<Polygon>
{
  if (radius_x <= 0.0) {
    throw std::invalid_argument("Radius X must be positive");
  }

  double ry = (radius_y < 0.0) ? radius_x : radius_y; // Default to symmetric
  if (ry <= 0.0) {
    throw std::invalid_argument("Radius Y must be positive");
  }

  // Create a triangle in the 3rd quadrant (negative X, negative Y)
  // Vertices: (0,0) -> (0, -ry) -> (-rx, 0) -> (0,0)
  std::vector<Point> points;
  points.reserve(4);
  points.emplace_back(0, 0, 0);
  points.emplace_back(0, -ry, 0);
  points.emplace_back(-radius_x, 0, 0);
  points.emplace_back(0, 0, 0);

  LineString ring(points);
  return std::make_unique<Polygon>(ring);
}

auto
create_fillet_profile(double radius, int segments) -> std::unique_ptr<Polygon>
{
  if (radius <= 0.0) {
    throw std::invalid_argument("Radius must be positive");
  }
  if (segments < 1) {
    throw std::invalid_argument("Segments must be at least 1");
  }

  // Center of the arc is at (-radius, -radius)
  double cx = -radius;
  double cy = -radius;

  std::vector<Point> points;
  points.reserve(segments + 3);

  points.emplace_back(0, 0, 0); // Origin

  // Arc from 0 degrees (0, -r) to 90 degrees (-r, 0) relative to center
  // In global coords, this goes from (0, -radius) to (-radius, 0)
  // We iterate backwards to keep consistent winding order if needed,
  // or just follow the path: (0,0) -> (0, -r) ... arc ... -> (-r, 0) -> (0,0)

  // Start of arc (0, -radius) corresponds to angle 0 relative to center (-r,
  // -r)? Center (-r, -r). Point (0, -r). x = -r + r*cos(0) = 0. y = -r +
  // r*sin(0) = -r. Correct.

  // End of arc (-r, 0).
  // x = -r + r*cos(90) = -r. y = -r + r*sin(90) = 0. Correct.

  // We want the points between angle 0 and 90 degrees.
  for (int i = 0; i <= segments; ++i) {
    double angle = (M_PI / 2.0) * (double(i) / segments);
    double x     = cx + radius * std::cos(angle);
    double y     = cy + radius * std::sin(angle);
    points.emplace_back(x, y, 0);
  }

  points.emplace_back(0, 0, 0); // Close back to origin

  LineString ring(points);
  return std::make_unique<Polygon>(ring);
}

} // namespace SFCGAL::algorithm
