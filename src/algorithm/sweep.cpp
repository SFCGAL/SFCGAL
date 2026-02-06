// Copyright (c) 2025-2026, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/algorithm/sweep.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/algorithm/isClosed.h>
#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/detail/tools/Registry.h>

#include <CGAL/Surface_mesh.h>

#include <cmath>
#include <stdexcept>

namespace SFCGAL::algorithm {

using Surface_mesh_3 = CGAL::Surface_mesh<Kernel::Point_3>;

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
 * @brief Check if path is closed using SFCGAL's isClosed algorithm
 */
auto
is_closed_path(const LineString &path) -> bool
{
  return static_cast<bool>(isClosed(path));
}

/**
 * @brief Compute initial frame from tangent direction
 */
auto
compute_initial_frame(const Kernel::Vector_3 &tangent) -> Frame
{
  Frame frame;
  frame.tangent = normalizeVector(tangent);

  // Choose helper vector perpendicular to tangent
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
 */
auto
propagate_frame_rmf(const Frame &prev_frame, const Kernel::Point_3 &prev_point,
                    const Kernel::Point_3  &curr_point,
                    const Kernel::Vector_3 &curr_tangent) -> Frame
{
  Frame frame;
  frame.tangent = normalizeVector(curr_tangent);

  // Vector from previous to current point
  Kernel::Vector_3 segment_vector    = curr_point - prev_point;
  Kernel::FT       segment_length_sq = segment_vector * segment_vector;

  if (segment_length_sq < Kernel::FT(1e-20)) {
    // Points are coincident, keep same frame
    frame.normal   = prev_frame.normal;
    frame.binormal = prev_frame.binormal;
    return frame;
  }

  // First reflection: reflect normal across plane perpendicular to
  // segment_vector
  Kernel::Vector_3 reflected_normal =
      prev_frame.normal - (Kernel::FT(2) / segment_length_sq) *
                              (segment_vector * prev_frame.normal) *
                              segment_vector;

  // Reflect previous tangent similarly
  Kernel::Vector_3 reflected_tangent =
      prev_frame.tangent - (Kernel::FT(2) / segment_length_sq) *
                               (segment_vector * prev_frame.tangent) *
                               segment_vector;

  // Second reflection to align with new tangent
  Kernel::Vector_3 tangent_diff    = frame.tangent - reflected_tangent;
  Kernel::FT       tangent_diff_sq = tangent_diff * tangent_diff;

  if (tangent_diff_sq < Kernel::FT(1e-20)) {
    // Tangents are parallel, use first reflection result
    frame.normal = reflected_normal;
  } else {
    frame.normal = reflected_normal - (Kernel::FT(2) / tangent_diff_sq) *
                                          (tangent_diff * reflected_normal) *
                                          tangent_diff;
  }

  frame.normal   = normalizeVector(frame.normal);
  frame.binormal = CGAL::cross_product(frame.tangent, frame.normal);

  return frame;
}

/**
 * @brief Correct holonomy for closed paths
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
 * @brief Compute RMF frames along path
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
 * @brief Compute frames using fixed up vector method
 *
 * This method maintains a constant "up" direction (typically Z-axis) for all
 * frames. This is ideal for planar paths where twist-free extrusion is desired.
 * The up vector is projected perpendicular to the tangent at each point.
 *
 * @param points Path points
 * @param up_vector Fixed up direction (default: Z-axis)
 * @return Vector of frames at each path point
 */
auto
compute_fixed_up_frames(const std::vector<Kernel::Point_3> &points,
                        const Kernel::Vector_3 &up_vector) -> std::vector<Frame>
{
  std::vector<Frame> frames;
  if (points.size() < 2) {
    return frames;
  }

  Kernel::Vector_3 up_normalized = normalizeVector(up_vector);

  for (size_t i = 0; i < points.size(); ++i) {
    Frame frame;

    // Compute tangent direction
    Kernel::Vector_3 tangent;
    if (i == 0) {
      // First point: use forward direction
      tangent = points[i + 1] - points[i];
    } else if (i == points.size() - 1) {
      // Last point: use backward direction
      tangent = points[i] - points[i - 1];
    } else {
      // Middle points: average of incoming and outgoing directions
      Kernel::Vector_3 incoming = points[i] - points[i - 1];
      Kernel::Vector_3 outgoing = points[i + 1] - points[i];
      tangent                   = incoming + outgoing;
    }

    frame.tangent = normalizeVector(tangent);

    // Project up_vector perpendicular to tangent to get binormal
    // binormal = up - (up · tangent) * tangent
    Kernel::FT       dot_product = up_normalized * frame.tangent;
    Kernel::Vector_3 projection  = dot_product * frame.tangent;
    Kernel::Vector_3 binormal    = up_normalized - projection;

    // Check if up_vector is too parallel to tangent
    Kernel::FT binormal_length_sq = binormal * binormal;
    if (binormal_length_sq < Kernel::FT(1e-10)) {
      // up_vector is nearly parallel to tangent, choose perpendicular vector
      Kernel::Vector_3 helper(1, 0, 0);
      if (std::abs(CGAL::to_double(frame.tangent * helper)) > 0.9) {
        helper = Kernel::Vector_3(0, 1, 0);
      }
      binormal = CGAL::cross_product(frame.tangent, helper);
    }

    frame.binormal = normalizeVector(binormal);
    frame.normal   = CGAL::cross_product(frame.tangent, frame.binormal);

    frames.push_back(frame);
  }

  return frames;
}

/**
 * @brief Compute bisector plane at a corner (miter join)
 *
 * Given three consecutive points p1, p2, p3 forming a corner at p2,
 * compute the plane perpendicular to the bisector of the angle.
 * This is used for creating sharp miter joins at corners.
 */
auto
compute_bisector_plane(const Kernel::Point_3 &p1, const Kernel::Point_3 &p2,
                       const Kernel::Point_3 &p3) -> Kernel::Plane_3
{
  Kernel::Vector_3 v1       = normalizeVector(p2 - p1);
  Kernel::Vector_3 v2       = normalizeVector(p3 - p2);
  Kernel::Vector_3 bisector = v1 + v2; // Bisector direction
  return Kernel::Plane_3(p2, bisector);
}

/**
 * @brief Project point onto plane along a direction
 *
 * Projects point p onto the given plane along the direction from p_prev to p.
 * This is used to create miter joins by intersecting swept profile segments
 * with bisector planes at corners.
 */
auto
project_onto_bisector_plane(const Kernel::Point_3 &p_prev,
                            const Kernel::Point_3 &p,
                            const Kernel::Plane_3 &plane) -> Kernel::Point_3
{
  Kernel::Vector_3 v = p - p_prev;

  // Compute intersection parameter t: plane.a*x + plane.b*y + plane.c*z +
  // plane.d = 0
  Kernel::FT numerator = -(plane.a() * p_prev.x() + plane.b() * p_prev.y() +
                           plane.c() * p_prev.z() + plane.d());
  Kernel::FT denominator =
      plane.a() * v.x() + plane.b() * v.y() + plane.c() * v.z();

  // Check if direction is parallel to plane
  if (CGAL::abs(denominator) < Kernel::FT(1e-10)) {
    return p; // Return original point if parallel
  }

  Kernel::FT t = numerator / denominator;
  return p_prev + t * v;
}

/**
 * @brief Compute frame perpendicular to segment axis
 *
 * Creates an orthonormal frame with tangent aligned to the segment axis.
 * This is used for SEGMENT_ALIGNED method where the frame is constant
 * along each segment (no twist).
 *
 * @param axis Normalized segment axis (direction from start to end)
 * @return Frame perpendicular to axis
 */
auto
compute_segment_frame(const Kernel::Vector_3 &axis) -> Frame
{
  Frame            frame;
  Kernel::Vector_3 normalized_axis = normalizeVector(axis);

  frame.tangent = normalized_axis;

  // Choose perpendicular vector by cross product with a reference vector
  // Try Z-axis first, if parallel try Y-axis
  Kernel::Vector_3 perpendicular =
      CGAL::cross_product(normalized_axis, Kernel::Vector_3(0, 0, 1));
  if (perpendicular * perpendicular < Kernel::FT(1e-10)) {
    // Axis is parallel to Z, use Y instead
    perpendicular =
        CGAL::cross_product(normalized_axis, Kernel::Vector_3(0, 1, 0));
  }
  frame.normal = normalizeVector(perpendicular);

  // Binormal is cross product of tangent and normal (right-handed system)
  frame.binormal =
      normalizeVector(CGAL::cross_product(frame.tangent, frame.normal));

  return frame;
}

/**
 * @brief Transform 2D profile point to 3D using frame
 * @param profile_pt 2D point in profile space
 * @param path_pt 3D point on the path
 * @param frame Orthonormal frame at path point
 * @param anchor_x X coordinate of anchor point in profile space
 * @param anchor_y Y coordinate of anchor point in profile space
 */
auto
transform_profile_point(const Kernel::Point_2 &profile_pt,
                        const Kernel::Point_3 &path_pt, const Frame &frame,
                        double anchor_x, double anchor_y) -> Kernel::Point_3
{
  // Profile coordinates: (X, Y) in 2D
  // Apply anchor point offset: subtract anchor from profile coordinates
  // This positions the anchor point on the path
  double x = CGAL::to_double(profile_pt.x()) - anchor_x;
  double y = CGAL::to_double(profile_pt.y()) - anchor_y;
  // Map X to normal direction, Y to binormal direction
  Kernel::Vector_3 offset = x * frame.normal + y * frame.binormal;

  return path_pt + offset;
}

/**
 * @brief Build sweep mesh by connecting consecutive profile instances
 * @return Vector of vertex rings (indices of vertices for each profile)
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

  // Add all vertices
  vertex_rings.reserve(n_path_points);

  for (const auto &profile : profiles) {
    std::vector<Surface_mesh_3::Vertex_index> ring;
    ring.reserve(profile.size());
    for (const auto &point : profile) {
      ring.push_back(mesh.add_vertex(point));
    }
    vertex_rings.push_back(std::move(ring));
  }

  // Connect consecutive rings with quads
  for (size_t i = 0; i < n_path_points - 1; ++i) {
    const auto &ring1 = vertex_rings[i];
    const auto &ring2 = vertex_rings[i + 1];

    for (size_t j = 0; j < n_profile_points; ++j) {
      size_t next_j =
          (j + 1) % n_profile_points; // Wrap around to close the tube

      // Create quad face with consistent orientation
      mesh.add_face(ring1[j], ring2[j], ring2[next_j], ring1[next_j]);
    }
  }

  return vertex_rings;
}

/**
 * @brief Add flat end caps to mesh
 */
void
add_flat_caps(
    Surface_mesh_3 &mesh, const std::vector<Kernel::Point_3> &path,
    const std::vector<Frame> & /*frames*/,
    const std::vector<std::vector<Surface_mesh_3::Vertex_index>> &vertex_rings)
{
  if (vertex_rings.empty()) {
    return;
  }

  // Get first and last vertex rings
  const auto &first_ring = vertex_rings.front();
  const auto &last_ring  = vertex_rings.back();

  size_t n_vertices = first_ring.size();

  // Add center vertices for triangular fan caps
  const Kernel::Point_3       &start_center     = path.front();
  Surface_mesh_3::Vertex_index start_center_idx = mesh.add_vertex(start_center);

  const Kernel::Point_3       &end_center     = path.back();
  Surface_mesh_3::Vertex_index end_center_idx = mesh.add_vertex(end_center);

  // Create start cap as triangular fan (inward-facing)
  for (size_t i = 0; i < n_vertices; ++i) {
    size_t next_i = (i + 1) % n_vertices;
    mesh.add_face(start_center_idx, first_ring[i], first_ring[next_i]);
  }

  // Create end cap as triangular fan (outward-facing)
  for (size_t i = 0; i < n_vertices; ++i) {
    size_t next_i = (i + 1) % n_vertices;
    mesh.add_face(end_center_idx, last_ring[next_i], last_ring[i]);
  }
}

} // anonymous namespace

// Public API implementation

auto
sweep(const LineString &path, const Geometry &profile,
      const SweepOptions &options) -> std::unique_ptr<PolyhedralSurface>
{
  // Validate inputs
  if (path.numPoints() < 2) {
    throw std::invalid_argument("Path must have at least 2 points, got: " +
                                std::to_string(path.numPoints()));
  }

  // Extract path points
  std::vector<Kernel::Point_3> path_points;
  path_points.reserve(path.numPoints());
  for (size_t i = 0; i < path.numPoints(); ++i) {
    const Point &point = path.pointN(i);
    path_points.emplace_back(point.x(), point.y(), point.z());
  }

  // Extract profile points
  std::vector<Kernel::Point_2> profile_points = extract_profile_points(profile);

  // Determine if path is closed
  bool closed = options.closed_path || is_closed_path(path);

  // For SEGMENT_ALIGNED method, use a completely different approach
  // This follows the old buffer3D.computeFlatBuffer() pattern:
  // - Process segment-by-segment (not point-by-point)
  // - Each segment has a constant frame perpendicular to its axis
  // - Profile is transformed at both ends of each segment with same frame
  // - Then vertices are projected onto bisector planes at corners for miter
  // joins
  // - Build mesh directly without intermediate storage to avoid multiple
  // transformations
  if (options.frame_method == SweepOptions::FrameMethod::SEGMENT_ALIGNED) {
    // Pre-compute bisector planes at intermediate points (corners)
    std::vector<Kernel::Plane_3> bisector_planes;
    if (path_points.size() >= 3) {
      bisector_planes.reserve(path_points.size() - 2);
      for (size_t i = 1; i < path_points.size() - 1; ++i) {
        bisector_planes.push_back(compute_bisector_plane(
            path_points[i - 1], path_points[i], path_points[i + 1]));
      }

      // For closed paths, add bisector plane at closing point
      // This is between last segment (pN-1, pN) and first segment (p0, p1)
      if (closed) {
        bisector_planes.push_back(compute_bisector_plane(
            path_points[path_points.size() - 2],
            path_points[path_points.size() - 1], path_points[1]));
      }
    }

    Surface_mesh_3                                         mesh;
    std::vector<std::vector<Surface_mesh_3::Vertex_index>> vertex_rings;

    // Previous segment's end ring - will be reused as next segment's start ring
    std::vector<Surface_mesh_3::Vertex_index> prev_end_ring;

    // For closed paths, save first segment's start ring to reuse as last
    // segment's end ring
    std::vector<Surface_mesh_3::Vertex_index> first_start_ring;

    // Process each segment independently
    for (size_t seg_idx = 0; seg_idx < path_points.size() - 1; ++seg_idx) {
      // 1. Compute axis for THIS segment (constant for entire segment)
      Kernel::Vector_3 axis = path_points[seg_idx + 1] - path_points[seg_idx];
      Frame            segment_frame = compute_segment_frame(axis);

      // 2. Handle start_ring: reuse previous end_ring or create new
      std::vector<Surface_mesh_3::Vertex_index> start_ring;

      if (seg_idx == 0) {
        // First segment: create start_ring from scratch
        std::vector<Kernel::Point_3> start_profile;
        start_profile.reserve(profile_points.size());

        for (const auto &profile_pt : profile_points) {
          start_profile.push_back(transform_profile_point(
              profile_pt, path_points[seg_idx], segment_frame, options.anchor_x,
              options.anchor_y));
        }

        // Add vertices to mesh (no projection needed at start of path)
        start_ring.reserve(profile_points.size());
        for (const auto &pt : start_profile) {
          start_ring.push_back(mesh.add_vertex(pt));
        }

        // Save first start ring for closed paths
        first_start_ring = start_ring;
      } else {
        // Subsequent segments: reuse previous segment's end_ring
        start_ring = prev_end_ring;
      }

      // 3. Create end_ring for this segment
      std::vector<Surface_mesh_3::Vertex_index> end_ring;

      // Note: For SEGMENT_ALIGNED, we don't reuse vertices for closed paths
      // because the projection logic requires independent vertex rings
      {
        // Create new end_ring
        std::vector<Kernel::Point_3> end_profile;
        end_profile.reserve(profile_points.size());

        for (const auto &profile_pt : profile_points) {
          end_profile.push_back(transform_profile_point(
              profile_pt, path_points[seg_idx + 1], segment_frame,
              options.anchor_x, options.anchor_y));
        }

        // 4. Project end profile onto bisector plane if needed and add to mesh
        end_ring.reserve(profile_points.size());

        for (size_t pt_idx = 0; pt_idx < profile_points.size(); ++pt_idx) {
          Kernel::Point_3 end_point = end_profile[pt_idx];

          // Project end point onto its bisector plane if it has one
          if (seg_idx < path_points.size() - 2) {
            // For projection, we need the line from start to end
            // Get start point from the mesh using start_ring index
            Kernel::Point_3 start_point = mesh.point(start_ring[pt_idx]);
            end_point = project_onto_bisector_plane(start_point, end_point,
                                                    bisector_planes[seg_idx]);
          }

          // Add vertex to mesh
          end_ring.push_back(mesh.add_vertex(end_point));
        }
      }

      // Save end_ring for next iteration
      prev_end_ring = end_ring;

      // 4. Connect faces between start and end profiles
      size_t n_pts = profile_points.size();
      for (size_t i = 0; i < n_pts - 1; ++i) {
        mesh.add_face(start_ring[i], start_ring[i + 1], end_ring[i + 1],
                      end_ring[i]);
      }

      // For closed profiles (Polygon or closed LineString), add wrap-around
      // face The profile extraction removes duplicate end point, so we need to
      // connect the last point back to the first to close the tube
      if (profile.geometryTypeId() == TYPE_POLYGON ||
          (profile.geometryTypeId() == TYPE_LINESTRING &&
           profile.as<LineString>().isClosed())) {
        mesh.add_face(start_ring[n_pts - 1], start_ring[0], end_ring[0],
                      end_ring[n_pts - 1]);
      }

      // Store rings for cap generation
      if (seg_idx == 0) {
        vertex_rings.push_back(start_ring);
      }
      if (seg_idx == path_points.size() - 2) {
        vertex_rings.push_back(end_ring);
      }
    }

    // Add end caps if requested and path is not closed
    if (!closed && (options.start_cap == SweepOptions::EndCapStyle::FLAT ||
                    options.end_cap == SweepOptions::EndCapStyle::FLAT)) {
      // For caps, we need frames at endpoints
      std::vector<Frame> endpoint_frames(2);

      // First segment frame for start cap
      Kernel::Vector_3 first_axis = path_points[1] - path_points[0];
      endpoint_frames[0]          = compute_segment_frame(first_axis);

      // Last segment frame for end cap
      Kernel::Vector_3 last_axis = path_points[path_points.size() - 1] -
                                   path_points[path_points.size() - 2];
      endpoint_frames[1] = compute_segment_frame(last_axis);

      add_flat_caps(mesh, path_points, endpoint_frames, vertex_rings);
    }

    return std::make_unique<PolyhedralSurface>(mesh);
  }

  // Compute frames along path for point-based methods
  std::vector<Frame> frames;
  if (options.frame_method == SweepOptions::FrameMethod::ROTATION_MINIMIZING) {
    frames = compute_rmf_frames(path_points, closed);
  } else if (options.frame_method == SweepOptions::FrameMethod::FIXED_UP) {
    frames = compute_fixed_up_frames(path_points, options.fixed_up_vector);
  } else {
    // TODO: Implement Frenet frame method
    throw std::runtime_error(
        "Only ROTATION_MINIMIZING, FIXED_UP and SEGMENT_ALIGNED frame methods "
        "are currently implemented");
  }

  // Transform profile at each path point using frames
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

  // Apply bisector plane projection for FIXED_UP method at corners
  // This creates miter joins for sharp corners (like AutoCAD)
  // Implementation follows the same logic as the old buffer3D computeFlatBuffer
  if (options.frame_method == SweepOptions::FrameMethod::FIXED_UP &&
      path_points.size() >= 3) {
    // Compute bisector planes at intermediate points (corners)
    std::vector<Kernel::Plane_3> bisector_planes;
    bisector_planes.reserve(path_points.size() - 2);

    for (size_t i = 1; i < path_points.size() - 1; ++i) {
      bisector_planes.push_back(compute_bisector_plane(
          path_points[i - 1], path_points[i], path_points[i + 1]));
    }

    // Process each segment and project vertices onto bisector planes
    // This is done segment-by-segment like the old buffer3D algorithm
    for (size_t i = 0; i < path_points.size() - 1; ++i) {
      for (size_t j = 0; j < profile_points.size(); ++j) {
        Kernel::Point_3 start_point = transformed_profiles[i][j];
        Kernel::Point_3 end_point   = transformed_profiles[i + 1][j];

        // Project start point onto its bisector plane if it has one
        if (i > 0) {
          start_point = project_onto_bisector_plane(start_point, end_point,
                                                    bisector_planes[i - 1]);
        }

        // Project end point onto its bisector plane if it has one
        if (i < path_points.size() - 2) {
          end_point = project_onto_bisector_plane(start_point, end_point,
                                                  bisector_planes[i]);
        }

        // Update the transformed profiles
        transformed_profiles[i][j]     = start_point;
        transformed_profiles[i + 1][j] = end_point;
      }
    }
  }

  // Build sweep mesh
  Surface_mesh_3 mesh;
  auto           vertex_rings = build_sweep_mesh(mesh, transformed_profiles);

  // Add end caps if requested
  if (!closed && (options.start_cap == SweepOptions::EndCapStyle::FLAT ||
                  options.end_cap == SweepOptions::EndCapStyle::FLAT)) {
    add_flat_caps(mesh, path_points, frames, vertex_rings);
  }

  return std::make_unique<PolyhedralSurface>(mesh);
}

auto
create_circular_profile(double radius, int segments)
    -> std::unique_ptr<LineString>
{
  if (radius <= 0.0) {
    throw std::invalid_argument("Radius must be positive, got: " +
                                std::to_string(radius));
  }
  if (segments < 3) {
    throw std::invalid_argument("Segments must be at least 3, got: " +
                                std::to_string(segments));
  }

  std::vector<Point> points;
  points.reserve(segments);

  for (int i = 0; i < segments; ++i) {
    double angle = 2.0 * M_PI * i / segments;
    double x     = radius * std::cos(angle);
    double y     = radius * std::sin(angle);
    points.emplace_back(x, y, 0);
  }
  // close the ring
  points.emplace_back(points[0]);

  return std::make_unique<LineString>(points);
}

auto
create_rectangular_profile(double width, double height)
    -> std::unique_ptr<Polygon>
{
  if (width <= 0.0) {
    throw std::invalid_argument("Width must be positive, got: " +
                                std::to_string(width));
  }
  if (height <= 0.0) {
    throw std::invalid_argument("Height must be positive, got: " +
                                std::to_string(height));
  }

  double half_width  = width / 2.0;
  double half_height = height / 2.0;

  std::vector<Point> points;
  points.reserve(5); // 4 corners + closure

  points.emplace_back(-half_width, -half_height, 0);
  points.emplace_back(half_width, -half_height, 0);
  points.emplace_back(half_width, half_height, 0);
  points.emplace_back(-half_width, half_height, 0);
  points.emplace_back(-half_width, -half_height, 0); // close the ring

  LineString ring(points);
  return std::make_unique<Polygon>(ring);
}

} // namespace SFCGAL::algorithm
