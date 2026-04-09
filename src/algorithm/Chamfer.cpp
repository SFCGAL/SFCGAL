// Copyright (c) 2025-2026, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/Chamfer.h"
#include "SFCGAL/algorithm/collectionHomogenize.h"
#include "SFCGAL/algorithm/difference.h"
#include "SFCGAL/algorithm/isClosed.h"
#include "SFCGAL/algorithm/sweep.h"
#include "SFCGAL/algorithm/union.h"

#include <SFCGAL/Geometry.h>
#include <SFCGAL/Kernel.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/squared_distance_3.h>

#include <SFCGAL/detail/tools/Log.h>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace SFCGAL::algorithm {

using Point_3        = Kernel::Point_3;
using Vector_3       = Kernel::Vector_3;
using Surface_mesh_3 = SFCGAL::Surface_mesh_3;

namespace PMP = CGAL::Polygon_mesh_processing;

namespace {

// Thresholds for numerical comparisons
constexpr double ZERO_LENGTH_SQ =
    1e-20; // Squared length below which a vector is zero
constexpr double NEAR_ZERO_LENGTH =
    1e-10; // Length below which a projection is degenerate
constexpr double EPS_SCALE = 1e-3; // Profile origin shift as fraction of radius
constexpr double MIN_OPENING_DEG =
    10.0; // Minimum supported opening angle (degrees)
constexpr double MAX_OPENING_DEG =
    170.0; // Maximum supported opening angle (degrees)

// ============================================================================
// Algorithm overview:
// 1. Convert solid to CGAL Surface_mesh (once for all edges)
// 2. Decompose multi-segment LineStrings into individual 2-point segments
// 3. For each segment: find halfedge, check convexity, compute dihedral angle,
//    build angle-adapted cutting profile, sweep along edge → cutter solid
// 4. Union all cutters
// 5. Boolean difference: solid - combined_cutter
// Edges that fail step 3 are skipped with a warning (not propagated).
// ============================================================================

// Normalize a vector to unit length.
// Returns zero vector if input length is below threshold.
// Uses CGAL::to_double for sqrt — unavoidable with EPECK kernel.
auto
normalize(const Vector_3 &v) -> Vector_3
{
  const Kernel::FT len_sq = v * v;
  if (len_sq < Kernel::FT(ZERO_LENGTH_SQ)) {
    return Vector_3(0, 0, 0);
  }
  return v / std::sqrt(CGAL::to_double(len_sq));
}

// Find the halfedge in mesh that matches edge (p1 -> p2)
auto
find_halfedge(const Surface_mesh_3 &mesh, const Point_3 &p1, const Point_3 &p2,
              double epsilon = 1e-6) -> Surface_mesh_3::Halfedge_index
{
  const double eps_sq = epsilon * epsilon;

  for (auto h : mesh.halfedges()) {
    const Point_3 &src = mesh.point(mesh.source(h));
    const Point_3 &tgt = mesh.point(mesh.target(h));

    const double d1 = CGAL::to_double(CGAL::squared_distance(src, p1));
    const double d2 = CGAL::to_double(CGAL::squared_distance(tgt, p2));

    if (d1 < eps_sq && d2 < eps_sq) {
      return h;
    }
  }

  return Surface_mesh_3::null_halfedge();
}

// Compute face surface directions (tangent to face, pointing toward solid
// interior) in the local (N, B) frame.
//
// In this frame, n1_perp is at angle 0 and n2_perp at angle theta_2.
// The tangent to each face is perpendicular to its normal. The correct
// tangent (pointing inward) is chosen by comparing with the inward bisector.
//
// Returns {f1_angle, f2_angle} in radians.
auto
compute_face_surface_angles(double theta_2) -> std::pair<double, double>
{
  const double sgn = (theta_2 > 0) ? 1.0 : -1.0;

  // Face 1 surface direction: ⊥ to n1_perp (angle 0), pointing inward
  const double f1_angle = -sgn * M_PI / 2.0;

  // Face 2 surface direction: ⊥ to n2_perp (angle θ₂), pointing inward
  const double f2_angle = theta_2 + sgn * M_PI / 2.0;

  return {f1_angle, f2_angle};
}

// Create chamfer triangle profile parameterized by the angle of n2
// in the local frame (N, B).
//
// Legs follow the face surfaces (tangent to faces, not along normals).
// This is correct for any dihedral angle, not just 90°.
auto
create_chamfer_profile_for_angle(double r1, double r2, double theta_2)
    -> std::unique_ptr<Polygon>
{
  auto [f1_angle, f2_angle] = compute_face_surface_angles(theta_2);

  // Leg endpoints at distance r along face surface directions
  const double leg1_x = r1 * std::cos(f1_angle);
  const double leg1_y = r1 * std::sin(f1_angle);

  const double leg2_x = r2 * std::cos(f2_angle);
  const double leg2_y = r2 * std::sin(f2_angle);

  // Small outward extension along the outward bisector of the two normals
  // to avoid coplanar faces in boolean ops.
  const double eps_shift =
      std::max(std::min(r1, r2) * EPS_SCALE, NEAR_ZERO_LENGTH);
  const double bis_angle = theta_2 / 2.0;
  const double origin_x  = eps_shift * std::cos(bis_angle);
  const double origin_y  = eps_shift * std::sin(bis_angle);

  std::vector<Point> points;
  points.reserve(4);
  points.emplace_back(origin_x, origin_y, 0);
  points.emplace_back(leg1_x, leg1_y, 0);
  points.emplace_back(leg2_x, leg2_y, 0);
  points.emplace_back(origin_x, origin_y, 0);

  return std::make_unique<Polygon>(LineString(points));
}

// Create fillet arc profile parameterized by the angle of n2.
//
// The fillet arc of radius R is tangent to both face surfaces.
// The tangent points (legs) are at distance R/tan(γ/2) from the edge,
// where γ is the opening angle. The arc center lies at the intersection
// of lines parallel to each face, offset inward by R.
auto
create_fillet_profile_for_angle(double radius, int segments, double theta_2)
    -> std::unique_ptr<Polygon>
{
  if (radius <= 0.0) {
    throw std::invalid_argument("Radius must be positive");
  }
  if (segments < 1) {
    throw std::invalid_argument("Segments must be at least 1");
  }

  auto [f1_angle, f2_angle] = compute_face_surface_angles(theta_2);

  // Opening angle γ = π - |θ₂|
  const double gamma       = M_PI - std::abs(theta_2);
  const double half_gamma  = gamma / 2.0;
  const double sin_theta_2 = std::sin(theta_2);
  const double tan_half    = std::tan(half_gamma);

  // Guard against degenerate angles (coplanar or knife-edge faces)
  if (std::abs(sin_theta_2) < NEAR_ZERO_LENGTH ||
      std::abs(tan_half) < NEAR_ZERO_LENGTH) {
    throw std::invalid_argument(
        "Degenerate dihedral angle for fillet: faces are nearly "
        "coplanar or form a knife-edge");
  }

  // Leg distance from edge = R / tan(γ/2)
  const double leg_dist = radius / tan_half;

  const double leg1_x = leg_dist * std::cos(f1_angle);
  const double leg1_y = leg_dist * std::sin(f1_angle);

  const double leg2_x = leg_dist * std::cos(f2_angle);
  const double leg2_y = leg_dist * std::sin(f2_angle);

  // Arc center: intersection of face lines offset inward by R.
  // In (N, B) frame:
  //   Face 1 equation: x = 0 → offset: x = -R
  //   Face 2 equation: x·cos(θ₂) + y·sin(θ₂) = 0 → offset: = -R
  // Solving: cx = -R, cy = R·(cos(θ₂) - 1) / sin(θ₂)
  const double cx = -radius;
  const double cy = radius * (std::cos(theta_2) - 1.0) / sin_theta_2;

  // Arc from leg1 to leg2 around center (cx, cy)
  double start_angle = std::atan2(leg1_y - cy, leg1_x - cx);
  double end_angle   = std::atan2(leg2_y - cy, leg2_x - cx);

  // Take the shorter arc through the exterior
  double delta = end_angle - start_angle;
  while (delta > M_PI) {
    delta -= 2.0 * M_PI;
  }
  while (delta < -M_PI) {
    delta += 2.0 * M_PI;
  }

  // Origin: small outward shift along the outward bisector
  const double eps_shift = std::max(radius * EPS_SCALE, NEAR_ZERO_LENGTH);
  const double bis_angle = theta_2 / 2.0;
  const double origin_x  = eps_shift * std::cos(bis_angle);
  const double origin_y  = eps_shift * std::sin(bis_angle);

  std::vector<Point> points;
  points.reserve(segments + 3);

  points.emplace_back(origin_x, origin_y, 0);

  for (int i = 0; i <= segments; ++i) {
    const double t     = double(i) / segments;
    const double angle = start_angle + t * delta;
    const double x     = cx + radius * std::cos(angle);
    const double y     = cy + radius * std::sin(angle);
    points.emplace_back(x, y, 0);
  }

  points.emplace_back(origin_x, origin_y, 0);

  return std::make_unique<Polygon>(LineString(points));
}

// Compute the signed angle of n2 relative to n1 in the plane perpendicular
// to the edge. Looking along the edge: n1's projection defines angle 0,
// and this returns the angle to n2's projection (counter-clockwise positive).
//   |theta_2| ≈ pi/2  → 90° corner (cube)
//   |theta_2| ≈ pi     → nearly flat (coplanar faces)
//   |theta_2| ≈ 0      → knife-edge
auto
compute_n2_angle(const Vector_3 &edge_dir, const Vector_3 &n1,
                 const Vector_3 &n2) -> double
{
  const Vector_3 T = normalize(edge_dir);

  // N = projection of n1 perpendicular to edge
  Vector_3 N     = n1 - (n1 * T) * T;
  double   N_len = std::sqrt(CGAL::to_double(N.squared_length()));

  if (N_len < NEAR_ZERO_LENGTH) {
    throw std::invalid_argument("Face normal n1 is parallel to edge direction");
  }
  N = N / N_len;

  // B completes the right-handed frame
  const Vector_3 B = normalize(CGAL::cross_product(T, N));

  // Project n2 perpendicular to edge
  Vector_3 n2_perp     = n2 - (n2 * T) * T;
  double   n2_perp_len = std::sqrt(CGAL::to_double(n2_perp.squared_length()));

  if (n2_perp_len < NEAR_ZERO_LENGTH) {
    throw std::invalid_argument("Face normal n2 is parallel to edge direction");
  }
  n2_perp = n2_perp / n2_perp_len;

  // Angle of n2_perp in (N, B) plane
  const double proj_N = CGAL::to_double(n2_perp * N);
  const double proj_B = CGAL::to_double(n2_perp * B);

  return std::atan2(proj_B, proj_N);
}

// Create the cutter solid for a single 2-point edge segment.
// Finds the halfedge, checks convexity, computes dihedral angle,
// builds the cutting profile, and sweeps it along the edge.
// Throws std::invalid_argument on validation failure.
auto
create_cutter_for_edge(const Surface_mesh_3 &mesh, const LineString &edge,
                       const ChamferOptions &options)
    -> std::unique_ptr<Geometry>
{
  if (edge.numPoints() < 2) {
    throw std::invalid_argument("Edge must have at least 2 points");
  }

  // Get first segment of edge for orientation
  const Point_3 p1 = edge.pointN(0).toPoint_3();
  const Point_3 p2 = edge.pointN(1).toPoint_3();

  // Find corresponding halfedge in mesh
  const auto hd = find_halfedge(mesh, p1, p2, options.epsilon);
  if (hd == Surface_mesh_3::null_halfedge()) {
    throw std::invalid_argument("Edge not found in solid mesh. "
                                "Ensure the edge coincides with a mesh edge.");
  }

  // Get the two incident faces
  const auto f1 = mesh.face(hd);
  const auto f2 = mesh.face(mesh.opposite(hd));

  if (f1 == Surface_mesh_3::null_face() || f2 == Surface_mesh_3::null_face()) {
    throw std::invalid_argument(
        "Edge is not shared by two faces (boundary edge?)");
  }

  // Compute face normals using CGAL PMP
  const Vector_3 n1 = PMP::compute_face_normal(f1, mesh);
  const Vector_3 n2 = PMP::compute_face_normal(f2, mesh);

  // Reject concave (reflex) edges: the opposite vertex of f1 must lie
  // on the interior side of f2's plane (dot with n2 < 0).
  {
    const auto     v_opp  = mesh.target(mesh.next(hd));
    const Vector_3 to_opp = mesh.point(v_opp) - p1;
    if (CGAL::to_double(to_opp * n2) >= 0.0) {
      throw std::invalid_argument("Edge is concave (reflex). "
                                  "Chamfer is only supported on convex edges.");
    }
  }

  // Compute dihedral angle between faces (angle between outward normals)
  const double dot_val       = CGAL::to_double(n1 * n2);
  const double alpha         = std::acos(std::clamp(dot_val, -1.0, 1.0));
  const double opening_angle = M_PI - alpha;

  // Validate opening angle range
  constexpr double MIN_OPENING = MIN_OPENING_DEG * M_PI / 180.0;
  constexpr double MAX_OPENING = MAX_OPENING_DEG * M_PI / 180.0;

  if (opening_angle < MIN_OPENING) {
    throw std::invalid_argument("Edge opening angle too small (" +
                                std::to_string(opening_angle * 180.0 / M_PI) +
                                "°). Minimum supported: 10°.");
  }
  if (opening_angle > MAX_OPENING) {
    throw std::invalid_argument("Edge opening angle too large (" +
                                std::to_string(opening_angle * 180.0 / M_PI) +
                                "°). Maximum supported: 170°.");
  }

  // Compute continuous angle of n2 in the local frame
  const Vector_3 edge_dir = normalize(p2 - p1);
  const double   theta_2  = compute_n2_angle(edge_dir, n1, n2);

  // Create profile based on actual dihedral angle
  std::unique_ptr<Polygon> profile;
  if (options.type == ChamferType::ROUND) {
    profile = create_fillet_profile_for_angle(options.radius, options.segments,
                                              theta_2);
  } else {
    const double ry =
        (options.radius_y < 0) ? options.radius : options.radius_y;
    profile = create_chamfer_profile_for_angle(options.radius, ry, theta_2);
  }

  // Sweep profile along edge to create cutter
  SweepOptions sweep_opts;
  sweep_opts.frame_method = SweepOptions::FrameMethod::SEGMENT_ALIGNED;
  sweep_opts.closed_path  = isClosed(edge);

  // Do NOT set reference_normal: the fallback in compute_segment_frame uses
  // tangent × Z_up which gives binormal=(0,0,-1) consistently for all
  // horizontal tangent directions (+X, -X, +Y, -Y), preventing frame flips.

  std::cout << "ici\n";
  auto cutter_surf = sweep(edge, *profile, sweep_opts);
  std::cout << cutter_surf->asText(8) << "\n";

  // Repair self-intersections that can arise at concave corners of the path
  // (miter joins extend the tube through the corner, creating face overlaps).
  Surface_mesh_3 repair_mesh = cutter_surf->toSurfaceMesh();
  if (!PMP::experimental::autorefine_and_remove_self_intersections(repair_mesh)) {
    SFCGAL_WARNING("Chamfer: could not fully repair self-intersections in cutter");
  }
  cutter_surf = std::make_unique<PolyhedralSurface>(repair_mesh);

  // if (!isClosed(*cutter_surf)) {
  //   throw std::invalid_argument(
  //       "Sweep produced an open surface (non-closed cutter)");
  // }

  return std::make_unique<Solid>(std::move(cutter_surf));
}

} // anonymous namespace

auto
chamfer(const Geometry &solid_geom, const Geometry &edge_geom,
        const ChamferOptions &options) -> std::unique_ptr<Geometry>
{
  // Validate inputs
  if (options.radius <= 0.0) {
    throw std::invalid_argument("Chamfer radius must be positive");
  }
  if (options.type == ChamferType::ROUND && options.segments < 1) {
    throw std::invalid_argument("Fillet segments must be at least 1");
  }

  // Convert solid to Surface_mesh_3 once for all edges
  Surface_mesh_3 mesh;
  if (solid_geom.geometryTypeId() == TYPE_SOLID) {
    mesh = solid_geom.as<Solid>().toSurfaceMesh();
  } else if (solid_geom.geometryTypeId() == TYPE_POLYHEDRALSURFACE) {
    mesh = solid_geom.as<PolyhedralSurface>().toSurfaceMesh();
  } else {
    throw std::invalid_argument(
        "Input geometry must be a Solid or PolyhedralSurface");
  }

  // Create all cutters first (before modifying the solid)
  std::vector<std::unique_ptr<Geometry>> cutters;

  if (edge_geom.geometryTypeId() == TYPE_LINESTRING) {
    try {
      auto cutter =
          create_cutter_for_edge(mesh, edge_geom.as<LineString>(), options);
      std::cout << cutter->asText(4) << "\n";
      cutters.push_back(std::move(cutter));
    } catch (const std::invalid_argument &e) {
      // Expected: edge not found, concave, angle out of range
      SFCGAL_WARNING(std::string("Chamfer: skipping edge - ") + e.what());
    }
  } else if (edge_geom.geometryTypeId() == TYPE_MULTILINESTRING) {
    const auto &multi = edge_geom.as<MultiLineString>();
    for (size_t i = 0; i < multi.numGeometries(); ++i) {
      try {
        auto cutter = create_cutter_for_edge(
            mesh, multi.geometryN(i).as<LineString>(), options);
        cutters.push_back(std::move(cutter));
      } catch (const std::invalid_argument &e) {
        // Expected: edge not found, concave, angle out of range
        SFCGAL_WARNING(std::string("Chamfer: skipping edge - ") + e.what());
      }
    }
  } else {
    throw std::invalid_argument(
        "Edge geometry must be LineString or MultiLineString");
  }

  if (cutters.empty()) {
    return solid_geom.clone();
  }

  // Combine all cutters using union
  std::unique_ptr<Geometry> combined_cutter = std::move(cutters[0]);
  for (size_t i = 1; i < cutters.size(); ++i) {
    combined_cutter = union3D(*combined_cutter, *cutters[i]);
  }

  // Apply single difference at the end
  auto result = difference3D(solid_geom, *combined_cutter);
  return collectionHomogenize(std::move(result));
}

} // namespace SFCGAL::algorithm
