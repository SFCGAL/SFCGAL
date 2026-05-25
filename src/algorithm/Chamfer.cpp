// Copyright (c) 2025-2026, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/Chamfer.h"
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
#include <CGAL/squared_distance_3.h>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace SFCGAL::algorithm {

using Kernel         = SFCGAL::Kernel;
using Point_3        = Kernel::Point_3;
using Vector_3       = Kernel::Vector_3;
using Surface_mesh_3 = SFCGAL::Surface_mesh_3;

namespace PMP = CGAL::Polygon_mesh_processing;

namespace {

// Helper to normalize vectors
auto
normalize(const Vector_3 &v) -> Vector_3
{
  const double len = std::sqrt(CGAL::to_double(v.squared_length()));
  if (len < 1e-12) {
    return Vector_3(0, 0, 0);
  }
  return v / len;
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

// Helper to create profile in a specific quadrant
auto
create_chamfer_in_quadrant(double r, double ry, int quadrant)
    -> std::unique_ptr<Polygon>
{
  double radius_y = (ry < 0) ? r : ry;

  std::vector<Point> points;
  points.reserve(4);

  switch (quadrant) {
  case 1: // +N, +B
    points.emplace_back(0, 0, 0);
    points.emplace_back(r, 0, 0);
    points.emplace_back(0, radius_y, 0);
    points.emplace_back(0, 0, 0);
    break;
  case 2: // -N, +B
    points.emplace_back(0, 0, 0);
    points.emplace_back(-r, 0, 0);
    points.emplace_back(0, radius_y, 0);
    points.emplace_back(0, 0, 0);
    break;
  case 3: // -N, -B
    points.emplace_back(0, 0, 0);
    points.emplace_back(-r, 0, 0);
    points.emplace_back(0, -radius_y, 0);
    points.emplace_back(0, 0, 0);
    break;
  case 4: // +N, -B
    points.emplace_back(0, 0, 0);
    points.emplace_back(r, 0, 0);
    points.emplace_back(0, -radius_y, 0);
    points.emplace_back(0, 0, 0);
    break;
  default:
    throw std::invalid_argument("Invalid quadrant: must be 1, 2, 3, or 4");
  }

  return std::make_unique<Polygon>(LineString(points));
}

// Helper to create fillet profile (quarter circle) in a specific quadrant
auto
create_fillet_in_quadrant(double radius, int segments, int quadrant)
    -> std::unique_ptr<Polygon>
{
  if (radius <= 0.0) {
    throw std::invalid_argument("Radius must be positive");
  }
  if (segments < 1) {
    throw std::invalid_argument("Segments must be at least 1");
  }

  std::vector<Point> points;
  points.reserve(segments + 3);

  points.emplace_back(0, 0, 0); // Origin

  // Center must be offset in the quadrant direction for concave arc
  // Pattern: center at (±radius, ±radius) depending on quadrant
  double cx, cy, start_angle, end_angle;

  switch (quadrant) {
  case 1: // +N, +B: arc from (r, 0) to (0, r)
    cx          = radius;
    cy          = radius;
    start_angle = 3.0 * M_PI / 2.0; // 270° → (r, 0)
    end_angle   = M_PI;             // 180° → (0, r)
    break;
  case 2: // -N, +B: arc from (-r, 0) to (0, r)
    cx          = -radius;
    cy          = radius;
    start_angle = 3.0 * M_PI / 2.0; // 270° → (-r, 0)
    end_angle   = 2.0 * M_PI;       // 360° → (0, r)
    break;
  case 3: // -N, -B: arc from (0, -r) to (-r, 0)
    cx          = -radius;
    cy          = -radius;
    start_angle = 0.0;        // 0° → (0, -r)
    end_angle   = M_PI / 2.0; // 90° → (-r, 0)
    break;
  case 4: // +N, -B: arc from (r, 0) to (0, -r)
    cx          = radius;
    cy          = -radius;
    start_angle = M_PI / 2.0; // 90° → (r, 0)
    end_angle   = M_PI;       // 180° → (0, -r)
    break;
  default:
    throw std::invalid_argument("Invalid quadrant: must be 1, 2, 3, or 4");
  }

  // Generate arc points
  for (int i = 0; i <= segments; ++i) {
    double t     = double(i) / segments;
    double angle = start_angle + t * (end_angle - start_angle);
    double x     = cx + radius * std::cos(angle);
    double y     = cy + radius * std::sin(angle);
    points.emplace_back(x, y, 0);
  }

  points.emplace_back(0, 0, 0); // Close back to origin

  return std::make_unique<Polygon>(LineString(points));
}

// Determine quadrant based on edge orientation and face normals
auto
determine_quadrant(const Vector_3 &edge_dir, const Vector_3 &n1,
                   const Vector_3 &n2) -> int
{
  // Bisectrice pointant vers l'intérieur du matériau
  const Vector_3 bisector_inward = normalize(-(n1 + n2));

  // Calculer le repère (T, N, B) exactement comme compute_segment_frame()
  const Vector_3 T = normalize(edge_dir);

  // Essayer Z-axis d'abord, si parallèle essayer Y-axis
  Vector_3 perpendicular = CGAL::cross_product(T, Vector_3(0, 0, 1));
  if (CGAL::to_double(perpendicular.squared_length()) < 1e-10) {
    perpendicular = CGAL::cross_product(T, Vector_3(0, 1, 0));
  }
  const Vector_3 N = normalize(perpendicular);
  const Vector_3 B = normalize(CGAL::cross_product(T, N));

  // Projeter la bisectrice sur le plan (N, B)
  const double proj_N = CGAL::to_double(bisector_inward * N);
  const double proj_B = CGAL::to_double(bisector_inward * B);

  // Determine quadrant
  if (proj_N > 0 && proj_B > 0)
    return 1;
  if (proj_N < 0 && proj_B > 0)
    return 2;
  if (proj_N < 0 && proj_B < 0)
    return 3;
  return 4; // proj_N > 0 && proj_B < 0
}

// Create cutter solid for one edge
auto
create_cutter_for_edge(const Geometry &solid_geom, const LineString &edge,
                       const ChamferOptions &options)
    -> std::unique_ptr<Geometry>
{
  if (edge.numPoints() < 2) {
    throw std::invalid_argument("Edge must have at least 2 points");
  }

  // Convert SOLID to Surface_mesh_3 for analysis
  Surface_mesh_3 mesh;
  if (solid_geom.geometryTypeId() == TYPE_SOLID) {
    mesh = solid_geom.as<Solid>().toSurfaceMesh();
  } else if (solid_geom.geometryTypeId() == TYPE_POLYHEDRALSURFACE) {
    mesh = solid_geom.as<PolyhedralSurface>().toSurfaceMesh();
  } else {
    throw std::invalid_argument(
        "Input geometry must be a Solid or PolyhedralSurface");
  }

  // Get first segment of edge for orientation
  const Point_3 p1 = edge.pointN(0).toPoint_3();
  const Point_3 p2 = edge.pointN(1).toPoint_3();

  // Find corresponding halfedge in mesh
  const auto hd = find_halfedge(mesh, p1, p2, options.epsilon);
  if (hd == Surface_mesh_3::null_halfedge()) {
    throw std::runtime_error("Edge not found in solid mesh. "
                             "Ensure the edge coincides with a mesh edge.");
  }

  // Get the two incident faces
  const auto f1 = mesh.face(hd);
  const auto f2 = mesh.face(mesh.opposite(hd));

  if (f1 == Surface_mesh_3::null_face() || f2 == Surface_mesh_3::null_face()) {
    throw std::runtime_error(
        "Edge is not shared by two faces (boundary edge?)");
  }

  // Compute face normals using CGAL PMP
  const Vector_3 n1 = PMP::compute_face_normal(f1, mesh);
  const Vector_3 n2 = PMP::compute_face_normal(f2, mesh);

  // Validate that faces are approximately perpendicular (90 degrees)
  const double dot                 = CGAL::to_double(n1 * n2);
  const double angle_between_faces = std::acos(std::clamp(dot, -1.0, 1.0));
  const double angle_diff          = std::abs(angle_between_faces - M_PI / 2.0);

  constexpr double TOLERANCE = 0.1; // ~5.7 degrees
  if (angle_diff > TOLERANCE) {
    throw std::invalid_argument(
        "Chamfer requires perpendicular faces (90°) in this version. "
        "Face angle: " +
        std::to_string(angle_between_faces * 180.0 / M_PI) + "°");
  }

  // Determine edge direction and quadrant
  const Vector_3 edge_dir = normalize(p2 - p1);
  const int      quadrant = determine_quadrant(edge_dir, n1, n2);

  // Create profile directly in the correct quadrant
  std::unique_ptr<Polygon> profile;
  if (options.type == ChamferType::ROUND) {
    // Create fillet profile in the correct quadrant
    profile =
        create_fillet_in_quadrant(options.radius, options.segments, quadrant);
  } else {
    profile =
        create_chamfer_in_quadrant(options.radius, options.radius_y, quadrant);
  }

  // Sweep profile along edge to create cutter
  SweepOptions sweep_opts;
  sweep_opts.frame_method = SweepOptions::FrameMethod::SEGMENT_ALIGNED;
  sweep_opts.closed_path  = isClosed(edge);

  auto cutter_surf  = sweep(edge, *profile, sweep_opts);
  auto cutter_solid = std::make_unique<Solid>(cutter_surf.release());

  // Return the cutter (difference will be applied later)
  return cutter_solid;
}

} // anonymous namespace

auto
chamfer(const Geometry &solid_geom, const Geometry &edge_geom,
        const ChamferOptions &options) -> std::unique_ptr<Geometry>
{
  // Collect edges from input geometry
  std::vector<const LineString *> edges;

  if (edge_geom.geometryTypeId() == TYPE_LINESTRING) {
    edges.push_back(&edge_geom.as<LineString>());
  } else if (edge_geom.geometryTypeId() == TYPE_MULTILINESTRING) {
    const auto &multi = edge_geom.as<MultiLineString>();
    for (size_t i = 0; i < multi.numGeometries(); ++i) {
      edges.push_back(&multi.lineStringN(i));
    }
  } else {
    throw std::invalid_argument(
        "Edge geometry must be LineString or MultiLineString");
  }

  // Create all cutters first (before modifying the solid)
  std::vector<std::unique_ptr<Geometry>> cutters;

  for (const LineString *edge : edges) {
    try {
      auto cutter = create_cutter_for_edge(solid_geom, *edge, options);
      cutters.push_back(std::move(cutter));
    } catch (const std::exception &e) {
      // Skip edges that cannot be chamfered
      std::cerr << "Warning: skipping edge - " << e.what() << std::endl;
    }
  }

  if (cutters.empty()) {
    // No valid cutters, return original geometry
    return solid_geom.clone();
  }

  // Combine all cutters using union
  std::unique_ptr<Geometry> combined_cutter = std::move(cutters[0]);
  for (size_t i = 1; i < cutters.size(); ++i) {
    combined_cutter = union3D(*combined_cutter, *cutters[i]);
  }

  // Apply single difference at the end
  return difference3D(solid_geom, *combined_cutter);
}

} // namespace SFCGAL::algorithm
