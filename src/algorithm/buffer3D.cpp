// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/buffer3D.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/algorithm/minkowskiSum3D.h"
#include "SFCGAL/algorithm/sweep.h"
#include "SFCGAL/algorithm/union.h"
#include "SFCGAL/detail/GeometrySet.h"
#include "SFCGAL/numeric.h"
#include "SFCGAL/primitive3d/Cylinder.h"
#include "SFCGAL/primitive3d/Sphere.h"
#include "SFCGAL/triangulate/triangulatePolygon.h"
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/minkowski_sum_3.h>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace SFCGAL::algorithm {

Buffer3D::Buffer3D(const Geometry &inputGeometry, double radius, int segments)
    : _radius(radius), _segments(segments)
{
  if (inputGeometry.is<Point>()) {
    _inputPoints.push_back(inputGeometry.as<Point>().toPoint_3());
  } else if (inputGeometry.is<LineString>()) {
    const auto &ls = inputGeometry.as<LineString>();
    for (size_t i = 0; i < ls.numPoints(); ++i) {
      _inputPoints.push_back(ls.pointN(i).toPoint_3());
    }
  } else {
    throw std::invalid_argument("Input geometry must be a Point or LineString");
  }
}

auto
Buffer3D::compute(BufferType type) const -> std::unique_ptr<PolyhedralSurface>
{
  if (_inputPoints.size() == 1) {
    return computePointBuffer();
  }

  switch (type) {
  case ROUND:
    return computeRoundBuffer();
  case CYLSPHERE:
    return computeCylSphereBuffer();
  case FLAT:
    return computeFlatBuffer();
  default:
    throw std::invalid_argument("Invalid buffer type");
  }
}

auto
Buffer3D::computePointBuffer() const -> std::unique_ptr<PolyhedralSurface>
{
  Kernel::Point_3 center(_inputPoints[0].x(), _inputPoints[0].y(),
                         _inputPoints[0].z());
  // Convert segments to subdivision level for icosahedron
  unsigned int subdivision_level =
      std::max(1U, static_cast<unsigned int>(_segments / 16));
  Sphere sphere(_radius, center, subdivision_level);
  return std::make_unique<PolyhedralSurface>(
      sphere.generatePolyhedralSurface());
}

auto
Buffer3D::computeRoundBuffer() const -> std::unique_ptr<PolyhedralSurface>
{
  using point_iterator = Point_3 *;
  using point_range    = std::pair<point_iterator, point_iterator>;
  using polyline       = std::list<point_range>;
  using Nef_polyhedron = CGAL::Nef_polyhedron_3<Kernel>;

  // Create sphere
  Point_3 center(0, 0, 0);
  // Convert segments to subdivision level for icosahedron
  unsigned int subdivision_level =
      std::max(1U, static_cast<unsigned int>(_segments / 16));
  SFCGAL::Sphere sphere(_radius, center, subdivision_level);

  // Generate polyhedron from sphere
  CGAL::Polyhedron_3<Kernel> spherePolyhedron = sphere.generatePolyhedron();

  // Convert Polyhedron to a Nef_polyhedron
  Nef_polyhedron N0(spherePolyhedron);

  // Create a polyline from _inputPoints
  polyline poly;
  if (!_inputPoints.empty()) {
    // Create a non-const copy of the points
    std::vector<Point_3> points_copy(_inputPoints.begin(), _inputPoints.end());
    poly.emplace_back(&points_copy.front(), &points_copy.back() + 1);

    // Create a Nef_polyhedron from the polyline
    Nef_polyhedron N1(poly.begin(), poly.end(),
                      Nef_polyhedron::Polylines_tag());

    // Perform Minkowski sum
    Nef_polyhedron result = CGAL::minkowski_sum_3(N0, N1);

    // Convert result to SFCGAL::PolyhedralSurface
    SFCGAL::detail::MarkedPolyhedron out;
    result.convert_to_polyhedron(out);
    return std::make_unique<PolyhedralSurface>(out);
  } // If _inputPoints is empty, return an empty PolyhedralSurface
  return std::make_unique<PolyhedralSurface>();
}

auto
Buffer3D::computeCylSphereBuffer() const -> std::unique_ptr<PolyhedralSurface>
{
  using Nef_polyhedron = CGAL::Nef_polyhedron_3<Kernel>;
  Nef_polyhedron result;

  // Add a sphere at the first point of the line
  if (!_inputPoints.empty()) {
    Kernel::Point_3 start_sphere_center(
        _inputPoints[0].x(), _inputPoints[0].y(), _inputPoints[0].z());
    Kernel::Vector_3 start_direction;
    if (_inputPoints.size() > 1) {
      start_direction = _inputPoints[1] - _inputPoints[0];
    } else {
      start_direction =
          Kernel::Vector_3(0, 0, 1); // Default direction if only one point
    }

    unsigned int subdivision_level =
        std::max(1U, static_cast<unsigned int>(_segments / 16));
    Sphere start_sphere(_radius, start_sphere_center, subdivision_level,
                        start_direction);
    CGAL::Polyhedron_3<Kernel> start_sphere_poly =
        start_sphere.generatePolyhedron();
    Nef_polyhedron start_sphere_nef(start_sphere_poly);

    result = start_sphere_nef;
  }

  // Create a cylinder and spheres for each segment of the line
  for (size_t i = 0; i < _inputPoints.size() - 1; ++i) {
    // Create a cylinder between each successive point
    Kernel::Vector_3 axis(_inputPoints[i + 1].x() - _inputPoints[i].x(),
                          _inputPoints[i + 1].y() - _inputPoints[i].y(),
                          _inputPoints[i + 1].z() - _inputPoints[i].z());
    Kernel::FT      height = CGAL::sqrt(CGAL::to_double(axis.squared_length()));
    Kernel::Point_3 base(_inputPoints[i].x(), _inputPoints[i].y(),
                         _inputPoints[i].z());
    Cylinder        cyl(base, axis, _radius, height, _segments);
    CGAL::Polyhedron_3<Kernel> cyl_poly = cyl.generatePolyhedron();
    Nef_polyhedron             cyl_nef(cyl_poly);

    result = result.join(cyl_nef);

    // Add a sphere at the junctions (rounded corners)
    if (i < _inputPoints.size() - 1) {
      Kernel::Point_3  sphereCenter(_inputPoints[i + 1].x(),
                                    _inputPoints[i + 1].y(),
                                    _inputPoints[i + 1].z());
      Kernel::Vector_3 sphere_direction;

      if (i < _inputPoints.size() - 2) {
        // For intermediate points, use the bisector of the two adjacent
        // segments
        Kernel::Vector_3 prev_dir = _inputPoints[i + 1] - _inputPoints[i];
        Kernel::Vector_3 next_dir = _inputPoints[i + 2] - _inputPoints[i + 1];
        prev_dir =
            prev_dir / CGAL::sqrt(CGAL::to_double(prev_dir.squared_length()));
        next_dir =
            next_dir / CGAL::sqrt(CGAL::to_double(next_dir.squared_length()));
        sphere_direction = prev_dir + next_dir;
      } else {
        // For the last point, use the direction of the last segment
        sphere_direction = _inputPoints[i + 1] - _inputPoints[i];
      }
      sphere_direction =
          sphere_direction /
          CGAL::sqrt(CGAL::to_double(sphere_direction.squared_length()));

      unsigned int subdivision_level =
          std::max(1U, static_cast<unsigned int>(_segments / 16));
      Sphere sphere(_radius, sphereCenter, subdivision_level, sphere_direction);
      CGAL::Polyhedron_3<Kernel> sphere_poly = sphere.generatePolyhedron();
      Nef_polyhedron             sphere_nef(sphere_poly);
      result = result.join(sphere_nef);
    }
  }

  // Convert the Nef_polyhedron to Polyhedron_3
  CGAL::Polyhedron_3<Kernel> merged_mesh;
  result.convert_to_polyhedron(merged_mesh);

  // Clean up the geometry
  PMP::remove_connected_components_of_negligible_size(merged_mesh);

  // Convert the merged mesh to PolyhedralSurface and return
  auto resultSurface = std::make_unique<PolyhedralSurface>();
  resultSurface->addPatches(merged_mesh);

  return resultSurface;
}

auto
Buffer3D::computeFlatBuffer() const -> std::unique_ptr<PolyhedralSurface>
{
  if (_inputPoints.size() < 2) {
    return std::make_unique<PolyhedralSurface>();
  }

  // Create LineString from input points
  std::vector<Point> points;
  points.reserve(_inputPoints.size());
  for (const auto &point : _inputPoints) {
    points.emplace_back(point.x(), point.y(), point.z());
  }
  LineString path(points);

  // Create circular profile using the sweep helper function
  auto profile = create_circular_profile(_radius, _segments);

  // Configure sweep options
  SweepOptions options;
  options.frame_method = SweepOptions::FrameMethod::ROTATION_MINIMIZING;
  options.start_cap    = SweepOptions::EndCapStyle::FLAT;
  options.end_cap      = SweepOptions::EndCapStyle::FLAT;
  options.closed_path  = false; // Will be auto-detected

  return sweep(path, *profile, options);
}

auto
Buffer3D::extend_point(const CGAL::Point_3<Kernel>  &point,
                       const CGAL::Vector_3<Kernel> &direction,
                       double distance) const -> CGAL::Point_3<Kernel>
{
  return point + direction * distance;
}

auto
Buffer3D::create_circle_points(const CGAL::Point_3<Kernel>  &center,
                               const CGAL::Vector_3<Kernel> &axis,
                               double radius, int segments) const
    -> std::vector<CGAL::Point_3<Kernel>>
{
  std::vector<CGAL::Point_3<Kernel>> points;
  CGAL::Vector_3<Kernel>             perpendicular =
      CGAL::cross_product(axis, CGAL::Vector_3<Kernel>(0, 0, 1));
  if (perpendicular == CGAL::NULL_VECTOR) {
    perpendicular = CGAL::cross_product(axis, Kernel::Vector_3(0, 1, 0));
  }
  perpendicular = normalizeVector(perpendicular);
  Kernel::Vector_3 perpendicular2 =
      normalizeVector(CGAL::cross_product(axis, perpendicular));

  for (int i = 0; i < segments; ++i) {
    double           angle  = 2.0 * M_PI * i / segments;
    Kernel::Vector_3 offset = radius * (std::cos(angle) * perpendicular +
                                        std::sin(angle) * perpendicular2);
    points.push_back(center + offset);
  }
  return points;
}

auto
Buffer3D::compute_bisector_plane(const CGAL::Point_3<Kernel> &p1,
                                 const CGAL::Point_3<Kernel> &p2,
                                 const CGAL::Point_3<Kernel> &p3) const
    -> CGAL::Plane_3<Kernel>
{
  CGAL::Vector_3<Kernel> v1       = normalizeVector(p2 - p1);
  CGAL::Vector_3<Kernel> v2       = normalizeVector(p3 - p2);
  CGAL::Vector_3<Kernel> bisector = v1 + v2;
  return {p2, bisector};
}

auto
Buffer3D::intersect_segment_plane(const CGAL::Point_3<Kernel> &p1,
                                  const CGAL::Point_3<Kernel> &p2,
                                  const CGAL::Plane_3<Kernel> &plane) const
    -> CGAL::Point_3<Kernel>
{
  CGAL::Vector_3<Kernel> v = p2 - p1;
  CGAL::Epeck::FT        t =
      -plane.a() * p1.x() - plane.b() * p1.y() - plane.c() * p1.z() - plane.d();
  t /= plane.a() * v.x() + plane.b() * v.y() + plane.c() * v.z();
  return p1 + t * v;
}

} // namespace SFCGAL::algorithm
