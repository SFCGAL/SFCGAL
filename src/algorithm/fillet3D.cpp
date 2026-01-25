// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/fillet3D.h"

#include "SFCGAL/algorithm/convexHull.h"
#include "SFCGAL/algorithm/difference.h"
#include "SFCGAL/algorithm/minkowskiSum3D.h"
#include "SFCGAL/algorithm/normal.h"
#include "SFCGAL/algorithm/translate.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"

#include <CGAL/number_utils.h>
#include <cmath>
#include <optional>
#include <vector>

namespace SFCGAL::algorithm {

namespace {

// Get incident face vectors for an edge segment
auto
get_incident_face_vectors(const SFCGAL::Solid      &solid,
                          const SFCGAL::LineString &segment)
    -> std::vector<SFCGAL::Kernel::Vector_3>
{
  std::vector<SFCGAL::Kernel::Vector_3> vectors;
  if (segment.numPoints() < 2)
    return vectors;

  const auto      &p1     = segment.pointN(0);
  const auto      &p2     = segment.pointN(1);
  Kernel::Vector_3 segVec = p2.toPoint_3() - p1.toPoint_3();

  auto process_surface = [&](const SFCGAL::PolyhedralSurface &surface) {
    for (size_t i = 0; i < surface.numPatches(); ++i) {
      const auto &poly   = surface.patchN(i);
      const auto &ring   = poly.exteriorRing();
      bool        has_p1 = false, has_p2 = false;
      size_t      idx_p1 = 0, idx_p2 = 0;
      for (size_t j = 0; j < ring.numPoints(); ++j) {
        if (ring.pointN(j).coordinate() == p1.coordinate()) {
          has_p1 = true;
          idx_p1 = j;
        }
        if (ring.pointN(j).coordinate() == p2.coordinate()) {
          has_p2 = true;
          idx_p2 = j;
        }
      }
      if (has_p1 && has_p2) {
        size_t n_eff = ring.numPoints() > 0 ? ring.numPoints() - 1 : 0;
        if (n_eff == 0)
          continue;
        bool consecutive = false, p1_to_p2 = false;
        if ((idx_p1 + 1) % n_eff == idx_p2 % n_eff) {
          consecutive = true;
          p1_to_p2    = true;
        } else if ((idx_p2 + 1) % n_eff == idx_p1 % n_eff) {
          consecutive = true;
          p1_to_p2    = false;
        }
        if (consecutive) {
          auto normal = SFCGAL::algorithm::normal3D<Kernel>(poly, true);
          Kernel::Vector_3 inward;
          if (p1_to_p2)
            inward = CGAL::cross_product(normal, segVec);
          else
            inward = CGAL::cross_product(segVec, normal);
          double len = std::sqrt(CGAL::to_double(inward.squared_length()));
          if (len > 1e-12)
            vectors.push_back(inward / len);
        }
      }
    }
  };
  process_surface(solid.exteriorShell());
  for (size_t i = 0; i < solid.numInteriorShells(); ++i)
    process_surface(solid.interiorShellN(i));

  // Deduplicate vectors
  if (vectors.size() > 2) {
    std::vector<Kernel::Vector_3> unique_vectors;
    for (const auto &v : vectors) {
      bool found = false;
      for (const auto &uv : unique_vectors)
        if (CGAL::to_double(CGAL::scalar_product(v, uv)) > 0.99) {
          found = true;
          break;
        }
      if (!found)
        unique_vectors.push_back(v);
    }
    return unique_vectors;
  }
  return vectors;
}

// Create arc points from P1 to P2 with center C using Rodrigues rotation
auto
create_arc_points(const Kernel::Point_3 &C, const Kernel::Point_3 &P1,
                  const Kernel::Point_3 &P2, int num_segments)
    -> std::vector<Kernel::Point_3>
{
  std::vector<Kernel::Point_3> arc_points;

  Kernel::Vector_3 v1 = P1 - C;
  Kernel::Vector_3 v2 = P2 - C;

  double arc_radius = std::sqrt(CGAL::to_double(v1.squared_length()));
  if (arc_radius < 1e-10) {
    arc_points.push_back(P1);
    return arc_points;
  }

  double cos_theta = CGAL::to_double(v1 * v2) / (arc_radius * arc_radius);
  cos_theta        = std::clamp(cos_theta, -1.0, 1.0);
  double theta     = std::acos(cos_theta);

  // Rotation axis (perpendicular to plane containing v1, v2)
  Kernel::Vector_3 rot_axis = CGAL::cross_product(v1, v2);
  double rot_axis_len = std::sqrt(CGAL::to_double(rot_axis.squared_length()));
  if (rot_axis_len < 1e-10) {
    // v1 and v2 are collinear, just return endpoints
    arc_points.push_back(P1);
    arc_points.push_back(P2);
    return arc_points;
  }
  rot_axis = rot_axis / rot_axis_len;

  // Sample arc using Rodrigues' rotation formula
  for (int i = 0; i <= num_segments; ++i) {
    double t     = static_cast<double>(i) / num_segments;
    double angle = t * theta;

    double cos_a   = std::cos(angle);
    double sin_a   = std::sin(angle);
    double k_dot_v = CGAL::to_double(rot_axis * v1);

    Kernel::Vector_3 v_rot = v1 * cos_a +
                             CGAL::cross_product(rot_axis, v1) * sin_a +
                             rot_axis * k_dot_v * (1.0 - cos_a);

    arc_points.push_back(C + v_rot);
  }

  return arc_points;
}

// Create fillet polygon cross-section (triangle with arc cutout)
auto
create_fillet_polygon(const Kernel::Point_3 &E, const Kernel::Vector_3 &u1,
                      const Kernel::Vector_3 &u2, double radius,
                      int arc_segments) -> std::unique_ptr<Polygon>
{
  // Arc endpoints
  Kernel::Point_3 P1 = E + u1 * radius;
  Kernel::Point_3 P2 = E + u2 * radius;

  // Arc center (for 90-degree case)
  // For general angle, should use bisector, but 90-degree is common
  Kernel::Point_3 C = E + u1 * radius + u2 * radius;

  // Generate arc points
  auto arc = create_arc_points(C, P1, P2, arc_segments);

  // Build polygon: E -> arc points -> back to E
  LineString ring;
  ring.addPoint(Point(E));
  for (const auto &pt : arc) {
    ring.addPoint(Point(pt));
  }
  ring.addPoint(Point(E)); // Close

  return std::make_unique<Polygon>(ring);
}

} // namespace

auto
fillet3D(const Geometry &solid_geom, const Geometry &edges, double radius,
         unsigned int num_subdivisions) -> std::unique_ptr<Geometry>
{
  try {
    // Convert solid input
    const Solid           *solid = nullptr;
    std::unique_ptr<Solid> temp_solid;
    if (solid_geom.geometryTypeId() == TYPE_SOLID)
      solid = static_cast<const Solid *>(&solid_geom);
    else if (solid_geom.geometryTypeId() == TYPE_POLYHEDRALSURFACE) {
      temp_solid = std::make_unique<Solid>(
          static_cast<const PolyhedralSurface &>(solid_geom));
      solid = temp_solid.get();
    }
    if (!solid)
      throw Exception("fillet3D requires a Solid or PolyhedralSurface");

    // Extract paths
    std::vector<const LineString *> paths;
    if (edges.geometryTypeId() == TYPE_LINESTRING)
      paths.push_back(&edges.as<LineString>());
    else if (edges.geometryTypeId() == TYPE_MULTILINESTRING) {
      const auto &mls = edges.as<MultiLineString>();
      for (size_t i = 0; i < mls.numGeometries(); ++i)
        paths.push_back(&(mls.geometryN(i).as<LineString>()));
    }

    // Arc segments
    int arc_segments = 8 * (num_subdivisions + 1);

    // Build fillet volumes using Minkowski sum approach
    std::vector<std::unique_ptr<Geometry>> parts;

    struct Interface {
      Point                        p_base;
      std::vector<Kernel::Point_3> arc_points;
      Kernel::Vector_3             u1, u2;
    };
    std::optional<Interface> prev_interface;

    for (const auto *path : paths) {
      if (path->isEmpty())
        continue;

      // Reset interface if paths are disconnected
      if (prev_interface &&
          path->pointN(0).coordinate() != prev_interface->p_base.coordinate())
        prev_interface = std::nullopt;

      for (size_t i = 0; i < path->numSegments(); ++i) {
        Point      p_start = path->pointN(i), p_end = path->pointN(i + 1);
        LineString segment;
        segment.addPoint(p_start);
        segment.addPoint(p_end);

        // Get incident face vectors
        auto dirs = get_incident_face_vectors(*solid, segment);
        if (dirs.size() < 2) {
          prev_interface = std::nullopt;
          continue;
        }

        Kernel::Vector_3 u1 = dirs[0], u2 = dirs[1];

        // Ensure consistency with previous segment
        if (prev_interface) {
          double dot11 = CGAL::to_double(prev_interface->u1 * u1);
          double dot12 = CGAL::to_double(prev_interface->u1 * u2);
          if (dot12 > dot11)
            std::swap(u1, u2);
        }

        // Create fillet polygon at start
        auto fillet_polygon = create_fillet_polygon(p_start.toPoint_3(), u1, u2,
                                                    radius, arc_segments);

        // Create path for Minkowski sum (segment as LineString from origin)
        Kernel::Vector_3 segVec = p_end.toPoint_3() - p_start.toPoint_3();
        LineString       path_relative;
        path_relative.addPoint(Point(0, 0, 0));
        path_relative.addPoint(Point(segVec.x(), segVec.y(), segVec.z()));

        // Translate fillet polygon to origin
        Kernel::Vector_3 to_origin = CGAL::ORIGIN - p_start.toPoint_3();
        std::unique_ptr<Geometry> poly_at_origin = fillet_polygon->clone();
        algorithm::translate(*poly_at_origin, to_origin);

        // Minkowski sum: sweep fillet polygon along segment
        auto prism = algorithm::minkowskiSum3D(*poly_at_origin, path_relative);

        // Translate back to actual position
        Kernel::Vector_3 from_origin = p_start.toPoint_3() - CGAL::ORIGIN;
        algorithm::translate(*prism, from_origin);

        if (!prism->isEmpty())
          parts.push_back(std::move(prism));

        // Corner handling via convex hull
        if (prev_interface) {
          // Generate current arc points
          Kernel::Point_3 C_start =
              p_start.toPoint_3() + u1 * radius + u2 * radius;
          Kernel::Point_3 P1_start = p_start.toPoint_3() + u1 * radius;
          Kernel::Point_3 P2_start = p_start.toPoint_3() + u2 * radius;
          auto            current_arc =
              create_arc_points(C_start, P1_start, P2_start, arc_segments);

          // Create MultiPoint with all corner points
          MultiPoint mp;
          mp.addGeometry(p_start);
          for (const auto &pt : prev_interface->arc_points)
            mp.addGeometry(Point(pt));
          for (const auto &pt : current_arc)
            mp.addGeometry(Point(pt));

          try {
            auto corner = algorithm::convexHull3D(mp);
            if (corner->geometryTypeId() == TYPE_POLYHEDRALSURFACE ||
                corner->geometryTypeId() == TYPE_SOLID) {
              if (corner->geometryTypeId() == TYPE_POLYHEDRALSURFACE)
                parts.push_back(std::make_unique<Solid>(
                    *static_cast<PolyhedralSurface *>(corner.get())));
              else
                parts.push_back(std::move(corner));
            }
          } catch (...) {
            // Convex hull failed, skip corner
          }
        }

        // Store interface for next segment
        Kernel::Point_3 C_end  = p_end.toPoint_3() + u1 * radius + u2 * radius;
        Kernel::Point_3 P1_end = p_end.toPoint_3() + u1 * radius;
        Kernel::Point_3 P2_end = p_end.toPoint_3() + u2 * radius;
        auto end_arc   = create_arc_points(C_end, P1_end, P2_end, arc_segments);
        prev_interface = Interface{p_end, end_arc, u1, u2};
      }
    }

    // Union all parts into a single geometry collection
    if (parts.empty())
      return solid_geom.clone();

    GeometryCollection fillet_volume;
    for (auto &part : parts)
      fillet_volume.addGeometry(part.release());

    // Subtract fillet volume from solid
    return difference3D_nef(solid_geom, fillet_volume);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] fillet3D exception: " << e.what() << std::endl;
    throw;
  } catch (...) {
    std::cerr << "[ERROR] fillet3D unknown exception" << std::endl;
    throw;
  }
}

auto
fillet3D_asymmetric(const Geometry &solid, const Geometry &edges,
                    double distance1, double distance2,
                    unsigned int num_subdivisions) -> std::unique_ptr<Geometry>
{
  // For asymmetric fillet, use average radius for now
  // TODO: Implement proper asymmetric fillet with different offsets
  double radius = std::sqrt(distance1 * distance1 + distance2 * distance2);
  return fillet3D(solid, edges, radius, num_subdivisions);
}

} // namespace SFCGAL::algorithm
