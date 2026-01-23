// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/orientedPrism.h"

#include "SFCGAL/algorithm/convexHull.h"
#include "SFCGAL/algorithm/minkowskiSum3D.h"
#include "SFCGAL/algorithm/normal.h"
#include "SFCGAL/algorithm/translate.h"
#include "SFCGAL/algorithm/union.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"

#include <CGAL/number_utils.h>
#include <iostream>
#include <optional>
#include <vector>

namespace SFCGAL::algorithm {

namespace {

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

} // namespace

auto
orientedPrism(const Geometry &solid_geom, const Geometry &path_geom, double d1,
              double d2) -> std::unique_ptr<Geometry>
{
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
    throw Exception("orientedPrism requires a Solid or PolyhedralSurface");

  std::vector<const LineString *> paths;
  if (path_geom.geometryTypeId() == TYPE_LINESTRING)
    paths.push_back(static_cast<const LineString *>(&path_geom));
  else if (path_geom.geometryTypeId() == TYPE_MULTILINESTRING) {
    const auto &mls = static_cast<const MultiLineString &>(path_geom);
    for (size_t i = 0; i < mls.numGeometries(); ++i)
      paths.push_back(&(mls.geometryN(i).as<LineString>()));
  }

  std::vector<std::unique_ptr<Geometry>> parts;
  struct Interface {
    Point            p_base, p1, p2;
    Kernel::Vector_3 u1, u2;
  };
  std::optional<Interface> prev_interface;

  for (const auto *path : paths) {
    if (path->isEmpty())
      continue;
    if (prev_interface &&
        path->pointN(0).coordinate() != prev_interface->p_base.coordinate())
      prev_interface = std::nullopt;
    for (size_t i = 0; i < path->numSegments(); ++i) {
      Point      p_start = path->pointN(i), p_end = path->pointN(i + 1);
      LineString segment;
      segment.addPoint(p_start);
      segment.addPoint(p_end);
      auto dirs = get_incident_face_vectors(*solid, segment);
      if (dirs.size() < 2) {
        prev_interface = std::nullopt;
        continue;
      }
      Kernel::Vector_3 u1 = dirs[0], u2 = dirs[1];
      if (prev_interface) {
        double dot11 = CGAL::to_double(prev_interface->u1 * u1);
        double dot12 = CGAL::to_double(prev_interface->u1 * u2);
        if (dot12 > dot11)
          std::swap(u1, u2);
      }
      Point p1_start(p_start.toPoint_3() + u1 * d1),
          p2_start(p_start.toPoint_3() + u2 * d2);
      Triangle         tri_start(p_start, p1_start, p2_start);
      Kernel::Vector_3 segVec = p_end.toPoint_3() - p_start.toPoint_3();
      LineString       path_relative;
      path_relative.addPoint(Point(0, 0, 0));
      path_relative.addPoint(Point(segVec.x(), segVec.y(), segVec.z()));
      Kernel::Vector_3          to_origin = CGAL::ORIGIN - p_start.toPoint_3();
      std::unique_ptr<Geometry> tri_at_origin = tri_start.clone();
      algorithm::translate(*tri_at_origin, to_origin);
      auto prism = algorithm::minkowskiSum3D(*tri_at_origin, path_relative);
      Kernel::Vector_3 from_origin = p_start.toPoint_3() - CGAL::ORIGIN;
      algorithm::translate(*prism, from_origin);
      if (!prism->isEmpty())
        parts.push_back(std::move(prism));
      if (prev_interface) {
        MultiPoint mp;
        mp.addGeometry(p_start);
        mp.addGeometry(prev_interface->p1);
        mp.addGeometry(prev_interface->p2);
        mp.addGeometry(p1_start);
        mp.addGeometry(p2_start);
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
        }
      }
      prev_interface = Interface{p_end, Point(p_end.toPoint_3() + u1 * d1),
                                 Point(p_end.toPoint_3() + u2 * d2), u1, u2};
    }
  }

  if (parts.empty())
    return std::make_unique<GeometryCollection>();

  // Return MultiSolid but avoid union3D loop for now
  auto multi = std::make_unique<MultiSolid>();
  for (auto &part : parts) {
    if (part->geometryTypeId() == TYPE_SOLID)
      multi->addGeometry(static_cast<Solid *>(part.release()));
    else if (part->geometryTypeId() == TYPE_POLYHEDRALSURFACE)
      multi->addGeometry(
          new Solid(*static_cast<PolyhedralSurface *>(part.get())));
  }
  return multi;
}

} // namespace SFCGAL::algorithm
