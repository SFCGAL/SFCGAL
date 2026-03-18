// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "split.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/detail/algorithm/split.h"

#include <memory>

namespace SFCGAL::algorithm {

SFCGAL_API auto
split(const Geometry &geometry, const Point &planePoint,
      const Kernel::Vector_3 &planeNormal, bool closeGeometries)
    -> std::unique_ptr<GeometryCollection>
{

  const CGAL::Plane_3<Kernel> plane(planePoint.toPoint_3(), planeNormal);

  switch (geometry.geometryTypeId()) {
  case TYPE_GEOMETRYCOLLECTION: {
    // call split for all the geometries in the collection
    auto result = std::make_unique<GeometryCollection>();
    for (const auto &internalGeometry : geometry.as<GeometryCollection>()) {
      try {
        std::unique_ptr<GeometryCollection> newGeoms =
            split(internalGeometry, planePoint, planeNormal, closeGeometries);
        for (const auto &newGeom : *newGeoms) {
          result->addGeometry(newGeom);
        }
      } catch (Exception &e) {
        std::cerr << e.what() << "\n";
      }
    }
    return result;
  }
  case TYPE_MULTIPOLYGON:
    return detail::split(geometry.as<MultiPolygon>(), plane, closeGeometries);
  case TYPE_MULTISOLID:
    return detail::split(geometry.as<MultiSolid>(), plane, closeGeometries);
  case TYPE_POLYGON:
    return detail::split(geometry.as<Polygon>(), plane, closeGeometries);
  case TYPE_POLYHEDRALSURFACE:
    return detail::split(geometry.as<PolyhedralSurface>(), plane,
                         closeGeometries);
  case TYPE_TRIANGULATEDSURFACE:
    return detail::split(geometry.as<TriangulatedSurface>(), plane,
                         closeGeometries);
  case TYPE_SOLID:
    // A solid is always closed
    return detail::split(geometry.as<Solid>(), plane, true);

  case TYPE_LINESTRING:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOINT:
  case TYPE_NURBSCURVE:
  case TYPE_POINT:
  case TYPE_TRIANGLE:
    break;
  }

  BOOST_THROW_EXCEPTION(Exception(
      (boost::format(
           "Unexpected geometry type (%s) in SFCGAL::algorithm::split") %
       geometry.geometryType())
          .str()));
}

} // namespace SFCGAL::algorithm
