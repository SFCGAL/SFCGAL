// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/split3D.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/detail/algorithm/split3D.h"

#include <memory>

namespace SFCGAL::algorithm {

SFCGAL_API auto
split3D(const Geometry &geometry, const Point &planePoint,
        const Kernel::Vector_3 &planeNormal, bool closeGeometries)
    -> std::unique_ptr<GeometryCollection>
{

  const CGAL::Plane_3<Kernel> plane(planePoint.toPoint_3(), planeNormal);

  switch (geometry.geometryTypeId()) {
  case TYPE_POLYHEDRALSURFACE:
    return detail::split3D(geometry.as<PolyhedralSurface>(), plane,
                           closeGeometries);
  case TYPE_TRIANGULATEDSURFACE:
    return detail::split3D(geometry.as<TriangulatedSurface>(), plane,
                           closeGeometries);
  case TYPE_SOLID:
    // A solid is always closed
    return detail::split3D(geometry.as<Solid>(), plane, true);

  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_LINESTRING:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOINT:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_NURBSCURVE:
  case TYPE_POINT:
  case TYPE_POLYGON:
  case TYPE_TRIANGLE:
    break;
  }

  BOOST_THROW_EXCEPTION(Exception(
      (boost::format(
           "Unexpected geometry type (%s) in SFCGAL::algorithm::split3D") %
       geometry.geometryType())
          .str()));
}

} // namespace SFCGAL::algorithm
