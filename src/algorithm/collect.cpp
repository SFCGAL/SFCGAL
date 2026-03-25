// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/collect.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include <memory>

namespace SFCGAL::algorithm {
/// @private
auto
collect(const Geometry &ga, const Geometry &gb) -> std::unique_ptr<Geometry>
{
  if (ga.geometryTypeId() == gb.geometryTypeId()) {
    if (ga.geometryTypeId() == TYPE_POINT) {
      auto mp = std::make_unique<MultiPoint>();
      mp->addGeometry(ga);
      mp->addGeometry(gb);
      return mp;
    }
    if (ga.geometryTypeId() == TYPE_LINESTRING) {
      auto mls = std::make_unique<MultiLineString>();
      mls->addGeometry(ga);
      mls->addGeometry(gb);
      return mls;
    }
    if (ga.geometryTypeId() == TYPE_POLYGON) {
      auto mp = std::make_unique<MultiPolygon>();
      mp->addGeometry(ga);
      mp->addGeometry(gb);
      return mp;
    }
    if (ga.geometryTypeId() == TYPE_SOLID) {
      auto mp = std::make_unique<MultiSolid>();
      mp->addGeometry(ga);
      mp->addGeometry(gb);
      return mp;
    }
  }

  // else
  auto coll = std::make_unique<GeometryCollection>();
  coll->addGeometry(ga);
  coll->addGeometry(gb);
  return coll;
}
} // namespace SFCGAL::algorithm
