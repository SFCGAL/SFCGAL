#ifndef SFCGAL_PRIMITIVE_SNAP_H_
#define SFCGAL_PRIMITIVE_SNAP_H_

#include <cmath>

#include "SFCGAL/Kernel.h"
#include "SFCGAL/Point.h"

namespace SFCGAL {

inline constexpr double SNAP_GRID = 1e8;

auto inline snapCoordinate(const Kernel::FT &value) -> Kernel::FT
{
  return {std::round(CGAL::to_double(value) * SNAP_GRID) / SNAP_GRID};
}

void inline snapPoint(Point &point)
{
  point = Point(snapCoordinate(point.x()), snapCoordinate(point.y()),
                snapCoordinate(point.z()), point.m());
}

} // namespace SFCGAL

#endif // SFCGAL_PRIMITIVE_SNAP_H_
