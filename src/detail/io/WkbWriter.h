// Copyright (c) 2023-2024, Oslandia.
// Copyright (c) 2024-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_IO_WKBWRITER_H_
#define SFCGAL_IO_WKBWRITER_H_

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

#include <boost/endian/conversion.hpp>

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"

namespace SFCGAL {
// Forward declarations
class Point;
class LineString;
class Polygon;
class Triangle;
class GeometryCollection;
class PolyhedralSurface;
class TriangulatedSurface;
class NURBSCurve;

// Use srid_t from SFCGAL/config.h
} // namespace SFCGAL

namespace SFCGAL::detail::io {

/**
 * Writer for WKB
 *
 */
class SFCGAL_API WkbWriter {
public:
  /**
   * @brief Construct WKB writer
   * @param outputStream Output stream to write WKB data to
   * @param asHexString If true, output as hexadecimal string, otherwise binary
   */
  WkbWriter(std::ostream &outputStream, bool asHexString = false)
      : _s(outputStream), _asHexString(asHexString) {};

  /**
   * write WKB for a geometry
   * wkbOrder is the native endianness by default.
   * @param geometry Geometry to write as WKB
   * @param wkbOrder Byte order for WKB output
   */
  void
  write(const Geometry &geometry, std::endian wkbOrder = std::endian::native);

  /**
   * write EWKB for a geometry
   * wkbOrder is the native endianness by default.
   * @param geometry Geometry to write as EWKB
   * @param srid Spatial reference identifier to include in EWKB
   * @param wkbOrder Byte order for WKB output
   */
  void
  write(const Geometry &geometry, const srid_t &srid,
        std::endian wkbOrder = std::endian::native);

private:
  /**
   * Dedicated method to write the geometry type into _wkb data
   */
  void
  writeGeometryType(const Geometry &geometry,
                    std::endian     wkbOrder = std::endian::native);

  /**
   * Dedicated method to write Point into _wkb data
   */
  void
  writeInner(const Point &geometry, std::endian wkbOrder = std::endian::native);

  /**
   * Dedicated method to write LineString into _wkb data
   */
  void
  writeInner(const LineString &geometry,
             std::endian       wkbOrder = std::endian::native);
  /**
   * Dedicated method to write Ring into _wkb data
   *
   * This method is shared by LineString and Polygon.
   */
  void
  writeInnerRing(const LineString &geometry,
                 std::endian       wkbOrder = std::endian::native);

  /**
   * Dedicated method to write Polygon into _wkb data
   */
  void
  writeInner(const Polygon &geometry,
             std::endian    wkbOrder = std::endian::native);

  /**
   * Dedicated method to write GeometryCollection into _wkb data
   */
  void
  writeInner(const GeometryCollection &geometry,
             std::endian               wkbOrder = std::endian::native);

  /**
   * Dedicated method to write PolyhedralSurface into _wkb data
   */
  void
  writeInner(const PolyhedralSurface &polyhedralSurface,
             std::endian              wkbOrder = std::endian::native);

  /**
   * Dedicated method to write TriangulatedSurface into _wkb data
   */
  void
  writeInner(const TriangulatedSurface &triangulatedSurface,
             std::endian                wkbOrder = std::endian::native);

  /**
   * Dedicated method to write Multi Geometries into _wkb data
   * Multi Geometries are: MultiPoint, MultiLineString, MultiPolygon.
   */
  template <typename M, typename G>
  void
  writeInner(const M &geometry, std::endian wkbOrder);

  /**
   * Dedicated method to write Triangle into _wkb data
   */
  void
  writeInner(const Triangle &geometry,
             std::endian     wkbOrder = std::endian::native);

  /**
   * Dedicated method to write NURBSCurve into _wkb data
   */
  void
  writeInner(const NURBSCurve &geometry,
             std::endian       wkbOrder = std::endian::native);

  /**
   * Dedicated method to write Point into _wkb data
   */
  void
  writeCoordinate(const Point &geometry,
                  std::endian  wkbOrder = std::endian::native);

  /**
   * Method to write Geometry into _wkb data
   * Only for recursive call use
   */
  void
  writeRec(const Geometry &geometry,
           std::endian     wkbOrder = std::endian::native);

  std::ostream &_s;

  /**
   * if false, will read the wkb as binary string, else will read the string as
   * hex ascii string (ie. 2 chars for 1 byte with values matching [0-9A-F]
   */
  bool _asHexString;

  template <std::size_t N>
  auto
  toStream(const std::array<std::byte, N> &arr) -> void
  {
    if (_asHexString) {
      for (const std::byte &byteVal : arr) {
        _s << _prefix << std::hex << std::setw(2) << std::setfill('0')
           << static_cast<int>(byteVal);
      }
    } else {
      for (const std::byte &byteVal : arr) {
        _s << static_cast<unsigned char>(byteVal);
      }
    }
  }

  template <typename T>
  auto
  toByte(const T value, std::endian byteOrder) -> void
  {
    T valueSwapped = value;
    if (std::endian::native != byteOrder) {
      boost::endian::endian_reverse_inplace(valueSwapped);
    }
    toStream(
        *reinterpret_cast<std::array<std::byte, sizeof(T)> *>(&valueSwapped));
  }

  srid_t _srid;

  bool        _useSrid = false;
  bool        _isEWKB  = false;
  std::string _prefix;
};

} // namespace SFCGAL::detail::io

#endif
