// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ENVELOPE_H_
#define SFCGAL_ENVELOPE_H_

#include <boost/assert.hpp>
#include <memory>

#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>

#include "SFCGAL/config.h"

#include "SFCGAL/Coordinate.h"
#include "SFCGAL/detail/Interval.h"

namespace SFCGAL {

class LineString;
class Polygon;
class Solid;
class PolyhedralSurface;

/**
 * Represents a bounding box
 * @todo add asText instead of "print"?
 * @todo add basic operations (especially intersects/intersects3D)
 */
class SFCGAL_API Envelope {
public:
  /**
   * default constructor (empty bounding box)
   */
  Envelope();
  /**
   * 2D box constructor with min,max values
   * @param xmin Minimum X value
   * @param xmax Maximum X value
   * @param ymin Minimum Y value
   * @param ymax Maximum Y value
   */
  Envelope(const double &xmin, const double &xmax, const double &ymin,
           const double &ymax);
  /**
   * 3D box constructor with min,max values
   * @param xmin Minimum X value
   * @param xmax Maximum X value
   * @param ymin Minimum Y value
   * @param ymax Maximum Y value
   * @param zmin Minimum Z value
   * @param zmax Maximum Z value
   */
  Envelope(const double &xmin, const double &xmax, const double &ymin,
           const double &ymax, const double &zmin, const double &zmax);
  /**
   * Constructor from single coordinate
   * @param coordinate The coordinate to create envelope from
   */
  Envelope(const Coordinate &coordinate);
  /**
   * Constructor from two coordinates
   * @param coordinate1 First coordinate
   * @param coordinate2 Second coordinate
   */
  Envelope(const Coordinate &coordinate1, const Coordinate &coordinate2);
  /**
   * copy constructor
   * @param other The envelope to copy from
   */
  Envelope(const Envelope &other);
  /**
   * assign operator
   * @param other The envelope to assign from
   * @return Reference to this envelope
   */
  auto
  operator=(const Envelope &other) -> Envelope &;

  /**
   * indicates if the bounding box is empty
   * @return true if empty, false otherwise
   */
  [[nodiscard]] auto
  isEmpty() const -> bool;
  /**
   * indicates if the bounding box has a 3D component
   * @return true if 3D, false otherwise
   */
  [[nodiscard]] auto
  is3D() const -> bool;

  /**
   * expand the box to include coordinate
   * @param coordinate The coordinate to include in the envelope
   */
  void
  expandToInclude(const Coordinate &coordinate);

  /**
   * @brief Get minimum X value
   * @return The minimum X coordinate
   */
  [[nodiscard]] auto
  xMin() const -> const double &
  {
    return _bounds[0].lower();
  }
  /**
   * @brief Get minimum Y value
   * @return The minimum Y coordinate
   */
  auto
  yMin() const -> const double &
  {
    return _bounds[1].lower();
  }
  /**
   * @brief Get minimum Z value
   * @return The minimum Z coordinate
   */
  auto
  zMin() const -> const double &
  {
    return _bounds[2].lower();
  }

  /**
   * @brief Get maximum X value
   * @return The maximum X coordinate
   */
  auto
  xMax() const -> const double &
  {
    return _bounds[0].upper();
  }
  /**
   * @brief Get maximum Y value
   * @return The maximum Y coordinate
   */
  auto
  yMax() const -> const double &
  {
    return _bounds[1].upper();
  }
  /**
   * @brief Get maximum Z value
   * @return The maximum Z coordinate
   */
  auto
  zMax() const -> const double &
  {
    return _bounds[2].upper();
  }

  /**
   * returns the n-th bound
   * @param n The index of the bound (0=X, 1=Y, 2=Z)
   * @return Reference to the nth interval bound
   */
  auto
  boundsN(const size_t &n) -> detail::Interval &
  {
    BOOST_ASSERT(n < 3);
    return _bounds[n];
  }
  /**
   * returns the n-th bound
   * @param n The index of the bound (0=X, 1=Y, 2=Z)
   * @return Const reference to the nth interval bound
   */
  auto
  boundsN(const size_t &n) const -> const detail::Interval &
  {
    BOOST_ASSERT(n < 3);
    return _bounds[n];
  }

  /**
   * Convenience function. Convert to CGAL::BBox_2
   * @return CGAL 2D bounding box
   */
  auto
  toBbox_2() const -> CGAL::Bbox_2
  {
    BOOST_ASSERT(!isEmpty());

    return {_bounds[0].lower(), _bounds[1].lower(), _bounds[0].upper(),
            _bounds[1].upper()};
  }

  /**
   * Convenience function. Convert to CGAL::BBox_3
   * @return CGAL 3D bounding box
   */
  auto
  toBbox_3() const -> CGAL::Bbox_3
  {
    if (is3D()) {
      return {_bounds[0].lower(), _bounds[1].lower(), _bounds[2].lower(),
              _bounds[0].upper(), _bounds[1].upper(), _bounds[2].upper()};
    }

    return {_bounds[0].lower(), _bounds[1].lower(), 0.0,
            _bounds[0].upper(), _bounds[1].upper(), 0.0};
  }

  /**
   * @brief Global binary operator on Envelopes. Test if A's bounding box
   * contains B's
   * @param envelopeA First envelope
   * @param envelopeB Second envelope
   * @return true if envelope envelopeA contains envelope envelopeB
   * FIXME: consider moving that outside of the class
   */
  static auto
  contains(const Envelope &envelopeA, const Envelope &envelopeB) -> bool;

  /**
   * @brief Global binary operator on Envelopes. Test if A's bounding box
   * overlaps B's
   * @param envelopeA First envelope
   * @param envelopeB Second envelope
   * @return true if envelope envelopeA overlaps envelope envelopeB
   */
  static auto
  overlaps(const Envelope &envelopeA, const Envelope &envelopeB) -> bool;

  //-- helpers

  /**
   * @brief convenience method to convert to 2D Polygon ring
   * @return A LineString representing the envelope boundary
   * @warning empty LineString for empty Envelope, may be X or Y collapsed
   */
  auto
  toRing() const -> std::unique_ptr<LineString>;
  /**
   * @brief convenience method to convert to 2D Polygon
   * @return A Polygon representing the envelope
   * @warning empty Polygon for empty Envelope, may be X or Y collapsed
   */
  auto
  toPolygon() const -> std::unique_ptr<Polygon>;

  /**
   * @brief convenience method to convert to 3D Shell
   * @return A PolyhedralSurface representing the envelope shell
   * @warning empty Solid for empty or non 3D Envelope, may be X, Y or Z
   * collapsed
   */
  auto
  toShell() const -> std::unique_ptr<PolyhedralSurface>;
  /**
   * @brief convenience method to convert to 3D Solid
   * @return A Solid representing the envelope
   * @warning empty Solid for empty or non 3D Envelope, may be X, Y or Z
   * collapsed
   */
  auto
  toSolid() const -> std::unique_ptr<Solid>;

  /**
   * @brief Display method
   * @param ostr Output stream to print to
   * @return Reference to the output stream
   */
  auto
  print(std::ostream &ostr) const -> std::ostream &;

private:
  /**
   * bounds of the interval ((xmin,xmax),(ymin,ymax),(zmin,zmax))
   */
  detail::Interval _bounds[3];
};

/**
 * @brief Global comparison operator on Envelope
 * @param lhs Left envelope to compare
 * @param rhs Right envelope to compare
 * @return true if envelopes are equal, false otherwise
 */
SFCGAL_API auto
operator==(const Envelope &lhs, const Envelope &rhs) -> bool;

} // namespace SFCGAL

#endif
