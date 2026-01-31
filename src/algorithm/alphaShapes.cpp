// Copyright (c) 2022-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/alphaShapes.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Polygon.h"

#include "SFCGAL/detail/GetPointsVisitor.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/Kernel.h"
#include <boost/format.hpp>

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>

#include <CGAL/Arr_non_caching_segment_basic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <memory>
#include <vector>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/** @cond */

using Vb              = CGAL::Alpha_shape_vertex_base_2<Kernel>;
using Fb              = CGAL::Alpha_shape_face_base_2<Kernel>;
using Tds             = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using Triangulation_2 = CGAL::Delaunay_triangulation_2<Kernel, Tds>;
using Alpha_shape_2   = CGAL::Alpha_shape_2<Triangulation_2>;

using Alpha_shape_edges_iterator = Alpha_shape_2::Alpha_shape_edges_iterator;

using Traits_2    = CGAL::Arr_non_caching_segment_basic_traits_2<Kernel>;
using Arrangement = CGAL::Arrangement_2<Traits_2>;

template <class OutputIterator>
void
alphaEdges(const Alpha_shape_2 &alphaShape, OutputIterator out)
{
  // keep only stable boundary of the alpha-shape
  auto it  = alphaShape.alpha_shape_edges_begin();
  auto end = alphaShape.alpha_shape_edges_end();
  for (; it != end; ++it) {
    if (alphaShape.classify(*it) == Alpha_shape_2::REGULAR) {
      *out++ = alphaShape.segment(*it);
    }
  }
}

static auto
computeAlpha(const Geometry &geometry, Alpha_shape_2 &alphaShape,
             double alpha = 0, size_t nbComponents = 1) -> double
{
  using CGAL::object_cast;
  double result = -1.0;

  if (geometry.isEmpty()) {
    return result;
  }

  SFCGAL::detail::GetPointsVisitor getPointVisitor;
  const_cast<Geometry &>(geometry).accept(getPointVisitor);

  // collect points
  if (getPointVisitor.points.size() < 4) {
    return result;
  }

  std::vector<Point_2> points;

  points.reserve(getPointVisitor.points.size());
  for (auto &point : getPointVisitor.points) {
    points.push_back(point->toPoint_2());
  }

  alphaShape.make_alpha_shape(points.begin(), points.end());
  alphaShape.set_alpha(Kernel::FT(alpha));

  // Ensure we compare the iterator from find_optimal_alpha(nbComponents)
  // against alpha_end() before dereferencing to avoid potential crash when
  // no valid alpha is found.
  auto it_alpha = alphaShape.find_optimal_alpha(nbComponents);
  if (it_alpha != alphaShape.alpha_end()) {
    result = CGAL::to_double(*it_alpha);
  } else {
    BOOST_THROW_EXCEPTION(Exception("Can't find alpha value."));
  }

  return result;
}

static auto
alphaToGeometry(const Alpha_shape_2 &alphaShape, bool allowHoles)
    -> std::unique_ptr<Geometry>
{
  std::vector<Segment_2> segments;
  alphaEdges(alphaShape, std::back_inserter(segments));

  Arrangement arr;

  CGAL::insert_non_intersecting_curves(arr, segments.begin(), segments.end());
  auto poly{std::make_unique<Polygon>()};
  for (auto f = arr.faces_begin(); f != arr.faces_end(); f++) {
    auto ring{std::make_unique<LineString>()};
    for (auto h = f->holes_begin(); h != f->holes_end(); h++) {
      auto he = *h;
      do {
        ring->addPoint(he->source()->point());
      } while (++he != *h);
    }

    if (ring->numPoints() > 3) {
      ring->addPoint(ring->startPoint());
      if (f->is_unbounded()) {
        poly->setExteriorRing(ring.release());
      } else if (allowHoles) {
        poly->addInteriorRing(ring.release());
      }
    }
  }

  std::unique_ptr<Geometry> result = std::move(poly);

  return result;
}

/** @endcond */

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

auto
optimalAlphaShapes(const Geometry &geometry, bool allowHoles,
                   size_t nbComponents) -> std::unique_ptr<Geometry>
{
  Alpha_shape_2 alphaShape;
  const double  optimalAlpha{
      computeAlpha(geometry, alphaShape, 10000, nbComponents)};
  if (optimalAlpha < 0) {
    return std::make_unique<GeometryCollection>();
  }

  alphaShape.set_alpha(optimalAlpha);

  return alphaToGeometry(alphaShape, allowHoles);
}

auto
alphaShapes(const Geometry &geometry, double alpha, bool allowHoles)
    -> std::unique_ptr<Geometry>
{
  using CGAL::object_cast;
  Alpha_shape_2 alphaShape;
  const double  optimalAlpha{computeAlpha(geometry, alphaShape, alpha)};
  if (optimalAlpha < 0) {
    return std::make_unique<GeometryCollection>();
  }

  return alphaToGeometry(alphaShape, allowHoles);
}

} // namespace SFCGAL::algorithm
