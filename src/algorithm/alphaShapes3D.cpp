// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/alphaShapes3D.h"
#include "SFCGAL/detail/GetPointsVisitor.h"

#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/** @cond */

using Vb3  = CGAL::Alpha_shape_vertex_base_3<Kernel>;
using Fb3  = CGAL::Alpha_shape_cell_base_3<Kernel>;
using Tds3 = CGAL::Triangulation_data_structure_3<Vb3, Fb3>;
using Delaunay3 =
    CGAL::Delaunay_triangulation_3<Kernel, Tds3, CGAL::Fast_location>;
using Alpha_shape_3 = CGAL::Alpha_shape_3<Delaunay3>;
using Point_3       = Kernel::Point_3;

//! Reads the geometry and compute alphaShape
auto
buildAlphaShape(const Geometry &geom) -> std::unique_ptr<Alpha_shape_3>
{
  if (geom.isEmpty()) {
    return nullptr;
  }

  // Collect points from geometry
  SFCGAL::detail::GetPointsVisitor getPointVisitor;
  const_cast<Geometry &>(geom).accept(getPointVisitor);

  // Need at least 4 points for 3D alpha shapes
  if (getPointVisitor.points.size() < 4) {
    return nullptr;
  }

  // Create points vector
  std::vector<Point_3> points;
  points.reserve(getPointVisitor.points.size());
  for (const auto &point : getPointVisitor.points) {
    points.push_back(point->toPoint_3());
  }

  // Create and configure alpha shape
  auto alphaShape =
      std::make_unique<Alpha_shape_3>(points.begin(), points.end());
  alphaShape->set_mode(Alpha_shape_3::REGULARIZED);
  // find_alpha_solid()
  // Set optimal alpha value for one solid component
  auto optimalAlphaValue = alphaShape->find_optimal_alpha(1);
  alphaShape->set_alpha(*optimalAlphaValue);

  return alphaShape;
}

//! Converts alphaShape to a PolyhedralSurface
// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
alphaShape3DToPolyhedralSurface(const Alpha_shape_3 &alphaShape,
                                bool                 allowCavities)
    -> std::unique_ptr<PolyhedralSurface>
{

  // Iterate through all facets of the alpha shape
  std::vector<Point_3>                                points;
  std::vector<std::array<std::size_t, 3>>             polygons;
  std::map<Alpha_shape_3::Vertex_handle, std::size_t> vertexIds;

  // helper to avoid adding duplicate points
  auto getVertexId =
      [&points,
       &vertexIds](Alpha_shape_3::Vertex_handle vertexHandle) -> unsigned long {
    auto insert_res = vertexIds.emplace(vertexHandle, points.size());
    if (insert_res.second) {
      points.push_back(vertexHandle->point());
    }
    return insert_res.first->second;
  };

  std::set<Alpha_shape_3::Cell_handle> nonBoundedComponent;

  if (!allowCavities) {
    std::vector<Alpha_shape_3::Cell_handle> stack;
    stack.push_back(alphaShape.infinite_cell());
    while (!stack.empty()) {
      auto cell = stack.back();
      stack.pop_back();
      if (!nonBoundedComponent.insert(cell).second) {
        continue;
      }

      for (int i = 0; i < 4; ++i) {
        auto cell_n = cell->neighbor(i);
        if (alphaShape.classify(cell_n) == Alpha_shape_3::EXTERIOR &&
            nonBoundedComponent.count(cell_n) == 0) {
          stack.push_back(cell_n);
        }
      }
    }
  }

  for (auto facetIterator = alphaShape.alpha_shape_facets_begin();
       facetIterator != alphaShape.alpha_shape_facets_end(); ++facetIterator) {

    if (alphaShape.classify(*facetIterator) == Alpha_shape_3::REGULAR) {
      auto facet = *facetIterator;

      // always look at the facet from the outside
      if (alphaShape.classify(facet.first) != Alpha_shape_3::EXTERIOR) {
        facet = alphaShape.mirror_facet(facet);
      }

      // Get facet vertices
      auto cell     = facet.first;
      int  facetIdx = facet.second;

      if (!allowCavities && nonBoundedComponent.count(cell) == 0) {
        continue;
      }

      polygons.push_back(
          CGAL::make_array(getVertexId(cell->vertex((facetIdx + 1) % 4)),
                           getVertexId(cell->vertex((facetIdx + 2) % 4)),
                           getVertexId(cell->vertex((facetIdx + 3) % 4))));

      // fix orientation of triangles
      if (facetIdx % 2 == 0) {
        std::swap(polygons.back()[0], polygons.back()[1]);
      }
    }
  }

  // we need to orient as solid does not imply manifold
  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
  CGAL::Surface_mesh<Point_3> surfaceMesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons,
                                                              surfaceMesh);

  return std::make_unique<PolyhedralSurface>(surfaceMesh);
}
// NOLINTEND(readability-function-cognitive-complexity)

/** @endcond */

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

auto
alphaShapes3D(const Geometry &geometry, bool allowCavities)
    -> std::unique_ptr<PolyhedralSurface>
{
  std::unique_ptr<Alpha_shape_3> alphaShape = buildAlphaShape(geometry);
  if (!alphaShape) {
    return std::make_unique<PolyhedralSurface>();
  }
  return alphaShape3DToPolyhedralSurface(*alphaShape, allowCavities);
}

} // namespace SFCGAL::algorithm
