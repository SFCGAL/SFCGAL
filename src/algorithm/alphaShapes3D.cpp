// Copyright (c) 2024-2024, SFCGAL Contributors and Oslandia
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/alphaShapes3D.h"
#include "SFCGAL/detail/GetPointsVisitor.h"

#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>


namespace SFCGAL::algorithm {

using Vb3  = CGAL::Alpha_shape_vertex_base_3<Kernel>;
using Fb3  = CGAL::Alpha_shape_cell_base_3<Kernel>;
using Tds3 = CGAL::Triangulation_data_structure_3<Vb3, Fb3>;
using Delaunay3 =
    CGAL::Delaunay_triangulation_3<Kernel, Tds3, CGAL::Fast_location>;
using Alpha_shape_3 = CGAL::Alpha_shape_3<Delaunay3>;
using Point_3       = Kernel::Point_3;

//! Reads the geometry and compute alphaShape
auto
buildAlphaShape(const Geometry &geom)
    -> std::unique_ptr<Alpha_shape_3>
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
auto
alphaShape3D_to_polyhedralSurface(const Alpha_shape_3 &alphaShape, bool allow_cavities)
    -> std::unique_ptr<PolyhedralSurface>
{

  // Iterate through all facets of the alpha shape
  std::vector<Point_3>                    points;
  std::vector<std::array<std::size_t, 3>> polygons;
  std::map<Alpha_shape_3::Vertex_handle, std::size_t> vids;

  auto getvid = [&points, &vids](Alpha_shape_3::Vertex_handle vh)
  {
    auto insert_res = vids.emplace(vh, points.size());
    if (insert_res.second)
    {
      points.push_back(vh->point());
    }
    return insert_res.first->second;
  };

  std::set<Alpha_shape_3::Cell_handle> non_bounded_component;

  if (!allow_cavities)
  {
    std::vector<Alpha_shape_3::Cell_handle> stack;
    stack.push_back(alphaShape.infinite_cell());
    while(!stack.empty())
    {
      auto cell = stack.back();
      stack.pop_back();
      if (!non_bounded_component.insert(cell).second)
        continue;
      for (int i=0;i<4; ++i)
      {
        auto cell_n=cell->neighbor(i);
        if (alphaShape.classify(cell_n)==Alpha_shape_3::EXTERIOR && non_bounded_component.count(cell_n)==0)
          stack.push_back(cell_n);
      }
    }
  }

  for (auto facetIterator = alphaShape.alpha_shape_facets_begin();
       facetIterator != alphaShape.alpha_shape_facets_end(); ++facetIterator) {

    if (alphaShape.classify(*facetIterator) == Alpha_shape_3::REGULAR) {
      auto facet = *facetIterator;

      // always look at the facet from the outside
      if (alphaShape.classify(facet.first)!=Alpha_shape_3::EXTERIOR)
      {
        facet = alphaShape.mirror_facet(facet);
      }

      // Get facet vertices
      auto cell     = facet.first;
      int  facetIdx = facet.second;

      if (!allow_cavities && non_bounded_component.count(cell)==0)
        continue;

      polygons.push_back(CGAL::make_array(
        getvid(cell->vertex((facetIdx + 1) % 4)),
        getvid(cell->vertex((facetIdx + 2) % 4)),
        getvid(cell->vertex((facetIdx + 3) % 4)))
      );

      // fix orientation of triangles
      if (facetIdx%2==0)
        std::swap(polygons.back()[0], polygons.back()[1]);
    }
  }

  // we need to orient as solid does not imply manifold
  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
  CGAL::Surface_mesh<Point_3> surfaceMesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons,
                                                              surfaceMesh);

  return std::make_unique<PolyhedralSurface>(PolyhedralSurface(surfaceMesh));
}

auto
alphaShapes3D(const Geometry &geom, bool allow_cavities)
    -> std::unique_ptr<PolyhedralSurface>
{
  std::unique_ptr<Alpha_shape_3> alphaShape = buildAlphaShape(geom);
  if (!alphaShape) {
    return std::make_unique<PolyhedralSurface>();
  }
  return alphaShape3D_to_polyhedralSurface(*alphaShape, allow_cavities);
}

} // namespace SFCGAL::algorithm
