// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <fstream>
#include <iostream>

/**
 * face information (depth)
 */
struct FaceInfo2 {
  FaceInfo2() = default;
  int nesting_level;

  [[nodiscard]] auto
  in_domain() const -> bool
  {
    return nesting_level % 2 == 1;
  }
};

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using triangulation_vertex_base = CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel>;
using triangulation_face_base = CGAL::Triangulation_face_base_with_info_2<FaceInfo2, Kernel>;
using constrained_triangulation_face_base = CGAL::Constrained_triangulation_face_base_2<Kernel,
                                                    triangulation_face_base>;
using triangulation_data_structure = CGAL::Triangulation_data_structure_2<
    triangulation_vertex_base, constrained_triangulation_face_base>;

// typedef CGAL::Exact_predicates_tag Itag;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<
    Kernel, triangulation_data_structure, CGAL::Exact_predicates_tag>;

using triangulation_point = CDT::Point;
using Point_2 = CGAL::Point_2<Kernel>;
using Polygon_2 = CGAL::Polygon_2<Kernel>;

void
mark_domains(CDT &ct, CDT::Face_handle start, int index,
             std::list<CDT::Edge> &border)
{
  if (start->info().nesting_level != -1) {
    return;
  }
  std::list<CDT::Face_handle> queue;
  queue.push_back(start);

  while (!queue.empty()) {
    CDT::Face_handle faceHandle = queue.front();
    queue.pop_front();
    if (faceHandle->info().nesting_level == -1) {
      faceHandle->info().nesting_level = index;
      for (int i = 0; i < 3; i++) {
        CDT::Edge        edge(faceHandle, i);
        CDT::Face_handle neighbor = faceHandle->neighbor(i);
        if (neighbor->info().nesting_level == -1) {
          if (ct.is_constrained(edge)) {
            border.push_back(edge);
          }
          else {
            queue.push_back(neighbor);
          }
        }
      }
    }
  }
}

// explore set of facets connected with non constrained edges,
// and attribute to each such set a nesting level.
// We start from facets incident to the infinite vertex, with a nesting
// level of 0. Then we recursively consider the non-explored facets incident
// to constrained edges bounding the former set and increase the nesting level
// by 1. Facets in the domain are those with an odd nesting level.
void
mark_domains(CDT &cdt)
{
  for (CDT::All_faces_iterator it = cdt.all_faces_begin();
       it != cdt.all_faces_end(); ++it) {
    it->info().nesting_level = -1;
  }

  int                  index = 0;
  std::list<CDT::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), index++, border);
  while (!border.empty()) {
    CDT::Edge edge = border.front();
    border.pop_front();
    CDT::Face_handle neighbor = edge.first->neighbor(edge.second);
    if (neighbor->info().nesting_level == -1) {
      mark_domains(cdt, neighbor, edge.first->info().nesting_level + 1, border);
    }
  }
}

void
insert_polygon(CDT &cdt, const Polygon_2 &polygon)
{
  if (polygon.is_empty()) {
    return;
  }
  CDT::Vertex_handle v_prev =
      cdt.insert(*CGAL::cpp0x::prev(polygon.vertices_end()));
  for (auto vit = polygon.vertices_begin();
       vit != polygon.vertices_end(); ++vit) {
    CDT::Vertex_handle vertexHandle = cdt.insert(*vit);
    cdt.insert_constraint(vertexHandle, v_prev);
    v_prev = vertexHandle;
  }
}

auto
main() -> int
{
  // construct two non-intersecting nested polygons
  Polygon_2 polygon1;
  polygon1.push_back(Point_2(0.0, 0.0));
  polygon1.push_back(Point_2(2.0, 0.0));
  polygon1.push_back(Point_2(1.7, 1.0));
  polygon1.push_back(Point_2(2.0, 2.0));
  polygon1.push_back(Point_2(0.0, 2.0));
  Polygon_2 polygon2;
  polygon2.push_back(Point_2(0.5, 0.5));
  polygon2.push_back(Point_2(1.5, 0.5));
  polygon2.push_back(Point_2(1.5, 1.5));
  polygon2.push_back(Point_2(0.5, 1.5));

  // Insert the polyons into a constrained triangulation
  CDT cdt;
  insert_polygon(cdt, polygon1);
  insert_polygon(cdt, polygon2);

  // Extract point and provide the an index
  std::vector<triangulation_point> points;
  for (CDT::Vertex_iterator it = cdt.vertices_begin(); it != cdt.vertices_end();
       ++it) {
    it->info() = points.size();
    points.push_back(it->point());
  }

  // Mark facets that are inside the domain bounded by the polygon
  mark_domains(cdt);

  //
  int count = 0;
  for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
       fit != cdt.finite_faces_end(); ++fit) {
    if (fit->info().in_domain()) {
      ++count;
    }
  }

  /*
   * export
   */

  std::ofstream ofs("polygon_triangulation2.obj");
  if (!ofs.good()) {
    std::cout << "can't open file\n";
    return 1;
  }

  //-- print vertices
  ofs << "# " << points.size() << " vertices\n";
  for (auto point : points) {
    ofs << "v " << point << " 0.0\n";
  }

  //-- print faces
  ofs << "# " << cdt.number_of_faces() << " faces\n";
  // warning : Delaunay_triangulation_2::All_faces_iterator iterator over
  // infinite faces
  for (CDT::Finite_faces_iterator it = cdt.finite_faces_begin();
       it != cdt.finite_faces_end(); ++it) {
    // ignore holes
    if (!it->info().in_domain()) {
      continue;
    }
    size_t infoA = it->vertex(0)->info();
    size_t infoB = it->vertex(1)->info();
    size_t infoC = it->vertex(2)->info();

    assert(it->is_valid());
    // assert ( ia < cdt.number_of_vertices() || ib < tri.number_of_vertices()
    // || ic < tri.number_of_vertices() ) ;

    ofs << "f " << (infoA + 1) << " " << (infoB + 1) << " " << (infoC + 1) << "\n";
  }

  return 0;
}
