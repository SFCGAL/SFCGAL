// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <CGAL/number_utils.h>
#include <utility>

#include "SFCGAL/numeric.h"
#include "SFCGAL/primitive3d/Cylinder.h"

#include <CGAL/Polygon_mesh_processing/transform.h>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace SFCGAL {

Cylinder::Cylinder(const Kernel::FT &radius, const Kernel::FT &height,
                   unsigned int num_radial)
{
  m_parameters["radius"]     = radius;
  m_parameters["height"]     = height;
  m_parameters["num_radial"] = num_radial;

  Cylinder::validateParameters(m_parameters);
}

auto
Cylinder::primitiveType() const -> std::string
{
  return "Cylinder";
}

auto
Cylinder::primitiveTypeId() const -> PrimitiveType
{
  return PrimitiveType::TYPE_CYLINDER;
}

void
Cylinder::setRadius(const Kernel::FT &radius)
{
  validateAndSetParameter("radius", radius);
}

void
Cylinder::setHeight(const Kernel::FT &height)
{
  validateAndSetParameter("height", height);
}

void
Cylinder::setNumRadial(unsigned int num)
{
  validateAndSetParameter("num_radial", num);
}

void
Cylinder::validateParameters(
    std::unordered_map<std::string, PrimitiveParameter> const &tempParameters)
    const
{
  const double radius =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("radius")));
  const double height =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("height")));
  const unsigned int num_radial =
      std::get<unsigned int>(tempParameters.at("num_radial"));

  if (radius <= 0.) {
    BOOST_THROW_EXCEPTION(Exception("Cylinder radius cannot be negative."));
  } else if (height <= 0.) {
    BOOST_THROW_EXCEPTION(Exception("Cylinder height cannot be negative."));
  } else if (num_radial < 3) {
    BOOST_THROW_EXCEPTION(
        Exception("Cylinder requires at least 3 radial segments."));
  }
}

void
Cylinder::invalidateCache()
{
  Primitive::invalidateCache();
  m_polyhedron.reset();
  m_surface_mesh.reset();
  m_polyhedral_surface.reset();
}

auto
Cylinder::generatePolyhedron() const -> Polyhedron_3
{
  if (m_polyhedron) {
    return *m_polyhedron;
  }

  Polyhedron_3   poly;
  Surface_mesh_3 sm = generateSurfaceMesh();
  CGAL::copy_face_graph(sm, poly);
  m_polyhedron = poly;
  return poly;
}

auto
Cylinder::generateSurfaceMesh() const -> Surface_mesh_3
{
  if (m_surface_mesh) {
    return *m_surface_mesh;
  }

  Surface_mesh_3 mesh;

  std::vector<Surface_mesh_3::Vertex_index> base_vertices;
  std::vector<Surface_mesh_3::Vertex_index> top_vertices;

  // Create vertices for the base and top
  for (unsigned int i = 0; i < numRadial(); ++i) {
    const double angle = 2.0 * M_PI * i / numRadial();
    const double x     = CGAL::to_double(radius()) * std::cos(angle);
    const double y     = CGAL::to_double(radius()) * std::sin(angle);

    base_vertices.push_back(mesh.add_vertex(Point_3(x, y, 0)));
    top_vertices.push_back(mesh.add_vertex(Point_3(x, y, height())));
  }

  // Add side faces
  for (unsigned int i = 0; i < numRadial(); ++i) {
    unsigned int next = (i + 1) % numRadial();
    mesh.add_face(base_vertices[i], top_vertices[next], top_vertices[i]);
    mesh.add_face(base_vertices[i], base_vertices[next], top_vertices[next]);
  }

  // Add base and top faces
  Surface_mesh_3::Vertex_index origin = mesh.add_vertex(Point_3(0, 0, 0));
  Surface_mesh_3::Vertex_index top_origin =
      mesh.add_vertex(Point_3(0, 0, height()));

  for (unsigned int i = 0; i < numRadial(); ++i) {
    unsigned int next = (i + 1) % numRadial();
    mesh.add_face(origin, base_vertices[next], base_vertices[i]);
    mesh.add_face(top_origin, top_vertices[i], top_vertices[next]);
  }

  // handle affine transformation
  if (m_transform != Kernel::Aff_transformation_3(CGAL::IDENTITY)) {
    PMP::transform(m_transform, mesh);
  }

  m_surface_mesh = mesh;
  return mesh;
}

auto
Cylinder::generatePolyhedralSurface() const -> PolyhedralSurface
{
  if (m_polyhedral_surface) {
    return *m_polyhedral_surface;
  }

  // generateSurfaceMesh already handles the affine transformation
  m_polyhedral_surface = PolyhedralSurface(generateSurfaceMesh());
  return *m_polyhedral_surface;
}

auto
Cylinder::volume(bool /*withDiscretization*/) const -> double
{
  return CGAL::to_double(radius() * radius() * height() * CGAL_PI);
}

auto
Cylinder::area3D(bool /*withDiscretization*/) const -> double
{
  return CGAL::to_double(2 * radius() * radius() * CGAL_PI +
                         2 * radius() * height() * CGAL_PI);
}

} // namespace SFCGAL
