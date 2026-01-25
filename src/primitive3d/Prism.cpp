// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Prism.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"

#include <cmath>

namespace SFCGAL {

// ----------------------------------------------------------------------------
// Default constructor: creates a unit equilateral triangular prism
// ----------------------------------------------------------------------------
Prism::Prism()
{
  // Create an equilateral triangle base at z=0
  // Height of equilateral triangle with side 1: h = sqrt(3)/2 ≈ 0.866
  const double triangleHeight = std::sqrt(3.0) / 2.0;

  _baseVertices[0] = Kernel::Point_3(0.0, 0.0, 0.0);
  _baseVertices[1] = Kernel::Point_3(1.0, 0.0, 0.0);
  _baseVertices[2] = Kernel::Point_3(0.5, triangleHeight, 0.0);

  // Top triangle at z=1 (unit height extrusion)
  _topVertices[0] = Kernel::Point_3(0.0, 0.0, 1.0);
  _topVertices[1] = Kernel::Point_3(1.0, 0.0, 1.0);
  _topVertices[2] = Kernel::Point_3(0.5, triangleHeight, 1.0);

  // Note: Default prism has no parameters in the m_parameters map
  // as the vertices are stored directly in member arrays.
  // This is different from Box which stores extents as parameters.
}

// ----------------------------------------------------------------------------
// Constructor from 6 explicit vertices
// ----------------------------------------------------------------------------
Prism::Prism(const Kernel::Point_3 &baseVertex0,
             const Kernel::Point_3 &baseVertex1,
             const Kernel::Point_3 &baseVertex2,
             const Kernel::Point_3 &topVertex0,
             const Kernel::Point_3 &topVertex1,
             const Kernel::Point_3 &topVertex2)
{
  _baseVertices[0] = baseVertex0;
  _baseVertices[1] = baseVertex1;
  _baseVertices[2] = baseVertex2;

  _topVertices[0] = topVertex0;
  _topVertices[1] = topVertex1;
  _topVertices[2] = topVertex2;

  // Validate that we have non-degenerate triangles
  if (areCollinear(baseVertex0, baseVertex1, baseVertex2)) {
    BOOST_THROW_EXCEPTION(
        Exception("Prism base vertices are collinear (degenerate triangle)."));
  }

  if (areCollinear(topVertex0, topVertex1, topVertex2)) {
    BOOST_THROW_EXCEPTION(
        Exception("Prism top vertices are collinear (degenerate triangle)."));
  }
}

// ----------------------------------------------------------------------------
// Constructor by extruding a triangle along a vector
// ----------------------------------------------------------------------------
Prism::Prism(const Triangle &baseTriangle,
             const Kernel::Vector_3 &extrusionVector)
{
  // Validate inputs
  if (baseTriangle.isEmpty()) {
    BOOST_THROW_EXCEPTION(
        Exception("Cannot create Prism from an empty triangle."));
  }

  // Check that extrusion vector is not zero
  if (extrusionVector.squared_length() == 0) {
    BOOST_THROW_EXCEPTION(
        Exception("Prism extrusion vector cannot have zero length."));
  }

  // Extract base vertices from the triangle
  _baseVertices[0] = baseTriangle.vertex(0).toPoint_3();
  _baseVertices[1] = baseTriangle.vertex(1).toPoint_3();
  _baseVertices[2] = baseTriangle.vertex(2).toPoint_3();

  // Validate that base triangle is not degenerate
  if (areCollinear(_baseVertices[0], _baseVertices[1], _baseVertices[2])) {
    BOOST_THROW_EXCEPTION(
        Exception("Prism base triangle is degenerate (collinear vertices)."));
  }

  // Compute top vertices by translating along extrusion vector
  _topVertices[0] = _baseVertices[0] + extrusionVector;
  _topVertices[1] = _baseVertices[1] + extrusionVector;
  _topVertices[2] = _baseVertices[2] + extrusionVector;
}

// ----------------------------------------------------------------------------
// Primitive type information
// ----------------------------------------------------------------------------
auto
Prism::primitiveType() const -> std::string
{
  return "Prism";
}

auto
Prism::primitiveTypeId() const -> PrimitiveType
{
  return PrimitiveType::TYPE_PRISM;
}

// ----------------------------------------------------------------------------
// Vertex accessors
// ----------------------------------------------------------------------------
auto
Prism::baseVertex(size_t index) const -> const Kernel::Point_3 &
{
  if (index >= 3) {
    BOOST_THROW_EXCEPTION(Exception("Prism vertex index out of range (0-2)."));
  }
  return _baseVertices[index];
}

auto
Prism::topVertex(size_t index) const -> const Kernel::Point_3 &
{
  if (index >= 3) {
    BOOST_THROW_EXCEPTION(Exception("Prism vertex index out of range (0-2)."));
  }
  return _topVertices[index];
}

auto
Prism::baseTriangle() const -> Triangle
{
  return Triangle(Point(_baseVertices[0]), Point(_baseVertices[1]),
                  Point(_baseVertices[2]));
}

auto
Prism::topTriangle() const -> Triangle
{
  return Triangle(Point(_topVertices[0]), Point(_topVertices[1]),
                  Point(_topVertices[2]));
}

auto
Prism::extrusionVector() const -> Kernel::Vector_3
{
  return _topVertices[0] - _baseVertices[0];
}

// ----------------------------------------------------------------------------
// Generate PolyhedralSurface representation
// ----------------------------------------------------------------------------
auto
Prism::generatePolyhedralSurface() const -> PolyhedralSurface
{
  // Check cache
  if (m_polyhedral_surface) {
    return *m_polyhedral_surface;
  }

  m_polyhedral_surface = PolyhedralSurface();

  // Convert vertices to SFCGAL Points for polygon construction
  std::array<Point, 3> basePoints = {Point(_baseVertices[0]),
                                     Point(_baseVertices[1]),
                                     Point(_baseVertices[2])};

  std::array<Point, 3> topPoints = {Point(_topVertices[0]),
                                    Point(_topVertices[1]),
                                    Point(_topVertices[2])};

  // -------------------------------------------------------------------------
  // Face construction follows outward-pointing normal convention:
  // - Base face: vertices in reverse order (normal points "down"/outward)
  // - Top face: vertices in original order (normal points "up"/outward)
  // - Side faces: form quadrilaterals connecting base to top edges
  // -------------------------------------------------------------------------

  // Base triangle face (normal pointing outward from base)
  // Reverse winding to get outward normal (assuming base is "bottom")
  m_polyhedral_surface->addPatch(
      LineString({basePoints[0], basePoints[2], basePoints[1], basePoints[0]}));

  // Top triangle face (normal pointing outward from top)
  m_polyhedral_surface->addPatch(
      LineString({topPoints[0], topPoints[1], topPoints[2], topPoints[0]}));

  // -------------------------------------------------------------------------
  // Three lateral quadrilateral faces
  // Each connects an edge of the base to the corresponding edge of the top.
  // The faces are defined to have outward-pointing normals.
  //
  // Side face i connects:
  //   baseVertices[i] -> baseVertices[(i+1)%3] -> topVertices[(i+1)%3] ->
  //   topVertices[i]
  // -------------------------------------------------------------------------

  // Side face 0: edge 0-1
  m_polyhedral_surface->addPatch(
      LineString({basePoints[0], basePoints[1], topPoints[1], topPoints[0],
                  basePoints[0]}));

  // Side face 1: edge 1-2
  m_polyhedral_surface->addPatch(
      LineString({basePoints[1], basePoints[2], topPoints[2], topPoints[1],
                  basePoints[1]}));

  // Side face 2: edge 2-0
  m_polyhedral_surface->addPatch(
      LineString({basePoints[2], basePoints[0], topPoints[0], topPoints[2],
                  basePoints[2]}));

  return *m_polyhedral_surface;
}

// ----------------------------------------------------------------------------
// Generate Solid from the prism
// ----------------------------------------------------------------------------
auto
Prism::generateSolid() const -> std::unique_ptr<Solid>
{
  PolyhedralSurface surface = generatePolyhedralSurface();
  return std::make_unique<Solid>(surface);
}

// ----------------------------------------------------------------------------
// Volume calculation
// ----------------------------------------------------------------------------
auto
Prism::volume(bool /*withDiscretization*/) const -> double
{
  // Volume of a triangular prism:
  // V = base_area * height
  //
  // For a general (possibly oblique) prism:
  // V = base_area * |extrusion_vector · base_normal| / |base_normal|
  //   = base_area * (extrusion component perpendicular to base)
  //
  // Since we want the actual volume (not signed), we take absolute value.

  // Compute base triangle area using cross product
  // Area = 0.5 * |edge1 x edge2|
  Kernel::Vector_3 edge1 = _baseVertices[1] - _baseVertices[0];
  Kernel::Vector_3 edge2 = _baseVertices[2] - _baseVertices[0];
  Kernel::Vector_3 crossProduct = CGAL::cross_product(edge1, edge2);

  // Base area = 0.5 * |cross_product|
  double crossLength =
      std::sqrt(CGAL::to_double(crossProduct.squared_length()));
  double baseArea = 0.5 * crossLength;

  // For the height, we need the component of extrusion vector along the normal
  // Normal direction = cross_product (not normalized)
  // Height = |extrusion · normal| / |normal|
  Kernel::Vector_3 extrusion = extrusionVector();

  // Dot product of extrusion with normal
  double dotProduct = CGAL::to_double(extrusion * crossProduct);

  // Height = |dot_product| / |cross_product|
  double height = std::abs(dotProduct) / crossLength;

  return baseArea * height;
}

// ----------------------------------------------------------------------------
// Surface area calculation
// ----------------------------------------------------------------------------
auto
Prism::area3D(bool /*withDiscretization*/) const -> double
{
  // Total surface area = 2 * (base triangle area) + sum of 3 lateral face areas

  // Helper lambda to compute triangle area from 3 points
  auto triangleArea = [](const Kernel::Point_3 &point0,
                         const Kernel::Point_3 &point1,
                         const Kernel::Point_3 &point2) -> double {
    Kernel::Vector_3 edge1 = point1 - point0;
    Kernel::Vector_3 edge2 = point2 - point0;
    Kernel::Vector_3 cross = CGAL::cross_product(edge1, edge2);
    return 0.5 * std::sqrt(CGAL::to_double(cross.squared_length()));
  };

  // Helper lambda to compute quadrilateral area (as sum of 2 triangles)
  auto quadArea = [&triangleArea](const Kernel::Point_3 &point0,
                                  const Kernel::Point_3 &point1,
                                  const Kernel::Point_3 &point2,
                                  const Kernel::Point_3 &point3) -> double {
    // Split quad into 2 triangles: (p0, p1, p2) and (p0, p2, p3)
    return triangleArea(point0, point1, point2) +
           triangleArea(point0, point2, point3);
  };

  double totalArea = 0.0;

  // Base and top triangles (same area for parallel faces)
  double baseArea =
      triangleArea(_baseVertices[0], _baseVertices[1], _baseVertices[2]);
  totalArea += 2.0 * baseArea;

  // Three lateral faces (quadrilaterals)
  // Face 0: base[0], base[1], top[1], top[0]
  totalArea += quadArea(_baseVertices[0], _baseVertices[1], _topVertices[1],
                        _topVertices[0]);

  // Face 1: base[1], base[2], top[2], top[1]
  totalArea += quadArea(_baseVertices[1], _baseVertices[2], _topVertices[2],
                        _topVertices[1]);

  // Face 2: base[2], base[0], top[0], top[2]
  totalArea += quadArea(_baseVertices[2], _baseVertices[0], _topVertices[0],
                        _topVertices[2]);

  return totalArea;
}

// ----------------------------------------------------------------------------
// Parameter validation
// ----------------------------------------------------------------------------
void
Prism::validateParameters(
    std::unordered_map<std::string, PrimitiveParameter> const
        & /*tempParameters*/) const
{
  // The Prism class stores vertices directly rather than using parameters,
  // so validation is done in the constructors.
  // This method is kept for interface compatibility with the Primitive base
  // class.
}

// ----------------------------------------------------------------------------
// Helper: check if three points are collinear
// ----------------------------------------------------------------------------
auto
Prism::areCollinear(const Kernel::Point_3 &point0,
                    const Kernel::Point_3 &point1,
                    const Kernel::Point_3 &point2) -> bool
{
  // Three points are collinear if the cross product of two edges is zero
  Kernel::Vector_3 edge1 = point1 - point0;
  Kernel::Vector_3 edge2 = point2 - point0;
  Kernel::Vector_3 cross = CGAL::cross_product(edge1, edge2);

  // Check if cross product has zero length (using exact arithmetic)
  return cross.squared_length() == 0;
}

} // namespace SFCGAL
