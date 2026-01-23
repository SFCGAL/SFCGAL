// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_PRISM_H_
#define SFCGAL_PRISM_H_

#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/export.h"
#include "SFCGAL/primitive3d/Primitive.h"

#include <array>

namespace SFCGAL {

/**
 * @class Prism
 * @brief Represents a triangular prism (wedge) in 3D space.
 *
 * A triangular prism is a polyhedron with 5 faces:
 * - 2 triangular faces (base and top)
 * - 3 quadrilateral faces (sides)
 *
 * @verbatim
 *
 *        topVertex2
 *           /\
 *          /  \
 *         /    \
 *   topVertex0--topVertex1
 *        |      |
 *        |      |        (extrusion direction)
 *        |      |
 *   baseVertex0--baseVertex1
 *         \    /
 *          \  /
 *           \/
 *       baseVertex2
 *
 * @endverbatim
 *
 * This primitive is particularly useful for 3D chamfer operations where
 * triangular prisms are subtracted from solids to create beveled edges.
 *
 * @note The prism is defined by two parallel triangular faces. The vertices
 *       must be provided in consistent order (both clockwise or both
 *       counter-clockwise when viewed from outside).
 *
 * @see chamfer3D algorithm which uses this primitive
 *
 * References:
 * - "Computational Geometry: Algorithms and Applications" by de Berg et al.
 * - ISO 19107:2019 - Geographic information - Spatial schema
 */
class SFCGAL_API Prism : public PrimitiveImpl<Prism, Primitive> {
public:
  /**
   * @brief Default constructor creating a unit triangular prism.
   *
   * Creates a prism with:
   * - Base triangle at z=0: (0,0,0), (1,0,0), (0.5, sqrt(3)/2, 0)
   * - Top triangle at z=1: same x,y coordinates with z=1
   */
  Prism();

  /**
   * @brief Constructs a Prism from 6 vertices.
   *
   * The vertices must define two parallel triangular faces with consistent
   * orientation. The base and top triangles should have the same winding order
   * when viewed from outside the prism.
   *
   * @param baseVertex0 First vertex of the base triangle
   * @param baseVertex1 Second vertex of the base triangle
   * @param baseVertex2 Third vertex of the base triangle
   * @param topVertex0 First vertex of the top triangle (corresponds to
   * baseVertex0)
   * @param topVertex1 Second vertex of the top triangle (corresponds to
   * baseVertex1)
   * @param topVertex2 Third vertex of the top triangle (corresponds to
   * baseVertex2)
   *
   * @throws SFCGAL::Exception if the triangles are degenerate (collinear
   * points)
   *
   * @note The correspondence between base and top vertices determines the
   *       lateral faces: face 0 connects baseVertex0-baseVertex1 to
   *       topVertex0-topVertex1, etc.
   */
  Prism(const Kernel::Point_3 &baseVertex0, const Kernel::Point_3 &baseVertex1,
        const Kernel::Point_3 &baseVertex2, const Kernel::Point_3 &topVertex0,
        const Kernel::Point_3 &topVertex1, const Kernel::Point_3 &topVertex2);

  /**
   * @brief Constructs a Prism by extruding a triangle along a vector.
   *
   * This constructor creates a prism by translating the base triangle
   * along the extrusion vector. This is the most common use case for
   * chamfer operations.
   *
   * @param baseTriangle The triangle forming the base of the prism
   * @param extrusionVector The vector defining the direction and length of
   * extrusion
   *
   * @throws SFCGAL::Exception if the base triangle is empty or degenerate
   * @throws SFCGAL::Exception if the extrusion vector has zero length
   *
   * @note The top triangle vertices are computed as:
   *       topVertex[i] = baseTriangle.vertex(i) + extrusionVector
   */
  Prism(const Triangle &baseTriangle, const Kernel::Vector_3 &extrusionVector);

  /**
   * @brief Copy constructor
   * @param other The prism to copy from
   */
  Prism(const Prism &other) = default;

  /**
   * @brief Destructor
   */
  ~Prism() override = default;

  /// @copydoc SFCGAL::Primitive::primitiveType
  [[nodiscard]] auto
  primitiveType() const -> std::string override;

  /// @copydoc SFCGAL::Primitive::primitiveTypeId
  [[nodiscard]] auto
  primitiveTypeId() const -> PrimitiveType override;

  /**
   * @brief Gets a base vertex by index.
   * @param index Vertex index (0, 1, or 2)
   * @return The base vertex at the specified index
   * @throws std::out_of_range if index >= 3
   */
  [[nodiscard]] auto
  baseVertex(size_t index) const -> const Kernel::Point_3 &;

  /**
   * @brief Gets a top vertex by index.
   * @param index Vertex index (0, 1, or 2)
   * @return The top vertex at the specified index
   * @throws std::out_of_range if index >= 3
   */
  [[nodiscard]] auto
  topVertex(size_t index) const -> const Kernel::Point_3 &;

  /**
   * @brief Gets the base triangle.
   * @return A Triangle object representing the base face
   */
  [[nodiscard]] auto
  baseTriangle() const -> Triangle;

  /**
   * @brief Gets the top triangle.
   * @return A Triangle object representing the top face
   */
  [[nodiscard]] auto
  topTriangle() const -> Triangle;

  /**
   * @brief Computes the extrusion vector from base to top.
   *
   * This is computed as topVertex(0) - baseVertex(0).
   *
   * @return The vector from base vertex 0 to top vertex 0
   */
  [[nodiscard]] auto
  extrusionVector() const -> Kernel::Vector_3;

  /**
   * @brief Generates a PolyhedralSurface representation of the prism.
   *
   * The generated surface consists of 5 polygons:
   * - 1 base triangle (with normal pointing outward, typically downward)
   * - 1 top triangle (with normal pointing outward, typically upward)
   * - 3 quadrilateral side faces (each triangulated into 2 triangles
   *   for compatibility with some algorithms)
   *
   * @return A PolyhedralSurface representing the prism boundary
   *
   * @note The face orientations follow the convention that normals point
   *       outward from the solid interior.
   */
  auto
  generatePolyhedralSurface() const -> PolyhedralSurface override;

  /**
   * @brief Generates a closed Solid from the prism.
   *
   * Creates a Solid geometry with the prism's PolyhedralSurface as its
   * exterior shell. This is useful for boolean operations like
   * difference3D used in chamfering.
   *
   * @return A unique_ptr to a Solid representing the prism volume
   */
  [[nodiscard]] auto
  generateSolid() const -> std::unique_ptr<Solid>;

  /**
   * @copydoc SFCGAL::Primitive::volume
   *
   * The volume of a triangular prism is computed as:
   * V = (1/2) * |base_edge1 x base_edge2| * height
   *
   * where height is the projection of the extrusion vector onto
   * the base triangle normal.
   *
   * @note For non-right prisms (oblique extrusion), the formula
   *       accounts for the oblique height.
   */
  [[nodiscard]] auto
  volume(bool withDiscretization = false) const -> double override;

  /**
   * @copydoc SFCGAL::Primitive::area3D
   *
   * The surface area is the sum of:
   * - 2 * (area of base triangle)
   * - 3 * (area of each lateral quadrilateral face)
   */
  [[nodiscard]] auto
  area3D(bool withDiscretization = false) const -> double override;

protected:
  /**
   * @brief Validates the prism parameters.
   *
   * Checks that:
   * - Base and top triangles are not degenerate (non-zero area)
   * - The prism has non-zero volume
   *
   * @param tempParameters Parameter map to validate
   * @throws SFCGAL::Exception if validation fails
   */
  void
  validateParameters(std::unordered_map<std::string, PrimitiveParameter> const
                         &tempParameters) const override;

private:
  /// Base triangle vertices
  std::array<Kernel::Point_3, 3> _baseVertices;

  /// Top triangle vertices
  std::array<Kernel::Point_3, 3> _topVertices;

  /**
   * @brief Helper to check if three points are collinear (degenerate triangle).
   * @param point0 First point
   * @param point1 Second point
   * @param point2 Third point
   * @return true if points are collinear, false otherwise
   */
  [[nodiscard]] static auto
  areCollinear(const Kernel::Point_3 &point0, const Kernel::Point_3 &point1,
               const Kernel::Point_3 &point2) -> bool;
};

} // namespace SFCGAL

#endif // SFCGAL_PRISM_H_
