// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/surfaceSimplification.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Surface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/isValid.h"

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#ifdef SFCGAL_WITH_EIGEN
  #include <CGAL/Cartesian_converter.h>
  #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
  #include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_plane_policies.h>
  #include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
  #include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#endif // SFCGAL_WITH_EIGEN

#include <map>
#include <stdexcept>

namespace SFCGAL::algorithm {

namespace SMS = CGAL::Surface_mesh_simplification;

namespace detail {

#ifdef SFCGAL_WITH_EIGEN
// Inexact kernel types for advanced strategies
using InexactKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using InexactMesh   = CGAL::Surface_mesh<InexactKernel::Point_3>;
using EK_to_IK      = CGAL::Cartesian_converter<Kernel, InexactKernel>;
using IK_to_EK      = CGAL::Cartesian_converter<InexactKernel, Kernel>;
#endif // SFCGAL_WITH_EIGEN

/**
 * @brief Apply edge collapse simplification using Edge Length cost and Midpoint
 * placement
 */
template <typename StopPredicate>
auto
simplifyMesh(Surface_mesh_3 &mesh, const StopPredicate &stop) -> size_t
{
  using Cost      = SMS::Edge_length_cost<Surface_mesh_3>;
  using Placement = SMS::Midpoint_placement<Surface_mesh_3>;

  return SMS::edge_collapse(
      mesh, stop,
      CGAL::parameters::get_cost(Cost()).get_placement(Placement()));
}

#ifdef SFCGAL_WITH_EIGEN
/**
 * @brief Apply edge collapse simplification using Garland-Heckbert strategy on
 * inexact mesh
 */
template <typename StopPredicate>
auto
simplifyMeshGarlandHeckbert(InexactMesh &mesh, const StopPredicate &stop)
    -> size_t
{
  using Policies =
      SMS::GarlandHeckbert_plane_policies<InexactMesh, InexactKernel>;
  Policies policies(mesh);

  return SMS::edge_collapse(mesh, stop,
                            CGAL::parameters::get_cost(policies.get_cost())
                                .get_placement(policies.get_placement()));
}

/**
 * @brief Apply edge collapse simplification using Lindstrom-Turk strategy on
 * inexact mesh
 */
template <typename StopPredicate>
auto
simplifyMeshLindstromTurk(InexactMesh &mesh, const StopPredicate &stop)
    -> size_t
{
  using Cost      = SMS::LindstromTurk_cost<InexactMesh>;
  using Placement = SMS::LindstromTurk_placement<InexactMesh>;

  return SMS::edge_collapse(
      mesh, stop,
      CGAL::parameters::get_cost(Cost()).get_placement(Placement()));
}
#endif // SFCGAL_WITH_EIGEN

/**
 * @brief Call a simplification function on a mesh (exact or inexact) depending
 * on a stop predicate.
 *
 * The function handles edge count as well as edge ratio stop predicates.
 */
template <typename Mesh, typename SimplifyFn>
auto
simplifyWithStopPredicate(Mesh &mesh, const SimplificationStopPredicate &stop,
                          SimplifyFn &&simplifyFn)
{
  size_t result;
  if (stop.type == SimplificationStopPredicate::PredicateType::EDGE_COUNT) {

    auto const initialEdgeCount = static_cast<size_t>(mesh.number_of_edges());
    auto const targetEdgeCount  = static_cast<size_t>(stop.value);

    if (targetEdgeCount >= initialEdgeCount) {
      return static_cast<size_t>(0);
    }

    size_t const edgesToRemove = initialEdgeCount - targetEdgeCount;
    SMS::Edge_count_stop_predicate<Mesh> predicate(edgesToRemove);
    result = simplifyFn(mesh, predicate);

  } else {

    double const ratio = stop.value;
    if (ratio <= 0.0 || ratio >= 1.0) {
      throw std::invalid_argument(
          "Edge count ratio must be in the range (0.0, 1.0)");
    }

    SMS::Edge_count_ratio_stop_predicate<Mesh> predicate(ratio);
    result = simplifyFn(mesh, predicate);
  }
  return result;
}

#ifdef SFCGAL_WITH_EIGEN
/**
 * @brief Simplify inexact meshes with Eigen strategies.
 */
auto
simplifyInexactMesh(Surface_mesh_3                    &mesh,
                    const SimplificationStopPredicate &stop,
                    SimplificationStrategy             strategy)
{
  constexpr size_t MAX_VERTICES = 500000;
  constexpr size_t MAX_FACES    = 1000000;

  if (mesh.number_of_vertices() > MAX_VERTICES ||
      mesh.number_of_faces() > MAX_FACES) {
    throw std::invalid_argument("Mesh too large for kernel conversion - use "
                                "EDGE_LENGTH strategy instead");
  }

  EK_to_IK    toInexact;
  IK_to_EK    toExact;
  InexactMesh inexactMesh;

  std::map<Surface_mesh_3::Vertex_index, InexactMesh::Vertex_index> vmap;

  for (auto vh : mesh.vertices()) {
    vmap[vh] = inexactMesh.add_vertex(toInexact(mesh.point(vh)));
  }

  for (auto fh : mesh.faces()) {
    std::vector<InexactMesh::Vertex_index> verts;
    for (auto vh : mesh.vertices_around_face(mesh.halfedge(fh))) {
      verts.push_back(vmap[vh]);
    }
    inexactMesh.add_face(verts);
  }

  auto simplifyFn = [&](auto &m, auto &predicate) {
    if (strategy == SimplificationStrategy::GARLAND_HECKBERT) {
      return simplifyMeshGarlandHeckbert(m, predicate);
    }
    return simplifyMeshLindstromTurk(m, predicate);
  };

  size_t result = simplifyWithStopPredicate(inexactMesh, stop, simplifyFn);

  mesh.clear();
  std::map<InexactMesh::Vertex_index, Surface_mesh_3::Vertex_index> backMap;

  for (auto vh : inexactMesh.vertices()) {
    backMap[vh] = mesh.add_vertex(toExact(inexactMesh.point(vh)));
  }

  for (auto fh : inexactMesh.faces()) {
    std::vector<Surface_mesh_3::Vertex_index> verts;
    for (auto vh : inexactMesh.vertices_around_face(inexactMesh.halfedge(fh))) {
      verts.push_back(backMap[vh]);
    }
    mesh.add_face(verts);
  }

  return result;
}
#endif

/**
 * @brief Entry point for simplifying surface meshes
 */
auto
simplifySurfaceMesh(Surface_mesh_3                    &mesh,
                    const SimplificationStopPredicate &stop,
                    SimplificationStrategy             strategy)
{
#ifdef SFCGAL_WITH_EIGEN
  if (strategy == SimplificationStrategy::GARLAND_HECKBERT ||
      strategy == SimplificationStrategy::LINDSTROM_TURK) {
    return simplifyInexactMesh(mesh, stop, strategy);
  }
#else
  (void)strategy;
#endif

  return simplifyWithStopPredicate(
      mesh, stop, [](Surface_mesh_3 &m, const auto &predicate) {
        return simplifyMesh(m, predicate);
      });
  ;
}

/**
 * @brief Simplify a TriangulatedSurface
 */
auto
simplifyTriangulatedSurface(const TriangulatedSurface         &surface,
                            const SimplificationStopPredicate &stopPredicate,
                            SimplificationStrategy             strategy)
    -> std::unique_ptr<TriangulatedSurface>
{
  if (surface.isEmpty()) {
    return std::make_unique<TriangulatedSurface>(surface);
  }

  // Convert to Surface_mesh
  Surface_mesh_3 mesh = surface.toSurfaceMesh();

  // Simplify
  simplifySurfaceMesh(mesh, stopPredicate, strategy);

  // Convert back to PolyhedralSurface first, then to TriangulatedSurface
  PolyhedralSurface polyhedral(mesh);
  return std::make_unique<TriangulatedSurface>(
      polyhedral.toTriangulatedSurface());
}

/**
 * @brief Simplify a PolyhedralSurface
 */
auto
simplifyPolyhedralSurface(const PolyhedralSurface           &surface,
                          const SimplificationStopPredicate &stopPredicate,
                          SimplificationStrategy             strategy)
    -> std::unique_ptr<PolyhedralSurface>
{
  if (surface.isEmpty()) {
    return std::make_unique<PolyhedralSurface>(surface);
  }

  // Convert to Surface_mesh
  Surface_mesh_3 mesh = surface.toSurfaceMesh();

  // Simplify
  simplifySurfaceMesh(mesh, stopPredicate, strategy);

  // Convert back to PolyhedralSurface
  return std::make_unique<PolyhedralSurface>(mesh);
}

/**
 * @brief Simplify a Solid (only exterior shell)
 */
auto
simplifySolid(const Solid                       &solid,
              const SimplificationStopPredicate &stopPredicate,
              SimplificationStrategy strategy) -> std::unique_ptr<Solid>
{
  if (solid.isEmpty()) {
    return std::make_unique<Solid>(solid);
  }

  // Simplify the exterior shell
  auto simplifiedExterior =
      simplifyPolyhedralSurface(solid.exteriorShell(), stopPredicate, strategy);

  // Create new solid with simplified exterior shell
  auto result = std::make_unique<Solid>(*simplifiedExterior);

  // Copy interior shells (not simplified)
  for (size_t i = 0; i < solid.numInteriorShells(); ++i) {
    result->addInteriorShell(solid.interiorShellN(i).clone());
  }

  return result;
}

/**
 * @brief Simplify a MultiSolid
 */
auto
simplifyMultiSolid(const MultiSolid                  &multiSolid,
                   const SimplificationStopPredicate &stopPredicate,
                   SimplificationStrategy             strategy)
    -> std::unique_ptr<MultiSolid>
{
  auto result = std::make_unique<MultiSolid>();

  for (size_t i = 0; i < multiSolid.numGeometries(); ++i) {
    const auto &solid = multiSolid.geometryN(i).as<Solid>();
    result->addGeometry(simplifySolid(solid, stopPredicate, strategy));
  }

  return result;
}

} // namespace detail

/**
 * @brief Simplify a surface mesh using CGAL edge collapse algorithm
 * @param geometry The input geometry to simplify (must be a surface or solid)
 * @param stopPredicate When to stop the simplification process
 * @param strategy The cost and placement strategy to use
 * @return A simplified copy of the input geometry
 * @pre The input geometry must be valid and non-empty
 * @pre For EDGE_COUNT_RATIO, the ratio must be in the range (0.0, 1.0)
 * @pre The geometry must be 3-dimensional
 */
auto
surfaceSimplification(const Geometry                    &geometry,
                      const SimplificationStopPredicate &stopPredicate,
                      SimplificationStrategy             strategy)
    -> std::unique_ptr<Geometry>
{
  // Validate geometry
  switch (geometry.geometryTypeId()) {
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(geometry);
    break;
  default:
    break;
  }

  std::unique_ptr<Geometry> result(surfaceSimplification(
      geometry, stopPredicate, strategy, NoValidityCheck()));
  propagateValidityFlag(*result, true);
  return result;
}

auto
surfaceSimplification(const Geometry                    &geometry,
                      const SimplificationStopPredicate &stopPredicate,
                      SimplificationStrategy             strategy,
                      NoValidityCheck /*unused*/) -> std::unique_ptr<Geometry>
{
  if (geometry.isEmpty()) {
    return geometry.clone();
  }

  // Validate ratio for EDGE_COUNT_RATIO
  if (stopPredicate.type ==
      SimplificationStopPredicate::PredicateType::EDGE_COUNT_RATIO) {
    if (stopPredicate.value <= 0.0 || stopPredicate.value >= 1.0) {
      throw std::invalid_argument(
          "Edge count ratio must be in the range (0.0, 1.0)");
    }
  }

  switch (geometry.geometryTypeId()) {
  case TYPE_TRIANGULATEDSURFACE:
    return detail::simplifyTriangulatedSurface(
        geometry.as<TriangulatedSurface>(), stopPredicate, strategy);

  case TYPE_POLYHEDRALSURFACE:
    return detail::simplifyPolyhedralSurface(geometry.as<PolyhedralSurface>(),
                                             stopPredicate, strategy);

  case TYPE_SOLID:
    return detail::simplifySolid(geometry.as<Solid>(), stopPredicate, strategy);

  case TYPE_MULTISOLID:
    return detail::simplifyMultiSolid(geometry.as<MultiSolid>(), stopPredicate,
                                      strategy);

  default:
    throw std::invalid_argument(
        "surfaceSimplification only supports TriangulatedSurface, "
        "PolyhedralSurface, Solid, and MultiSolid geometries");
  }
}

} // namespace SFCGAL::algorithm
