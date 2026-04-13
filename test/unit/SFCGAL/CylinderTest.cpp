// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Cylinder.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/io/wkt.h"
#include <boost/test/unit_test.hpp>
#include <cmath>

using namespace SFCGAL;

#include "../../test_config.h"

BOOST_AUTO_TEST_SUITE(CylinderTests)

BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Cylinder cyl;
  BOOST_CHECK_CLOSE(cyl.radius(), 1.0, 1e-6);
  BOOST_CHECK_CLOSE(cyl.height(), 1.0, 1e-6);
  BOOST_CHECK_EQUAL(cyl.numRadial(), 32);
  BOOST_CHECK_EQUAL(cyl.baseCenter(), Point_3(0, 0, 0));
  BOOST_CHECK_EQUAL(cyl.axis(), Vector_3(0, 0, 1));
}

BOOST_AUTO_TEST_CASE(testCustomConstructor)
{
  Point_3  base(1, 2, 3);
  Vector_3 axis(0, 1, 0);
  Cylinder cyl(base, axis, 2.0, 5.0, 16);
  BOOST_CHECK_CLOSE(cyl.radius(), 2.0, 1e-6);
  BOOST_CHECK_CLOSE(cyl.height(), 5.0, 1e-6);
  BOOST_CHECK_EQUAL(cyl.numRadial(), 16);
  BOOST_CHECK_EQUAL(cyl.baseCenter(), base);
  BOOST_CHECK_EQUAL(cyl.axis(), axis);
}

BOOST_AUTO_TEST_CASE(testCopy)
{
  Point_3  base(1, 2, 3);
  Vector_3 axis(0, 1, 0);
  Cylinder cyl(base, axis, 2.0, 5.0, 16);
  BOOST_CHECK_CLOSE(cyl.radius(), 2.0, 1e-6);
  BOOST_CHECK_CLOSE(cyl.height(), 5.0, 1e-6);
  BOOST_CHECK_EQUAL(cyl.numRadial(), 16);
  BOOST_CHECK_EQUAL(cyl.baseCenter(), base);
  BOOST_CHECK_EQUAL(cyl.axis(), axis);

  std::string expectedWkt(SFCGAL_TEST_DIRECTORY);
  expectedWkt += "/data/cylinder_expected_cpp1.wkt";
  std::ifstream efs(expectedWkt.c_str());
  BOOST_REQUIRE(efs.good());
  std::getline(efs, expectedWkt);

  // Create PolyhedralSurface from Polyhedron and check WKT output
  auto polyhedral_surface = cyl.generatePolyhedralSurface();
  BOOST_CHECK(algorithm::isValid(polyhedral_surface));
  BOOST_CHECK_EQUAL(polyhedral_surface.asText(1), expectedWkt);

  auto              surface_mesh = cyl.generateSurfaceMesh();
  PolyhedralSurface poly_surface(surface_mesh);
  BOOST_CHECK(algorithm::isValid(poly_surface));
  BOOST_CHECK_EQUAL(poly_surface.asText(1), expectedWkt);

  Cylinder cyl2;
  cyl2 = cyl;

  // Create PolyhedralSurface from Polyhedron and check WKT output
  auto polyhedral_surface2 = cyl2.generatePolyhedralSurface();
  BOOST_CHECK(algorithm::isValid(polyhedral_surface2));
  BOOST_CHECK_EQUAL(polyhedral_surface2.asText(1), expectedWkt);

  auto              surface_mesh2 = cyl2.generateSurfaceMesh();
  PolyhedralSurface poly_surface2(surface_mesh2);
  BOOST_CHECK(algorithm::isValid(poly_surface2));
  BOOST_CHECK_EQUAL(poly_surface2.asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testSetters)
{
  Cylinder cyl;
  cyl.setRadius(3.0);
  cyl.setHeight(4.0);
  cyl.setNumRadial(24);
  cyl.setBaseCenter(Point_3(1, 1, 1));
  cyl.setAxis(Vector_3(1, 1, 1));

  BOOST_CHECK_CLOSE(cyl.radius(), 3.0, 1e-6);
  BOOST_CHECK_CLOSE(cyl.height(), 4.0, 1e-6);
  BOOST_CHECK_EQUAL(cyl.numRadial(), 24);
  BOOST_CHECK_EQUAL(cyl.baseCenter(), Point_3(1, 1, 1));
  BOOST_CHECK_EQUAL(cyl.axis(), Vector_3(1, 1, 1));

  // radius and height cannot be negative
  BOOST_CHECK_THROW(cyl.setHeight(-1.2), Exception);
  BOOST_CHECK_THROW(cyl.setRadius(-3.2), Exception);
}

BOOST_AUTO_TEST_CASE(testGenerateSurfaceMesh)
{
  Cylinder cyl(Point_3(0, 0, 0), Vector_3(0, 0, 1), 1.0, 2.0, 4);
  auto     mesh = cyl.generateSurfaceMesh();

  BOOST_CHECK_EQUAL(mesh.number_of_vertices(), cyl.numRadial() * 2 + 2);
  BOOST_CHECK_EQUAL(mesh.number_of_edges(), cyl.numRadial() * 6);
  BOOST_CHECK_EQUAL(mesh.number_of_faces(), cyl.numRadial() * 4);
}

BOOST_AUTO_TEST_CASE(testVolume)
{
  Cylinder cyl(Point_3(0, 0, 0), Vector_3(0, 0, 1), 2.0, 5.0, 32);
  double   volume          = cyl.volume();
  double   expected_volume = M_PI * 2.0 * 2.0 * 5.0;
  BOOST_CHECK_CLOSE(volume, expected_volume, 0.01);
}

BOOST_AUTO_TEST_CASE(testSurfaceArea)
{
  Cylinder cyl(Point_3(0, 0, 0), Vector_3(0, 0, 1), 2.0, 5.0, 32);
  double   area          = cyl.area3D();
  double   expected_area = 2 * M_PI * 2.0 * 2.0 + 2 * M_PI * 2.0 * 5.0;
  BOOST_CHECK_CLOSE(area, expected_area, 0.01);
}

BOOST_AUTO_TEST_CASE(testTiltedCylinder)
{
  Point_3  base(1, 1, 1);
  Vector_3 axis(1, 1, 1);
  Cylinder cyl(base, axis, 1.0, std::sqrt(3.0), 16);
  auto     mesh = cyl.generateSurfaceMesh();

  // Check that the top center is where we expect it to be
  Point_3 expected_top(2, 2, 2);
  bool    found_top = false;
  for (auto v : mesh.vertices()) {
    if (CGAL::squared_distance(mesh.point(v), expected_top) < 1e-10) {
      found_top = true;
      break;
    }
  }
  BOOST_CHECK(found_top);

  // Create PolyhedralSurface and output WKT for visual inspection
  PolyhedralSurface poly_surface(mesh);

  std::string expectedWkt(SFCGAL_TEST_DIRECTORY);
  expectedWkt += "/data/cylinder_expected_cpp2.wkt";
  std::ifstream efs(expectedWkt.c_str());
  BOOST_REQUIRE(efs.good());
  std::getline(efs, expectedWkt);

  BOOST_CHECK_EQUAL(poly_surface.asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testPolyhedron)
{
  Cylinder cyl(Point_3(0, 0, 0), Vector_3(0, 0, 1), 1.0, 2.0, 8);
  auto     polyhedron = cyl.generatePolyhedron();

  BOOST_CHECK_EQUAL(polyhedron.size_of_vertices(), cyl.numRadial() * 2 + 2);
  BOOST_CHECK_EQUAL(polyhedron.size_of_facets(), cyl.numRadial() * 4);

  // Create PolyhedralSurface from Polyhedron and check WKT output
  PolyhedralSurface poly_surface(polyhedron);

  std::string expectedWkt(SFCGAL_TEST_DIRECTORY);
  expectedWkt += "/data/cylinder_expected_cpp3.wkt";
  std::ifstream efs(expectedWkt.c_str());
  BOOST_REQUIRE(efs.good());
  std::getline(efs, expectedWkt);

  BOOST_CHECK_EQUAL(poly_surface.asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testClone)
{
  Point_3                   base(1, 2, 3);
  Vector_3                  axis(0, 1, 0);
  Cylinder                  cylinder(base, axis, 2.0, 5.0, 16);
  std::unique_ptr<Cylinder> cylinderCloned = cylinder.clone();

  BOOST_CHECK_EQUAL(cylinder, *cylinderCloned);
  BOOST_CHECK(algorithm::covers3D(cylinder.generatePolyhedralSurface(),
                                  cylinderCloned->generatePolyhedralSurface()));
}

BOOST_AUTO_TEST_SUITE_END()
