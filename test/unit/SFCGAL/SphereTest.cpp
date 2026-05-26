// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Sphere.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/io/wkt.h"
#include <boost/test/unit_test.hpp>
#include <cmath>

using namespace SFCGAL;

#include "../../test_config.h"

BOOST_AUTO_TEST_SUITE(SphereTests)

BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Sphere sphere;
  BOOST_CHECK_CLOSE(sphere.radius(), 1.0, 1e-6);
  BOOST_CHECK_EQUAL(sphere.numSubdivisions(), 2);
}

BOOST_AUTO_TEST_CASE(testCustomConstructor)
{
  Sphere sphere(2.0, 3);
  BOOST_CHECK_CLOSE(sphere.radius(), 2.0, 1e-6);
  BOOST_CHECK_EQUAL(sphere.numSubdivisions(), 3);
}

BOOST_AUTO_TEST_CASE(testSetters)
{
  Sphere sphere;
  sphere.setRadius(3.0);
  sphere.setNumSubdivisions(4);

  BOOST_CHECK_CLOSE(sphere.radius(), 3.0, 1e-6);
  BOOST_CHECK_EQUAL(sphere.numSubdivisions(), 4);

  // test setter and getter from Primitive
  sphere.setParameter(std::string("radius"), 12.3);
  BOOST_CHECK_CLOSE(sphere.radius(), 12.3, 1e-6);

  // sphere does not have a parameter called foo
  BOOST_CHECK_THROW(sphere.setParameter(std::string("foo"), 12.3), Exception);
  BOOST_CHECK_THROW(static_cast<void>(sphere.parameter("foo")), Exception);
}

BOOST_AUTO_TEST_CASE(testGeneratePolyhedron)
{
  Sphere sphere(1.0, 2); // 2 subdivisions
  auto   polyhedron = sphere.generatePolyhedron();

  // For icosahedron with 2 subdivisions: 162 vertices, 320 faces
  BOOST_CHECK_EQUAL(polyhedron.size_of_vertices(), 162);
  BOOST_CHECK_EQUAL(polyhedron.size_of_facets(), 320);
}

BOOST_AUTO_TEST_CASE(testGeneratePoints)
{
  Sphere sphere(1.0, 2); // 2 subdivisions
  auto   points = sphere.generatePoints();

  // The number of points should match the number of vertices in polyhedron
  BOOST_CHECK_EQUAL(points.size(), 162);

  // Check that all points are at a distance of 1 from the center
  for (const auto &point : points) {
    BOOST_CHECK_CLOSE(CGAL::sqrt(CGAL::to_double(
                          CGAL::squared_distance(Point_3(0, 0, 0), point))),
                      1.0, 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(testVolume)
{
  Sphere sphere(2.0, 3);
  double volume          = sphere.volume();
  double expected_volume = 4.0 / 3.0 * M_PI * 2.0 * 2.0 * 2.0;
  BOOST_CHECK_CLOSE(volume, expected_volume, 0.1); // 0.1% tolerance
}

BOOST_AUTO_TEST_CASE(testSurfaceArea)
{
  Sphere sphere(2.0, 3);
  double area          = sphere.area3D();
  double expected_area = 4.0 * M_PI * 2.0 * 2.0;
  BOOST_CHECK_CLOSE(area, expected_area, 0.1); // 0.1% tolerance
}

BOOST_AUTO_TEST_CASE(testWKT)
{
  Sphere            sphere(1.0, 1);
  auto              polyhedron = sphere.generatePolyhedron();
  PolyhedralSurface surface(polyhedron);
  std::string       wkt = surface.asText(1);

  BOOST_TEST_MESSAGE("Sphere WKT: " << wkt);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") == 0);
}

BOOST_AUTO_TEST_CASE(testClone)
{
  Sphere                  sphere(3.2, 6);
  std::unique_ptr<Sphere> sphereCloned = sphere.clone();

  BOOST_CHECK_EQUAL(sphere, *sphereCloned);
}

BOOST_AUTO_TEST_CASE(testTransform)
{
  Sphere                  sphere(3.0, 1);
  std::unique_ptr<Sphere> sphereTranslated = sphere.clone();
  BOOST_CHECK_EQUAL(sphere, *sphereTranslated);

  sphereTranslated->translate(Kernel::Vector_3(0.001, 0, 0));
  BOOST_CHECK_NE(sphere, *sphereTranslated);
  BOOST_CHECK(sphere.almostEqual(*sphereTranslated, 0.001));
  BOOST_CHECK(!sphere.almostEqual(*sphereTranslated, 0.0001));

  std::unique_ptr<Sphere> sphereTranslated2 = sphere.clone();
  sphereTranslated2->translate(Kernel::Vector_3(100, 10, 20));
  PolyhedralSurface polyhedral_surface =
      sphereTranslated2->generatePolyhedralSurface();

  std::string expectedWkt(SFCGAL_TEST_DIRECTORY);
  expectedWkt += "/data/sphere_translated_expected.wkt";
  std::ifstream efs(expectedWkt.c_str());
  BOOST_REQUIRE(efs.good());
  std::getline(efs, expectedWkt);

  BOOST_CHECK_EQUAL(polyhedral_surface.asText(1), expectedWkt);
}

BOOST_AUTO_TEST_SUITE_END()
