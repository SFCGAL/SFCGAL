// Copyright (c) 2025-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Cone.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/volume.h"
#include <boost/test/unit_test.hpp>

using namespace SFCGAL;

#include "../../test_config.h"

BOOST_AUTO_TEST_SUITE(ConeTests)

BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Cone cone = Cone();
  BOOST_CHECK_EQUAL(cone.height(), 1.0);
  BOOST_CHECK_EQUAL(cone.bottomRadius(), 1.0);
  BOOST_CHECK_EQUAL(cone.topRadius(), 0.0);
  BOOST_CHECK_EQUAL(cone.numRadial(), 32);
}

BOOST_AUTO_TEST_CASE(testCustomConstructor)
{
  Cone cone(57, 13, 44, 100);
  BOOST_CHECK_EQUAL(cone.bottomRadius(), 57);
  BOOST_CHECK_EQUAL(cone.topRadius(), 13);
  BOOST_CHECK_EQUAL(cone.height(), 44);
  BOOST_CHECK_EQUAL(cone.numRadial(), 100);
}

BOOST_AUTO_TEST_CASE(testSetters)
{
  Cone cone;
  cone.setBottomRadius(15.0);
  cone.setNumRadial(12);

  BOOST_CHECK_CLOSE(cone.topRadius(), 0.0, 1e-6);
  BOOST_CHECK_CLOSE(cone.bottomRadius(), 15.0, 1e-6);
  BOOST_CHECK_EQUAL(cone.numRadial(), 12);

  // test setter and getter from Primitive
  cone.setParameter(std::string("bottom_radius"), 5.3);
  BOOST_CHECK_CLOSE(cone.bottomRadius(), 5.3, 1e-6);

  cone.setParameter(std::string("top_radius"), 2.3);
  BOOST_CHECK_CLOSE(cone.topRadius(), 2.3, 1e-6);

  // cone does not have a parameter called foo
  BOOST_CHECK_THROW(cone.setParameter(std::string("foo"), 12.3), Exception);
  BOOST_CHECK_THROW(static_cast<void>(cone.parameter("foo")), Exception);

  // height, bottom and top radius cannot be negartive
  BOOST_CHECK_THROW(cone.setParameter(std::string("height"), -2.3), Exception);
  BOOST_CHECK_THROW(cone.setParameter(std::string("bottom_radius"), -20.38),
                    Exception);
  BOOST_CHECK_THROW(cone.setParameter(std::string("top_radius"), -5.2),
                    Exception);
}

BOOST_AUTO_TEST_CASE(testVolume)
{
  // small base, big top
  Cone   cone(23, 59, 68, 45);
  double expected_volume             = 382181.029;
  double expected_discretized_volume = 380940.437;
  BOOST_CHECK_CLOSE(cone.volume(), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(false), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(true), expected_discretized_volume, 1e-6);

  // big base, small top
  cone.setBottomRadius(59);
  cone.setTopRadius(23);
  BOOST_CHECK_CLOSE(cone.volume(), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(false), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(true), expected_discretized_volume, 1e-6);

  // same base and top (cylinder)
  cone.setBottomRadius(33);
  cone.setTopRadius(33);
  expected_volume             = 232641.219;
  expected_discretized_volume = 231886.045;
  BOOST_CHECK_CLOSE(cone.volume(), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(false), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(true), expected_discretized_volume, 1e-6);

  // top radius 0, real cone
  cone.setBottomRadius(36);
  cone.setTopRadius(0);
  expected_volume             = 92287.426;
  expected_discretized_volume = 91987.853;
  BOOST_CHECK_CLOSE(cone.volume(), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(false), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(true), expected_discretized_volume, 1e-6);

  // bottom radius 0, real cone upside-down
  cone.setBottomRadius(0);
  cone.setTopRadius(36);
  BOOST_CHECK_CLOSE(cone.volume(), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(false), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(true), expected_discretized_volume, 1e-6);

  // all radiuses 0
  cone.setBottomRadius(0);
  cone.setTopRadius(0);
  expected_volume             = 0;
  expected_discretized_volume = 0;
  BOOST_CHECK_EQUAL(cone.volume(), expected_volume);
  BOOST_CHECK_EQUAL(cone.volume(false), expected_volume);
  BOOST_CHECK_EQUAL(cone.volume(true), expected_discretized_volume);
}

BOOST_AUTO_TEST_CASE(testSurfaceArea)
{
  // small base, big top
  Cone   cone(23, 59, 68, 45);
  double expected_area             = 32418.741;
  double expected_discretized_area = 32351.199;
  BOOST_CHECK_CLOSE(cone.area3D(), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(false), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(true), expected_discretized_area, 1e-5);

  // big base, small top
  cone.setBottomRadius(59);
  cone.setTopRadius(23);
  BOOST_CHECK_CLOSE(cone.area3D(), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(false), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(true), expected_discretized_area, 1e-5);

  // same base and top (cylinder)
  cone.setBottomRadius(32);
  expected_area             = 16730.913;
  expected_discretized_area = 16704.954;
  BOOST_CHECK_CLOSE(cone.area3D(), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(false), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(true), expected_discretized_area, 1e-5);

  // top radius 0, real cone
  cone.setBottomRadius(87);
  cone.setTopRadius(0);
  cone.setHeight(13);
  cone.setNumRadial(30);
  expected_area             = 47821.428;
  expected_discretized_area = 47475.459;
  BOOST_CHECK_CLOSE(cone.area3D(), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(false), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(true), expected_discretized_area, 1e-5);

  // bottom radius 0, real cone upside-down
  cone.setTopRadius(32);
  cone.setBottomRadius(0);
  expected_area             = 6689.313;
  expected_discretized_area = 6643.212;
  BOOST_CHECK_CLOSE(cone.area3D(), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(false), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(true), expected_discretized_area, 1e-5);

  // all radiuses 0
  cone.setBottomRadius(0);
  cone.setTopRadius(0);
  expected_area             = 0;
  expected_discretized_area = 0;
  BOOST_CHECK_EQUAL(cone.area3D(), expected_area);
  BOOST_CHECK_EQUAL(cone.area3D(false), expected_area);
  BOOST_CHECK_EQUAL(cone.area3D(true), expected_discretized_area);
}

BOOST_AUTO_TEST_CASE(testPolyhedralSurface)
{
  // Create PolyhedralSurface from Polyhedron and check WKT output
  Cone cone(54, 21, 87, 7);
  auto polyhedral_surface = cone.generatePolyhedralSurface();

  BOOST_CHECK(algorithm::isValid(polyhedral_surface));

  BOOST_CHECK_EQUAL(
      polyhedral_surface.asText(1),
      "POLYHEDRALSURFACE Z (((54.0 0.0 -43.5,33.7 -42.2 -43.5,-12.0 -52.6 "
      "-43.5,-48.7 -23.4 -43.5,-48.7 23.4 -43.5,-12.0 52.6 -43.5,33.7 42.2 "
      "-43.5,54.0 0.0 -43.5)),((21.0 0.0 43.5,13.1 16.4 43.5,-4.7 20.5 "
      "43.5,-18.9 9.1 43.5,-18.9 -9.1 43.5,-4.7 -20.5 43.5,13.1 -16.4 "
      "43.5,21.0 0.0 43.5)),((54.0 0.0 -43.5,21.0 0.0 43.5,13.1 -16.4 "
      "43.5,33.7 -42.2 -43.5,54.0 0.0 -43.5)),((33.7 -42.2 -43.5,13.1 -16.4 "
      "43.5,-4.7 -20.5 43.5,-12.0 -52.6 -43.5,33.7 -42.2 -43.5)),((-12.0 -52.6 "
      "-43.5,-4.7 -20.5 43.5,-18.9 -9.1 43.5,-48.7 -23.4 -43.5,-12.0 -52.6 "
      "-43.5)),((-48.7 -23.4 -43.5,-18.9 -9.1 43.5,-18.9 9.1 43.5,-48.7 23.4 "
      "-43.5,-48.7 -23.4 -43.5)),((-48.7 23.4 -43.5,-18.9 9.1 43.5,-4.7 20.5 "
      "43.5,-12.0 52.6 -43.5,-48.7 23.4 -43.5)),((-12.0 52.6 -43.5,-4.7 20.5 "
      "43.5,13.1 16.4 43.5,33.7 42.2 -43.5,-12.0 52.6 -43.5)),((33.7 42.2 "
      "-43.5,13.1 16.4 43.5,21.0 0.0 43.5,54.0 0.0 -43.5,33.7 42.2 -43.5)))");

  BOOST_CHECK_CLOSE(CGAL::to_double(algorithm::volume(
                        Solid(cone.generatePolyhedralSurface()))),
                    cone.volume(true), 1e-6);

  BOOST_CHECK_CLOSE(
      CGAL::to_double(algorithm::area3D(cone.generatePolyhedralSurface())),
      cone.area3D(true), 1e-6);
}

BOOST_AUTO_TEST_CASE(testGetSetHeight)
{
  Cone cone(6, 7, 8, 9);
  BOOST_CHECK_EQUAL(cone.height(), 8);

  cone.setHeight(6.904);
  BOOST_CHECK_EQUAL(cone.height(), 6.904);
}

BOOST_AUTO_TEST_CASE(testGetSetBottomRadius)
{
  Cone cone(6, 7, 8, 9);
  BOOST_CHECK_EQUAL(cone.bottomRadius(), 6);

  cone.setBottomRadius(5.904);
  BOOST_CHECK_EQUAL(cone.bottomRadius(), 5.904);
}

BOOST_AUTO_TEST_CASE(testGetSetTopRadius)
{
  Cone cone(6, 7, 8, 9);
  BOOST_CHECK_EQUAL(cone.topRadius(), 7);

  cone.setTopRadius(6.904);
  BOOST_CHECK_EQUAL(cone.topRadius(), 6.904);
}

BOOST_AUTO_TEST_CASE(testGetSetNumRadial)
{
  Cone cone(6, 7, 8, 9);
  BOOST_CHECK_EQUAL(cone.numRadial(), 9);

  cone.setNumRadial(54);
  BOOST_CHECK_EQUAL(cone.numRadial(), 54);
}

BOOST_AUTO_TEST_CASE(testClone)
{
  Cone                  cone(3.5, 34.5, 3.3);
  std::unique_ptr<Cone> coneCloned = cone.clone();

  BOOST_CHECK_EQUAL(cone, *coneCloned);
  BOOST_CHECK(algorithm::covers3D(cone.generatePolyhedralSurface(),
                                  coneCloned->generatePolyhedralSurface()));
}

BOOST_AUTO_TEST_CASE(testTransform)
{
  Cone                  cone(3.0, 0.1, 5.0);
  std::unique_ptr<Cone> coneTranslated = cone.clone();
  BOOST_CHECK_EQUAL(cone, *coneTranslated);

  coneTranslated->translate(Kernel::Vector_3(0.001, 0, 0));
  BOOST_CHECK_NE(cone, *coneTranslated);
  BOOST_CHECK(cone.almostEqual(*coneTranslated, 0.001));
  BOOST_CHECK(!cone.almostEqual(*coneTranslated, 0.0001));

  // translation
  std::unique_ptr<Cone> coneTranslated2 = cone.clone();
  coneTranslated2->translate(Kernel::Vector_3(100, 10, 20));
  PolyhedralSurface polyhedral_surface_translated =
      coneTranslated2->generatePolyhedralSurface();

  std::string expectedWktTranslated(SFCGAL_TEST_DIRECTORY);
  expectedWktTranslated += "/data/cone_translated_expected.wkt";
  std::ifstream efsT(expectedWktTranslated.c_str());
  BOOST_REQUIRE(efsT.good());
  std::getline(efsT, expectedWktTranslated);

  BOOST_CHECK_EQUAL(polyhedral_surface_translated.asText(1),
                    expectedWktTranslated);

  // rotation
  std::unique_ptr<Cone> coneRotated = coneTranslated2->clone();
  coneRotated->rotate(45 * M_PI / 180., Kernel::Vector_3(0, 1, 0));
  PolyhedralSurface polyhedral_surface_rotated =
      coneRotated->generatePolyhedralSurface();

  std::string expectedWktRotated(SFCGAL_TEST_DIRECTORY);
  expectedWktRotated += "/data/cone_rotated_expected.wkt";
  std::ifstream efsR(expectedWktRotated.c_str());
  BOOST_REQUIRE(efsR.good());
  std::getline(efsR, expectedWktRotated);

  BOOST_CHECK_EQUAL(polyhedral_surface_rotated.asText(1), expectedWktRotated);

  // scale
  std::unique_ptr<Cone> coneScaled = cone.clone();
  BOOST_CHECK_CLOSE(coneScaled->volume(), 48.74705, 1e-5);
  BOOST_CHECK_CLOSE(coneScaled->area3D(), 84.59815, 1e-5);
  coneScaled->scale(Kernel::Vector_3(2, 2, 2));
  BOOST_CHECK_CLOSE(coneScaled->volume(), 48.74705 * 8.0, 1e-5);
  BOOST_CHECK_CLOSE(coneScaled->area3D(), 84.59815 * 4.0, 1e-5);
  PolyhedralSurface polyhedral_surface_scaled =
      coneScaled->generatePolyhedralSurface();

  std::string expectedWktScaled(SFCGAL_TEST_DIRECTORY);
  expectedWktScaled += "/data/cone_scaled_expected.wkt";
  std::ifstream efsS(expectedWktScaled.c_str());
  BOOST_REQUIRE(efsS.good());
  std::getline(efsS, expectedWktScaled);

  BOOST_CHECK_EQUAL(polyhedral_surface_scaled.asText(1), expectedWktScaled);

  std::unique_ptr<Cone> coneScaled2 = cone.clone();
  coneScaled2->scale(Kernel::Vector_3(2, 1, 3));
  BOOST_CHECK_CLOSE(coneScaled2->volume(), 48.74705 * 6.0, 1e-5);
  BOOST_CHECK_CLOSE(coneScaled2->area3D(), 288.40878, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()
