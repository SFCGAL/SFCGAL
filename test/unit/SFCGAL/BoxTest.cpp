// Copyright (c) 2025-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Box.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/isValid.h"
#include <boost/test/unit_test.hpp>

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(BoxTests)

BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Box box;
  BOOST_CHECK_EQUAL(box.xExtent(), 1.0);
  BOOST_CHECK_EQUAL(box.yExtent(), 1.0);
  BOOST_CHECK_EQUAL(box.zExtent(), 1.0);
}

BOOST_AUTO_TEST_CASE(testCustomConstructor)
{
  Box box(14, 10, 57);
  BOOST_CHECK_EQUAL(box.xExtent(), 14);
  BOOST_CHECK_EQUAL(box.yExtent(), 10);
  BOOST_CHECK_EQUAL(box.zExtent(), 57);
}

BOOST_AUTO_TEST_CASE(testSetters)
{
  Box box;
  box.setXExtent(3.0);
  box.setYExtent(5.0);
  box.setZExtent(7.0);

  BOOST_CHECK_CLOSE(box.xExtent(), 3.0, 1e-6);
  BOOST_CHECK_CLOSE(box.yExtent(), 5.0, 1e-6);
  BOOST_CHECK_CLOSE(box.zExtent(), 7.0, 1e-6);

  // extent cannot be negative
  BOOST_CHECK_THROW(box.setXExtent(-1.2), Exception);
  BOOST_CHECK_THROW(box.setYExtent(-2.2), Exception);
  BOOST_CHECK_THROW(box.setZExtent(-5.4), Exception);
}

BOOST_AUTO_TEST_CASE(testVolume)
{
  Box    box(23, 9, 91);
  double volume = box.volume();
  BOOST_CHECK_EQUAL(volume, 18837);
}

BOOST_AUTO_TEST_CASE(testSurfaceArea)
{
  Box    box(14, 4, 96);
  double area = box.area3D();
  BOOST_CHECK_EQUAL(area, 3568);
}

BOOST_AUTO_TEST_CASE(testPolyhedralSurface)
{
  // Create PolyhedralSurface from Polyhedron and check WKT output
  Box  box(24, 10, 82);
  auto polyhedral_surface = box.generatePolyhedralSurface();

  BOOST_CHECK(algorithm::isValid(polyhedral_surface));

  BOOST_CHECK_EQUAL(
      polyhedral_surface.asText(1),
      "POLYHEDRALSURFACE Z (((-12.0 -5.0 -41.0,-12.0 5.0 -41.0,12.0 5.0 "
      "-41.0,12.0 -5.0 -41.0,-12.0 -5.0 -41.0)),((-12.0 -5.0 41.0,12.0 -5.0 "
      "41.0,12.0 5.0 41.0,-12.0 5.0 41.0,-12.0 -5.0 41.0)),((-12.0 -5.0 "
      "-41.0,12.0 -5.0 -41.0,12.0 -5.0 41.0,-12.0 -5.0 41.0,-12.0 -5.0 "
      "-41.0)),((-12.0 5.0 -41.0,-12.0 5.0 41.0,12.0 5.0 41.0,12.0 5.0 "
      "-41.0,-12.0 5.0 -41.0)),((12.0 -5.0 -41.0,12.0 5.0 -41.0,12.0 5.0 "
      "41.0,12.0 -5.0 41.0,12.0 -5.0 -41.0)),((-12.0 -5.0 -41.0,-12.0 -5.0 "
      "41.0,-12.0 5.0 41.0,-12.0 5.0 -41.0,-12.0 -5.0 -41.0)))");
}

BOOST_AUTO_TEST_CASE(testGetSetExtents)
{
  Box box(3.5, 62.78, 1.621);
  BOOST_CHECK_EQUAL(box.xExtent(), 3.5);
  BOOST_CHECK_EQUAL(box.yExtent(), 62.78);
  BOOST_CHECK_EQUAL(box.zExtent(), 1.621);

  box.setXExtent(6.904);
  box.setYExtent(9.44);
  box.setZExtent(543.8209);

  BOOST_CHECK_EQUAL(box.xExtent(), 6.904);
  BOOST_CHECK_EQUAL(box.yExtent(), 9.44);
  BOOST_CHECK_EQUAL(box.zExtent(), 543.8209);
}

BOOST_AUTO_TEST_CASE(testClone)
{
  Box                  box(3.5, 34.5, 3.3);
  std::unique_ptr<Box> boxCloned = box.clone();

  BOOST_CHECK_EQUAL(box, *boxCloned);
  BOOST_CHECK(algorithm::covers3D(box.generatePolyhedralSurface(),
                                  boxCloned->generatePolyhedralSurface()));

  BOOST_CHECK_EQUAL(box.transformation(),
                    Kernel::Aff_transformation_3(CGAL::IDENTITY));
  Kernel::Aff_transformation_3 translation(CGAL::TRANSLATION,
                                           SFCGAL::Kernel::Vector_3(3, 2, 1));
  box.setTransformation(translation);
  BOOST_CHECK_EQUAL(box.transformation(), translation);
  std::unique_ptr<Box> box2 = box.clone();
  BOOST_CHECK_EQUAL(box, *box2);
  BOOST_CHECK_EQUAL(box.transformation(), box2->transformation());
}

BOOST_AUTO_TEST_CASE(testTransform)
{
  Box                  box(3.0, 4.0, 5.0);
  std::unique_ptr<Box> boxTranslated = box.clone();
  BOOST_CHECK_EQUAL(box, *boxTranslated);

  boxTranslated->translate(Kernel::Vector_3(0.001, 0, 0));
  BOOST_CHECK_NE(box, *boxTranslated);
  BOOST_CHECK(box.almostEqual(*boxTranslated, 0.001));
  BOOST_CHECK(!box.almostEqual(*boxTranslated, 0.0001));

  // translation
  std::unique_ptr<Box> boxTranslated2 = box.clone();
  boxTranslated2->translate(Kernel::Vector_3(100, 10, 20));
  PolyhedralSurface polyhedral_surface_translated =
      boxTranslated2->generatePolyhedralSurface();

  BOOST_CHECK_EQUAL(
      polyhedral_surface_translated.asText(1),
      "POLYHEDRALSURFACE Z (((98.5 8.0 17.5,98.5 12.0 17.5,101.5 12.0 "
      "17.5,101.5 8.0 17.5,98.5 8.0 17.5)),((98.5 8.0 22.5,101.5 8.0 "
      "22.5,101.5 12.0 22.5,98.5 12.0 22.5,98.5 8.0 22.5)),((98.5 8.0 "
      "17.5,101.5 8.0 17.5,101.5 8.0 22.5,98.5 8.0 22.5,98.5 8.0 17.5)),((98.5 "
      "12.0 17.5,98.5 12.0 22.5,101.5 12.0 22.5,101.5 12.0 17.5,98.5 12.0 "
      "17.5)),((101.5 8.0 17.5,101.5 12.0 17.5,101.5 12.0 22.5,101.5 8.0 "
      "22.5,101.5 8.0 17.5)),((98.5 8.0 17.5,98.5 8.0 22.5,98.5 12.0 22.5,98.5 "
      "12.0 17.5,98.5 8.0 17.5)))");
}

BOOST_AUTO_TEST_SUITE_END()
