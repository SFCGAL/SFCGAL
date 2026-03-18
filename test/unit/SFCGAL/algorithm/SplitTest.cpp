// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <memory>

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/split.h"
#include "SFCGAL/io/wkt.h"
#include "SFCGAL/primitive3d/Cube.h"

using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_Split)

BOOST_AUTO_TEST_CASE(testSplit_Empty)
{
  const PolyhedralSurface emptyGeom;
  BOOST_CHECK(emptyGeom.isEmpty());
  std::unique_ptr<GeometryCollection> geomSplit = algorithm::split(
      emptyGeom, Point(0, 0, 0), Kernel::Vector_3(1, 0, -1), false);
  BOOST_CHECK(geomSplit->isEmpty());
}

BOOST_AUTO_TEST_CASE(testSplit_Cube)
{
  Cube              cube(57);
  PolyhedralSurface phs = cube.generatePolyhedralSurface();
  BOOST_CHECK(!phs.isEmpty());

  std::unique_ptr<GeometryCollection> phsSplit =
      algorithm::split(phs, Point(0, 0, 0), Kernel::Vector_3(1, 0, -1), false);
  BOOST_CHECK_EQUAL(phsSplit->numGeometries(), 2);

  const std::string expectedGeom1Wkt =
      "POLYHEDRALSURFACE Z ("
      "((0.0 0.0 0.0,0.0 57.0 0.0,57.0 57.0 0.0,0.0 0.0 0.0)),"
      "((57.0 57.0 0.0,57.0 0.0 0.0,0.0 0.0 0.0,57.0 57.0 0.0)),"
      "((28.5 0.0 28.5,0.0 0.0 0.0,57.0 0.0 0.0,28.5 0.0 28.5)),"
      "((57.0 57.0 57.0,57.0 57.0 0.0,0.0 57.0 0.0,57.0 57.0 57.0)),"
      "((57.0 57.0 57.0,57.0 0.0 57.0,57.0 57.0 0.0,57.0 57.0 57.0)),"
      "((57.0 57.0 0.0,57.0 0.0 57.0,57.0 0.0 0.0,57.0 57.0 0.0)),"
      "((57.0 0.0 57.0,28.5 0.0 28.5,57.0 0.0 0.0,57.0 0.0 57.0)))";

  const std::unique_ptr<Geometry> expectedGeom1 = io::readWkt(expectedGeom1Wkt);
  BOOST_CHECK(algorithm::covers3D(phsSplit->geometryN(0), *expectedGeom1));

  const std::string expectedGeom2Wkt =
      "POLYHEDRALSURFACE Z ("
      "((57.0 0.0 57.0,57.0 57.0 57.0,0.0 57.0 57.0,57.0 0.0 57.0)),"
      "((57.0 0.0 57.0,0.0 57.0 57.0,0.0 0.0 57.0,57.0 0.0 57.0)),"
      "((28.5 0.0 28.5,57.0 0.0 57.0,0.0 0.0 57.0,28.5 0.0 28.5)),"
      "((57.0 57.0 57.0,0.0 57.0 0.0,0.0 57.0 57.0,57.0 57.0 57.0)),"
      "((0.0 57.0 57.0,0.0 0.0 0.0,0.0 0.0 57.0,0.0 57.0 57.0)),"
      "((0.0 57.0 57.0,0.0 57.0 0.0,0.0 0.0 0.0,0.0 57.0 57.0)),"
      "((0.0 0.0 0.0,28.5 0.0 28.5,0.0 0.0 57.0,0.0 0.0 0.0)))";

  const std::unique_ptr<Geometry> expectedGeom2 = io::readWkt(expectedGeom2Wkt);
  BOOST_CHECK(algorithm::covers3D(phsSplit->geometryN(1), *expectedGeom2));
}

BOOST_AUTO_TEST_CASE(testSplit_UShape)
{
  const std::string inputWkt(
      "PolyhedralSurface Z("
      "((0 0 0,0 10 0,2 10 0,2 2 0,8 2 0,8 10 0,10 10 0,10 0 0,0 0 0)),"
      "((0 0 2,10 0 2,10 10 2,8 10 2,8 2 2,2 2 2,2 10 2,0 10 2,0 0 2)),"
      "((0 0 0,0 0 2,0 10 2,0 10 0,0 0 0)),"
      "((0 10 0,0 10 2,2 10 2,2 10 0,0 10 0)),"
      "((2 10 0,2 10 2,2 2 2,2 2 0,2 10 0)),"
      "((2 2 0,2 2 2,8 2 2,8 2 0,2 2 0)),"
      "((8 2 0,8 2 2,8 10 2,8 10 0,8 2 0)),"
      "((8 10 0,8 10 2,10 10 2,10 10 0,8 10 0)),"
      "((10 10 0,10 10 2,10 0 2,10 0 0,10 10 0)),"
      "((10 0 0,10 0 2,0 0 2,0 0 0,10 0 0)))");
  const std::unique_ptr<Geometry> inputGeom(io::readWkt(inputWkt));
  BOOST_CHECK(!inputGeom->isEmpty());
  BOOST_CHECK_EQUAL(inputGeom->geometryTypeId(),
                    SFCGAL::GeometryType::TYPE_POLYHEDRALSURFACE);

  std::unique_ptr<GeometryCollection> phsSplit = algorithm::split(
      *inputGeom, Point(1, 6, 0), Kernel::Vector_3(0, 1, 0), false);

  const std::string expectedGeom1Wkt =
      "POLYHEDRALSURFACE Z (((1 6 0,2 6 0,2 2 0,1 6 0)),((2 2 0,8 2 0,0 0 0,2 "
      "2 0)),((8 2 0,10 0 0,0 0 0,8 2 0)),((9 6 0,10 6 0,8 2 0,9 6 0)),((10 0 "
      "2,8 2 2,0 0 2,10 0 2)),((8 2 2,2 2 2,0 0 2,8 2 2)),((1 6 2,0 6 2,2 2 "
      "2,1 6 2)),((0 6 2,0 6 1,0 0 2,0 6 2)),((2 6 0,2 6 1,2 2 0,2 6 0)),((8 2 "
      "2,2 2 0,2 2 2,8 2 2)),((8 2 2,8 2 0,2 2 0,8 2 2)),((8 6 2,8 6 1,8 2 2,8 "
      "6 2)),((10 6 0,10 6 1,10 0 0,10 6 0)),((0 0 2,10 0 0,10 0 2,0 0 2)),((0 "
      "0 2,0 0 0,10 0 0,0 0 2)),((0 0 2,0 6 1,0 0 0,0 0 2)),((10 6 1,10 6 2,10 "
      "0 2,10 6 1)),((2 6 1,2 6 2,2 2 2,2 6 1)),((0 6 1,0 6 0,0 0 0,0 6 "
      "1)),((1 6 0,0 0 0,0 6 0,1 6 0)),((0 0 0,1 6 0,2 2 0,0 0 0)),((2 2 0,2 6 "
      "1,2 2 2,2 2 0)),((8 6 0,9 6 0,8 2 0,8 6 0)),((8 2 0,10 6 0,10 0 0,8 2 "
      "0)),((8 6 1,8 6 0,8 2 0,8 6 1)),((10 0 0,10 6 1,10 0 2,10 0 0)),((9 6 "
      "2,8 6 2,8 2 2,9 6 2)),((8 2 2,8 6 1,8 2 0,8 2 2)),((2 6 2,1 6 2,2 2 2,2 "
      "6 2)),((9 6 2,10 0 2,10 6 2,9 6 2)),((10 0 2,9 6 2,8 2 2,10 0 2)),((2 2 "
      "2,0 6 2,0 0 2,2 2 2)))";

  const std::unique_ptr<Geometry> expectedGeom1 = io::readWkt(expectedGeom1Wkt);
  BOOST_CHECK(algorithm::covers3D(phsSplit->geometryN(0), *expectedGeom1));

  const std::string expectedGeom2Wkt =
      "POLYHEDRALSURFACE Z (((1 6 0,0 6 0,0 10 0,1 6 0)),((1 6 2,2 6 2,0 10 "
      "2,1 6 2)),((0 6 0,0 6 1,0 10 0,0 6 0)),((2 10 2,0 10 0,0 10 2,2 10 "
      "2)),((2 10 2,2 10 0,0 10 0,2 10 2)),((2 6 2,2 6 1,2 10 2,2 6 2)),((0 6 "
      "1,0 6 2,0 10 2,0 6 1)),((2 10 2,2 6 1,2 10 0,2 10 2)),((1 6 0,2 10 0,2 "
      "6 0,1 6 0)),((2 10 0,1 6 0,0 10 0,2 10 0)),((0 10 0,0 6 1,0 10 2,0 10 "
      "0)),((2 6 1,2 6 0,2 10 0,2 6 1)),((0 10 2,2 6 2,2 10 2,0 10 2)),((0 6 "
      "2,1 6 2,0 10 2,0 6 2)))";

  const std::unique_ptr<Geometry> expectedGeom2 = io::readWkt(expectedGeom2Wkt);
  BOOST_CHECK(algorithm::covers3D(phsSplit->geometryN(1), *expectedGeom2));

  const std::string expectedGeom3Wkt =
      "POLYHEDRALSURFACE Z (((9 6 0,8 6 0,10 10 0,9 6 0)),((8 6 2,9 6 2,8 10 "
      "2,8 6 2)),((9 6 2,10 6 2,10 10 2,9 6 2)),((8 6 0,8 6 1,8 10 0,8 6 "
      "0)),((10 10 2,8 10 0,8 10 2,10 10 2)),((10 10 2,10 10 0,8 10 0,10 10 "
      "2)),((10 6 2,10 6 1,10 10 2,10 6 2)),((10 10 2,10 6 1,10 10 0,10 10 "
      "2)),((10 10 0,8 6 0,8 10 0,10 10 0)),((10 6 0,9 6 0,10 10 0,10 6 "
      "0)),((8 10 0,8 6 1,8 10 2,8 10 0)),((10 6 1,10 6 0,10 10 0,10 6 1)),((8 "
      "10 2,9 6 2,10 10 2,8 10 2)),((8 6 1,8 6 2,8 10 2,8 6 1)))";

  const std::unique_ptr<Geometry> expectedGeom3 = io::readWkt(expectedGeom3Wkt);
  BOOST_CHECK(algorithm::covers3D(phsSplit->geometryN(2), *expectedGeom3));
}

BOOST_AUTO_TEST_CASE(testSplit_Tin)
{
  const std::string inputWkt(
      "TIN Z (((0 0 0,5 0 0,2 2 3,0 0 0)),((5 0 0,5 5 0,2 2 3,5 0 "
      "0)),((5 5 0,0 5 0,2 2 3,5 5 0)),((0 5 0,0 0 0,2 2 3,0 5 0)),((0 0 0,5 0 "
      "0,5 5 0,0 0 0)),((0 0 0,5 5 0,0 5 0,0 0 0)))");
  const std::unique_ptr<Geometry> inputGeom(io::readWkt(inputWkt));
  BOOST_CHECK(!inputGeom->isEmpty());
  BOOST_CHECK_EQUAL(inputGeom->geometryTypeId(),
                    SFCGAL::GeometryType::TYPE_TRIANGULATEDSURFACE);

  std::unique_ptr<GeometryCollection> geomSplit = algorithm::split(
      *inputGeom, Point(2, 2, 1.2), Kernel::Vector_3(0, 0, 1), false);
  BOOST_CHECK_EQUAL(geomSplit->numGeometries(), 2);

  const std::string expectedGeom1Wkt =
      "POLYHEDRALSURFACE Z (((3.8 0.8 1.2,0.8 0.8 1.2,5.0 0.0 0.0,3.8 0.8 "
      "1.2)),((0.8 0.8 1.2,0.8 3.8 1.2,0.0 0.0 0.0,0.8 0.8 1.2)),((3.8 3.8 "
      "1.2,0.0 5.0 0.0,0.8 3.8 1.2,3.8 3.8 1.2)),((0.0 5.0 0.0,3.8 3.8 1.2,5.0 "
      "5.0 0.0,0.0 5.0 0.0)),((0.0 0.0 0.0,0.8 3.8 1.2,0.0 5.0 0.0,0.0 0.0 "
      "0.0)),((3.8 0.8 1.2,5.0 5.0 0.0,3.8 3.8 1.2,3.8 0.8 1.2)),((5.0 5.0 "
      "0.0,3.8 0.8 1.2,5.0 0.0 0.0,5.0 5.0 0.0)),((5.0 0.0 0.0,0.8 0.8 1.2,0.0 "
      "0.0 0.0,5.0 0.0 0.0)))";

  const std::unique_ptr<Geometry> expectedGeom1 = io::readWkt(expectedGeom1Wkt);
  BOOST_CHECK(algorithm::covers3D(geomSplit->geometryN(0), *expectedGeom1));

  const std::string expectedGeom2Wkt =
      "POLYHEDRALSURFACE Z (((3.8 0.8 1.2,3.8 3.8 1.2,2.0 2.0 3.0,3.8 0.8 "
      "1.2)),((3.8 3.8 1.2,0.8 3.8 1.2,2.0 2.0 3.0,3.8 3.8 1.2)),((0.8 3.8 "
      "1.2,0.8 0.8 1.2,2.0 2.0 3.0,0.8 3.8 1.2)),((0.8 0.8 1.2,3.8 0.8 1.2,2.0 "
      "2.0 3.0,0.8 0.8 1.2)))";

  const std::unique_ptr<Geometry> expectedGeom2 = io::readWkt(expectedGeom2Wkt);
  BOOST_CHECK(algorithm::covers3D(geomSplit->geometryN(1), *expectedGeom2));
}

BOOST_AUTO_TEST_CASE(testSplit_GeometryCollection)
{
  GeometryCollection geomCollection;
  geomCollection.addGeometry(std::make_unique<Point>(2.0, 3.0, 5.0));
  geomCollection.addGeometry(io::readWkt(
      "POLYHEDRALSURFACE Z (((0 0 0,0 10 0,10 10 0,10 0 0,0 0 0)),((0 0 0,10 0 "
      "0,10 0 5,0 0 5,0 0 0)),((10 0 0,10 10 0,10 10 5,10 0 5,10 0 0)),((10 10 "
      "0,0 10 0,0 10 5,10 10 5,10 10 0)),((0 10 0,0 0 0,0 0 5,0 10 5,0 10 "
      "0)),((0 0 5,10 0 5,5 0 8,0 0 5)),((10 10 5,0 10 5,5 10 8,10 10 5)),((0 "
      "0 5,5 0 8,5 10 8,0 10 5,0 0 5)),((10 0 5,10 10 5,5 10 8,5 0 8,10 0 "
      "5)))"));

  // Do not close geometries
  std::unique_ptr<GeometryCollection> geomSplitOpen = algorithm::split(
      geomCollection, Point(0, 0, 2.5), Kernel::Vector_3(0, 0, 1), false);

  BOOST_CHECK_EQUAL(geomSplitOpen->numGeometries(), 2);

  const std::string expectedGeom1Wkt =
      "POLYHEDRALSURFACE Z (((10.0 10.0 0.0,0.0 0.0 0.0,0.0 10.0 0.0,10.0 10.0 "
      "0.0)),((10.0 10.0 0.0,10.0 0.0 0.0,0.0 0.0 0.0,10.0 10.0 0.0)),((5.0 "
      "0.0 2.5,0.0 0.0 2.5,10.0 0.0 0.0,5.0 0.0 2.5)),((0.0 0.0 2.5,0.0 5.0 "
      "2.5,0.0 0.0 0.0,0.0 0.0 2.5)),((0.0 10.0 2.5,0.0 0.0 0.0,0.0 5.0 "
      "2.5,0.0 10.0 2.5)),((0.0 0.0 0.0,0.0 10.0 2.5,0.0 10.0 0.0,0.0 0.0 "
      "0.0)),((0.0 10.0 2.5,5.0 10.0 2.5,0.0 10.0 0.0,0.0 10.0 2.5)),((10.0 "
      "10.0 2.5,0.0 10.0 0.0,5.0 10.0 2.5,10.0 10.0 2.5)),((0.0 10.0 0.0,10.0 "
      "10.0 2.5,10.0 10.0 0.0,0.0 10.0 0.0)),((10.0 10.0 2.5,10.0 5.0 2.5,10.0 "
      "10.0 0.0,10.0 10.0 2.5)),((10.0 0.0 2.5,10.0 10.0 0.0,10.0 5.0 2.5,10.0 "
      "0.0 2.5)),((10.0 10.0 0.0,10.0 0.0 2.5,10.0 0.0 0.0,10.0 10.0 "
      "0.0)),((10.0 0.0 0.0,0.0 0.0 2.5,0.0 0.0 0.0,10.0 0.0 0.0)),((10.0 0.0 "
      "2.5,5.0 0.0 2.5,10.0 0.0 0.0,10.0 0.0 2.5)))";

  const std::unique_ptr<Geometry> expectedGeom1 = io::readWkt(expectedGeom1Wkt);
  BOOST_CHECK(algorithm::covers3D(geomSplitOpen->geometryN(0), *expectedGeom1));

  const std::string expectedGeom2Wkt =
      "POLYHEDRALSURFACE Z (((5.0 0.0 2.5,10.0 0.0 2.5,0.0 0.0 5.0,5.0 0.0 "
      "2.5)),((10.0 5.0 2.5,10.0 10.0 2.5,10.0 0.0 5.0,10.0 5.0 2.5)),((10.0 "
      "0.0 2.5,10.0 5.0 2.5,10.0 0.0 5.0,10.0 0.0 2.5)),((5.0 10.0 2.5,0.0 "
      "10.0 2.5,10.0 10.0 5.0,5.0 10.0 2.5)),((10.0 10.0 2.5,5.0 10.0 2.5,10.0 "
      "10.0 5.0,10.0 10.0 2.5)),((0.0 10.0 2.5,0.0 5.0 2.5,0.0 10.0 5.0,0.0 "
      "10.0 2.5)),((5.0 0.0 8.0,0.0 0.0 5.0,10.0 0.0 5.0,5.0 0.0 8.0)),((5.0 "
      "10.0 8.0,10.0 10.0 5.0,0.0 10.0 5.0,5.0 10.0 8.0)),((5.0 10.0 8.0,0.0 "
      "0.0 5.0,5.0 0.0 8.0,5.0 10.0 8.0)),((5.0 10.0 8.0,0.0 10.0 5.0,0.0 0.0 "
      "5.0,5.0 10.0 8.0)),((5.0 10.0 8.0,5.0 0.0 8.0,10.0 10.0 5.0,5.0 10.0 "
      "8.0)),((10.0 10.0 5.0,5.0 0.0 8.0,10.0 0.0 5.0,10.0 10.0 5.0)),((10.0 "
      "10.0 5.0,0.0 10.0 2.5,0.0 10.0 5.0,10.0 10.0 5.0)),((10.0 0.0 5.0,10.0 "
      "10.0 2.5,10.0 10.0 5.0,10.0 0.0 5.0)),((0.0 0.0 2.5,0.0 10.0 5.0,0.0 "
      "5.0 2.5,0.0 0.0 2.5)),((0.0 10.0 5.0,0.0 0.0 2.5,0.0 0.0 5.0,0.0 10.0 "
      "5.0)),((0.0 0.0 2.5,5.0 0.0 2.5,0.0 0.0 5.0,0.0 0.0 2.5)),((0.0 0.0 "
      "5.0,10.0 0.0 2.5,10.0 0.0 5.0,0.0 0.0 5.0)))";

  const std::unique_ptr<Geometry> expectedGeom2 = io::readWkt(expectedGeom2Wkt);
  BOOST_CHECK(algorithm::covers3D(geomSplitOpen->geometryN(1), *expectedGeom2));

  // Close geometries
  std::unique_ptr<GeometryCollection> geomSplitClose = algorithm::split(
      geomCollection, Point(0, 0, 2.5), Kernel::Vector_3(0, 0, 1), true);

  BOOST_CHECK_EQUAL(geomSplitClose->numGeometries(), 2);

  const std::string expectedGeom3Wkt =
      "POLYHEDRALSURFACE Z (((10.0 10.0 0.0,0.0 0.0 0.0,0.0 10.0 0.0,10.0 10.0 "
      "0.0)),((10.0 10.0 0.0,10.0 0.0 0.0,0.0 0.0 0.0,10.0 10.0 0.0)),((0.0 "
      "0.0 2.5,10.0 0.0 0.0,5.0 0.0 2.5,0.0 0.0 2.5)),((0.0 0.0 2.5,0.0 5.0 "
      "2.5,0.0 0.0 0.0,0.0 0.0 2.5)),((0.0 10.0 2.5,0.0 0.0 0.0,0.0 5.0 "
      "2.5,0.0 10.0 2.5)),((0.0 0.0 0.0,0.0 10.0 2.5,0.0 10.0 0.0,0.0 0.0 "
      "0.0)),((0.0 10.0 2.5,5.0 10.0 2.5,0.0 10.0 0.0,0.0 10.0 2.5)),((5.0 0.0 "
      "2.5,0.0 5.0 2.5,0.0 0.0 2.5,5.0 0.0 2.5)),((10.0 10.0 2.5,0.0 10.0 "
      "0.0,5.0 10.0 2.5,10.0 10.0 2.5)),((0.0 10.0 0.0,10.0 10.0 2.5,10.0 10.0 "
      "0.0,0.0 10.0 0.0)),((10.0 10.0 2.5,10.0 5.0 2.5,10.0 10.0 0.0,10.0 10.0 "
      "2.5)),((10.0 0.0 2.5,10.0 5.0 2.5,5.0 0.0 2.5,10.0 0.0 2.5)),((10.0 5.0 "
      "2.5,5.0 10.0 2.5,5.0 0.0 2.5,10.0 5.0 2.5)),((5.0 10.0 2.5,10.0 5.0 "
      "2.5,10.0 10.0 2.5,5.0 10.0 2.5)),((10.0 0.0 2.5,10.0 10.0 0.0,10.0 5.0 "
      "2.5,10.0 0.0 2.5)),((10.0 10.0 0.0,10.0 0.0 2.5,10.0 0.0 0.0,10.0 10.0 "
      "0.0)),((5.0 10.0 2.5,0.0 5.0 2.5,5.0 0.0 2.5,5.0 10.0 2.5)),((10.0 0.0 "
      "0.0,0.0 0.0 2.5,0.0 0.0 0.0,10.0 0.0 0.0)),((10.0 0.0 2.5,5.0 0.0 "
      "2.5,10.0 0.0 0.0,10.0 0.0 2.5)),((5.0 10.0 2.5,0.0 10.0 2.5,0.0 5.0 "
      "2.5,5.0 10.0 2.5)))";

  const std::unique_ptr<Geometry> expectedGeom3 = io::readWkt(expectedGeom3Wkt);
  BOOST_CHECK(
      algorithm::covers3D(geomSplitClose->geometryN(0), *expectedGeom3));

  const std::string expectedGeom4Wkt =
      "POLYHEDRALSURFACE Z (((10.0 0.0 2.5,0.0 0.0 5.0,5.0 0.0 2.5,10.0 0.0 "
      "2.5)),((10.0 10.0 2.5,10.0 0.0 5.0,10.0 5.0 2.5,10.0 10.0 2.5)),((10.0 "
      "0.0 2.5,10.0 5.0 2.5,10.0 0.0 5.0,10.0 0.0 2.5)),((0.0 10.0 2.5,10.0 "
      "10.0 5.0,5.0 10.0 2.5,0.0 10.0 2.5)),((10.0 10.0 2.5,5.0 10.0 2.5,10.0 "
      "10.0 5.0,10.0 10.0 2.5)),((0.0 10.0 2.5,0.0 5.0 2.5,0.0 10.0 5.0,0.0 "
      "10.0 2.5)),((5.0 0.0 8.0,0.0 0.0 5.0,10.0 0.0 5.0,5.0 0.0 8.0)),((5.0 "
      "10.0 8.0,10.0 10.0 5.0,0.0 10.0 5.0,5.0 10.0 8.0)),((5.0 10.0 8.0,0.0 "
      "0.0 5.0,5.0 0.0 8.0,5.0 10.0 8.0)),((5.0 10.0 8.0,0.0 10.0 5.0,0.0 0.0 "
      "5.0,5.0 10.0 8.0)),((5.0 10.0 8.0,5.0 0.0 8.0,10.0 10.0 5.0,5.0 10.0 "
      "8.0)),((10.0 10.0 5.0,5.0 0.0 8.0,10.0 0.0 5.0,10.0 10.0 5.0)),((10.0 "
      "10.0 5.0,0.0 10.0 2.5,0.0 10.0 5.0,10.0 10.0 5.0)),((0.0 5.0 2.5,5.0 "
      "0.0 2.5,0.0 0.0 2.5,0.0 5.0 2.5)),((0.0 10.0 2.5,5.0 10.0 2.5,0.0 5.0 "
      "2.5,0.0 10.0 2.5)),((10.0 0.0 5.0,10.0 10.0 2.5,10.0 10.0 5.0,10.0 0.0 "
      "5.0)),((0.0 0.0 2.5,0.0 10.0 5.0,0.0 5.0 2.5,0.0 0.0 2.5)),((0.0 10.0 "
      "5.0,0.0 0.0 2.5,0.0 0.0 5.0,0.0 10.0 5.0)),((5.0 10.0 2.5,10.0 5.0 "
      "2.5,5.0 0.0 2.5,5.0 10.0 2.5)),((10.0 10.0 2.5,10.0 5.0 2.5,5.0 10.0 "
      "2.5,10.0 10.0 2.5)),((0.0 0.0 2.5,5.0 0.0 2.5,0.0 0.0 5.0,0.0 0.0 "
      "2.5)),((5.0 10.0 2.5,5.0 0.0 2.5,0.0 5.0 2.5,5.0 10.0 2.5)),((10.0 5.0 "
      "2.5,10.0 0.0 2.5,5.0 0.0 2.5,10.0 5.0 2.5)),((0.0 0.0 5.0,10.0 0.0 "
      "2.5,10.0 0.0 5.0,0.0 0.0 5.0)))";

  const std::unique_ptr<Geometry> expectedGeom4 = io::readWkt(expectedGeom4Wkt);
  BOOST_CHECK(
      algorithm::covers3D(geomSplitClose->geometryN(1), *expectedGeom4));
}

BOOST_AUTO_TEST_CASE(testSplit_Solid)
{
  std::unique_ptr<Geometry> houseSolid = io::readWkt(
      "SOLID Z ((((0 0 0,0 10 0,10 10 0,10 0 0,0 0 0)),((0 0 0,10 0 0,10 0 5,0 "
      "0 5,0 0 0)),((10 0 0,10 10 0,10 10 5,10 0 5,10 0 0)),((10 10 0,0 10 0,0 "
      "10 5,10 10 5,10 10 0)),((0 10 0,0 0 0,0 0 5,0 10 5,0 10 0)),((0 0 5,10 "
      "0 5,5 0 8,0 0 5)),((10 10 5,0 10 5,5 10 8,10 10 5)),((0 0 5,5 0 8,5 10 "
      "8,0 10 5,0 0 5)),((10 0 5,10 10 5,5 10 8,5 0 8,10 0 5))))");
  BOOST_CHECK(!houseSolid->isEmpty());
  BOOST_CHECK_EQUAL(houseSolid->geometryTypeId(),
                    SFCGAL::GeometryType::TYPE_SOLID);

  std::unique_ptr<GeometryCollection> geomSplitOpen = algorithm::split(
      *houseSolid, Point(0, 0, 2.5), Kernel::Vector_3(0, 0, 1), false);
  std::unique_ptr<GeometryCollection> geomSplitClose = algorithm::split(
      *houseSolid, Point(0, 0, 2.5), Kernel::Vector_3(0, 0, 1), true);

  const unsigned int nrSplit = 2;
  BOOST_CHECK_EQUAL(geomSplitOpen->numGeometries(), nrSplit);
  BOOST_CHECK_EQUAL(geomSplitClose->numGeometries(), nrSplit);
  for (unsigned int i = 0; i < nrSplit; ++i) {
    BOOST_CHECK(algorithm::covers3D(geomSplitOpen->geometryN(i),
                                    geomSplitClose->geometryN(i)));
  }

  const std::string expectedGeom1Wkt =
      "SOLID Z ((((10.0 10.0 0.0,0.0 0.0 0.0,0.0 10.0 0.0,10.0 10.0 "
      "0.0)),((10.0 10.0 0.0,10.0 0.0 0.0,0.0 0.0 0.0,10.0 10.0 0.0)),((0.0 "
      "0.0 2.5,10.0 0.0 0.0,5.0 0.0 2.5,0.0 0.0 2.5)),((0.0 0.0 2.5,0.0 5.0 "
      "2.5,0.0 0.0 0.0,0.0 0.0 2.5)),((0.0 10.0 2.5,0.0 0.0 0.0,0.0 5.0 "
      "2.5,0.0 10.0 2.5)),((0.0 0.0 0.0,0.0 10.0 2.5,0.0 10.0 0.0,0.0 0.0 "
      "0.0)),((0.0 10.0 2.5,5.0 10.0 2.5,0.0 10.0 0.0,0.0 10.0 2.5)),((5.0 0.0 "
      "2.5,0.0 5.0 2.5,0.0 0.0 2.5,5.0 0.0 2.5)),((10.0 10.0 2.5,0.0 10.0 "
      "0.0,5.0 10.0 2.5,10.0 10.0 2.5)),((0.0 10.0 0.0,10.0 10.0 2.5,10.0 10.0 "
      "0.0,0.0 10.0 0.0)),((10.0 10.0 2.5,10.0 5.0 2.5,10.0 10.0 0.0,10.0 10.0 "
      "2.5)),((10.0 0.0 2.5,10.0 5.0 2.5,5.0 0.0 2.5,10.0 0.0 2.5)),((10.0 5.0 "
      "2.5,5.0 10.0 2.5,5.0 0.0 2.5,10.0 5.0 2.5)),((5.0 10.0 2.5,10.0 5.0 "
      "2.5,10.0 10.0 2.5,5.0 10.0 2.5)),((10.0 0.0 2.5,10.0 10.0 0.0,10.0 5.0 "
      "2.5,10.0 0.0 2.5)),((10.0 10.0 0.0,10.0 0.0 2.5,10.0 0.0 0.0,10.0 10.0 "
      "0.0)),((5.0 10.0 2.5,0.0 5.0 2.5,5.0 0.0 2.5,5.0 10.0 2.5)),((10.0 0.0 "
      "0.0,0.0 0.0 2.5,0.0 0.0 0.0,10.0 0.0 0.0)),((10.0 0.0 2.5,5.0 0.0 "
      "2.5,10.0 0.0 0.0,10.0 0.0 2.5)),((5.0 10.0 2.5,0.0 10.0 2.5,0.0 5.0 "
      "2.5,5.0 10.0 2.5))))";

  const std::unique_ptr<Geometry> expectedGeom1 = io::readWkt(expectedGeom1Wkt);
  BOOST_CHECK(
      algorithm::covers3D(geomSplitClose->geometryN(0), *expectedGeom1));

  const std::string expectedGeom2Wkt =
      "SOLID Z ((((10.0 0.0 2.5,0.0 0.0 5.0,5.0 0.0 2.5,10.0 0.0 2.5)),((10.0 "
      "10.0 2.5,10.0 0.0 5.0,10.0 5.0 2.5,10.0 10.0 2.5)),((10.0 0.0 2.5,10.0 "
      "5.0 2.5,10.0 0.0 5.0,10.0 0.0 2.5)),((0.0 10.0 2.5,10.0 10.0 5.0,5.0 "
      "10.0 2.5,0.0 10.0 2.5)),((10.0 10.0 2.5,5.0 10.0 2.5,10.0 10.0 5.0,10.0 "
      "10.0 2.5)),((0.0 10.0 2.5,0.0 5.0 2.5,0.0 10.0 5.0,0.0 10.0 2.5)),((5.0 "
      "0.0 8.0,0.0 0.0 5.0,10.0 0.0 5.0,5.0 0.0 8.0)),((5.0 10.0 8.0,10.0 10.0 "
      "5.0,0.0 10.0 5.0,5.0 10.0 8.0)),((5.0 10.0 8.0,0.0 0.0 5.0,5.0 0.0 "
      "8.0,5.0 10.0 8.0)),((5.0 10.0 8.0,0.0 10.0 5.0,0.0 0.0 5.0,5.0 10.0 "
      "8.0)),((5.0 10.0 8.0,5.0 0.0 8.0,10.0 10.0 5.0,5.0 10.0 8.0)),((10.0 "
      "10.0 5.0,5.0 0.0 8.0,10.0 0.0 5.0,10.0 10.0 5.0)),((10.0 10.0 5.0,0.0 "
      "10.0 2.5,0.0 10.0 5.0,10.0 10.0 5.0)),((0.0 5.0 2.5,5.0 0.0 2.5,0.0 0.0 "
      "2.5,0.0 5.0 2.5)),((0.0 10.0 2.5,5.0 10.0 2.5,0.0 5.0 2.5,0.0 10.0 "
      "2.5)),((10.0 0.0 5.0,10.0 10.0 2.5,10.0 10.0 5.0,10.0 0.0 5.0)),((0.0 "
      "0.0 2.5,0.0 10.0 5.0,0.0 5.0 2.5,0.0 0.0 2.5)),((0.0 10.0 5.0,0.0 0.0 "
      "2.5,0.0 0.0 5.0,0.0 10.0 5.0)),((5.0 10.0 2.5,10.0 5.0 2.5,5.0 0.0 "
      "2.5,5.0 10.0 2.5)),((10.0 10.0 2.5,10.0 5.0 2.5,5.0 10.0 2.5,10.0 10.0 "
      "2.5)),((0.0 0.0 2.5,5.0 0.0 2.5,0.0 0.0 5.0,0.0 0.0 2.5)),((5.0 10.0 "
      "2.5,5.0 0.0 2.5,0.0 5.0 2.5,5.0 10.0 2.5)),((10.0 5.0 2.5,10.0 0.0 "
      "2.5,5.0 0.0 2.5,10.0 5.0 2.5)),((0.0 0.0 5.0,10.0 0.0 2.5,10.0 0.0 "
      "5.0,0.0 0.0 5.0))))";

  const std::unique_ptr<Geometry> expectedGeom2 = io::readWkt(expectedGeom2Wkt);
  BOOST_CHECK(
      algorithm::covers3D(geomSplitClose->geometryN(1), *expectedGeom2));
}

BOOST_AUTO_TEST_CASE(testSplit_NoIntersection)
{
  std::unique_ptr<Geometry> lshape = io::readWkt(
      "PolyhedralSurface Z (((0 0 0,0 5 0,1 5 0,1 1 0,5 1 0,5 0 0,0 0 0)),((0 "
      "0 2,5 0 2,5 1 2,1 1 2,1 5 2,0 5 2,0 0 2)),((0 0 0,0 0 2,0 5 2,0 5 0,0 0 "
      "0)),((0 5 0,0 5 2,1 5 2,1 5 0,0 5 0)),((1 5 0,1 5 2,1 1 2,1 1 0,1 5 "
      "0)),((1 1 0,1 1 2,5 1 2,5 1 0,1 1 0)),((5 1 0,5 1 2,5 0 2,5 0 0,5 1 "
      "0)),((5 0 0,5 0 2,0 0 2,0 0 0,5 0 0)))");

  // change the plane position and check intersection
  std::vector<std::tuple<Point, int>> tests = {
      {Point(2, 1.8, 0), 3},   // intersection
      {Point(3.2, 3.3, 0), 0}, // bbox intersect, no intersection
      {Point(5.2, 4.8, 0), 0}  // plane and geometry bbox do not intersect
  };

  for (const auto &[position, nrGeoms] : tests) {
    std::unique_ptr<GeometryCollection> geomSplit =
        algorithm::split(*lshape, position, Kernel::Vector_3(1, 1, 0), false);
    BOOST_CHECK_EQUAL(geomSplit->numGeometries(), nrGeoms);
  }
}

BOOST_AUTO_TEST_SUITE_END()
