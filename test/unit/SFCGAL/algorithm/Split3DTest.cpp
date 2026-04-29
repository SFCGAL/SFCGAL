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
#include "SFCGAL/algorithm/split3D.h"
#include "SFCGAL/io/wkt.h"
#include "SFCGAL/primitive3d/Cube.h"
#include "SFCGAL/version.h"

using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_Split3D)

BOOST_AUTO_TEST_CASE(testSplit_Empty)
{
  const PolyhedralSurface emptyGeom;
  BOOST_CHECK(emptyGeom.isEmpty());
  std::unique_ptr<GeometryCollection> geomSplit = algorithm::split3D(
      emptyGeom, Point(0, 0, 0), Kernel::Vector_3(1, 0, -1), false);
  BOOST_CHECK(geomSplit->isEmpty());
}

BOOST_AUTO_TEST_CASE(testSplit_Cube)
{
  Cube              cube(57);
  PolyhedralSurface phs = cube.generatePolyhedralSurface();
  BOOST_CHECK(!phs.isEmpty());

  std::unique_ptr<GeometryCollection> phsSplit = algorithm::split3D(
      phs, Point(0, 0, 0), Kernel::Vector_3(1, 0, -1), false);
  const unsigned int nrSplit = 2;
  BOOST_CHECK_EQUAL(phsSplit->numGeometries(), nrSplit);

  for (unsigned int i = 0; i < nrSplit; ++i) {
    std::string componentData(SFCGAL_TEST_DIRECTORY);
    componentData += "/data/split3D/cubeComponent" + std::to_string(i) + ".wkt";
    std::ifstream ifs(componentData.c_str());
    BOOST_REQUIRE(ifs.good());

    std::string expectedWkt;
    std::getline(ifs, expectedWkt);

    const std::unique_ptr<Geometry> expectedGeom = io::readWkt(expectedWkt);
    BOOST_CHECK(algorithm::covers3D(phsSplit->geometryN(i), *expectedGeom));
  }
}

BOOST_AUTO_TEST_CASE(testSplit_UShape)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/split3D/ushape.wkt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string inputWkt;
  std::getline(ifs, inputWkt);
  const std::unique_ptr<Geometry> inputGeom(io::readWkt(inputWkt));
  BOOST_CHECK(!inputGeom->isEmpty());
  BOOST_CHECK_EQUAL(inputGeom->geometryTypeId(),
                    SFCGAL::GeometryType::TYPE_POLYHEDRALSURFACE);

  std::unique_ptr<GeometryCollection> phsSplit = algorithm::split3D(
      *inputGeom, Point(1, 6, 0), Kernel::Vector_3(0, 1, 0), false);
  const unsigned int nrSplit = 3;
  BOOST_CHECK_EQUAL(phsSplit->numGeometries(), nrSplit);

  for (unsigned int i = 0; i < nrSplit; ++i) {
    std::string componentData(SFCGAL_TEST_DIRECTORY);

#if SFCGAL_CGAL_VERSION_MAJOR >= 6 and SFCGAL_CGAL_VERSION_MINOR >= 1
    const std::string cgal_version = "61";
#elif SFCGAL_CGAL_VERSION_MAJOR >= 6 and SFCGAL_CGAL_VERSION_MINOR == 0
    const std::string cgal_version = "60";
#else
    const std::string cgal_version = "57";
#endif

    componentData += "/data/split3D/ushapeComponent" + std::to_string(i) +
                     "_cgal" + cgal_version + ".wkt";
    std::ifstream ifs(componentData.c_str());
    BOOST_REQUIRE(ifs.good());

    std::string expectedWkt;
    std::getline(ifs, expectedWkt);

    const std::unique_ptr<Geometry> expectedGeom = io::readWkt(expectedWkt);
    BOOST_CHECK(algorithm::covers3D(phsSplit->geometryN(i), *expectedGeom));
  }
}

BOOST_AUTO_TEST_CASE(testSplit_Tin)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/split3D/tin.wkt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string inputWkt;
  std::getline(ifs, inputWkt);
  const std::unique_ptr<Geometry> inputGeom(io::readWkt(inputWkt));
  BOOST_CHECK(!inputGeom->isEmpty());
  BOOST_CHECK_EQUAL(inputGeom->geometryTypeId(),
                    SFCGAL::GeometryType::TYPE_TRIANGULATEDSURFACE);

  std::unique_ptr<GeometryCollection> tinSplit = algorithm::split3D(
      *inputGeom, Point(2, 2, 1.2), Kernel::Vector_3(0, 0, 1), false);
  const unsigned int nrSplit = 2;
  BOOST_CHECK_EQUAL(tinSplit->numGeometries(), nrSplit);

  for (unsigned int i = 0; i < nrSplit; ++i) {
    std::string componentData(SFCGAL_TEST_DIRECTORY);
    componentData += "/data/split3D/tinComponent" + std::to_string(i) + ".wkt";
    std::ifstream ifs(componentData.c_str());
    BOOST_REQUIRE(ifs.good());

    std::string expectedWkt;
    std::getline(ifs, expectedWkt);

    const std::unique_ptr<Geometry> expectedGeom = io::readWkt(expectedWkt);
    BOOST_CHECK(algorithm::covers3D(tinSplit->geometryN(i), *expectedGeom));
  }
}

BOOST_AUTO_TEST_CASE(testSplit_Solid)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/split3D/houseSolid.wkt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string inputWkt;
  std::getline(ifs, inputWkt);
  const std::unique_ptr<Geometry> houseSolid(io::readWkt(inputWkt));
  BOOST_CHECK(!houseSolid->isEmpty());
  BOOST_CHECK_EQUAL(houseSolid->geometryTypeId(),
                    SFCGAL::GeometryType::TYPE_SOLID);

  std::unique_ptr<GeometryCollection> geomSplitOpen = algorithm::split3D(
      *houseSolid, Point(0, 0, 2.5), Kernel::Vector_3(0, 0, 1), false);
  std::unique_ptr<GeometryCollection> geomSplitClose = algorithm::split3D(
      *houseSolid, Point(0, 0, 2.5), Kernel::Vector_3(0, 0, 1), true);

  const unsigned int nrSplit = 2;
  BOOST_CHECK_EQUAL(geomSplitOpen->numGeometries(), nrSplit);
  BOOST_CHECK_EQUAL(geomSplitClose->numGeometries(), nrSplit);
  for (unsigned int i = 0; i < nrSplit; ++i) {
    BOOST_CHECK(algorithm::covers3D(geomSplitOpen->geometryN(i),
                                    geomSplitClose->geometryN(i)));
  }

  for (unsigned int i = 0; i < nrSplit; ++i) {
    std::string componentData(SFCGAL_TEST_DIRECTORY);
    componentData +=
        "/data/split3D/houseSolidComponent" + std::to_string(i) + ".wkt";
    std::ifstream ifs(componentData.c_str());
    BOOST_REQUIRE(ifs.good());

    std::string expectedWkt;
    std::getline(ifs, expectedWkt);

    const std::unique_ptr<Geometry> expectedGeom = io::readWkt(expectedWkt);
    BOOST_CHECK(
        algorithm::covers3D(geomSplitClose->geometryN(i), *expectedGeom));
  }
}

BOOST_AUTO_TEST_CASE(testSplit_NoIntersection)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/split3D/lshape.wkt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string inputWkt;
  std::getline(ifs, inputWkt);
  const std::unique_ptr<Geometry> lshape(io::readWkt(inputWkt));

  // change the plane position and check intersection
  std::vector<std::tuple<Point, int>> tests = {
      {Point(2, 1.8, 0), 3},   // intersection
      {Point(3.2, 3.3, 0), 1}, // bbox intersect, no intersection
      {Point(5.2, 4.8, 0), 1}  // plane and geometry bbox do not intersect
  };

  for (const auto &[position, nrGeoms] : tests) {
    std::unique_ptr<GeometryCollection> geomSplit =
        algorithm::split3D(*lshape, position, Kernel::Vector_3(1, 1, 0), false);
    BOOST_CHECK_EQUAL(geomSplit->numGeometries(), nrGeoms);
  }
}

BOOST_AUTO_TEST_SUITE_END()
