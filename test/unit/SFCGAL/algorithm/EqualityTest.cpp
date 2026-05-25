// Copyright (c) 2024-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <memory>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/equality.h"
#include "SFCGAL/io/wkt.h"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;
using namespace SFCGAL::algorithm;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_EqualityTest)

// --- EqualityStrictness flag manipulation ---

BOOST_AUTO_TEST_CASE(testDefaultStrictness)
{
  EqualityStrictness strictness;
  BOOST_CHECK(!(strictness & EqualityStrictness::CheckCoverOrPoint));
  BOOST_CHECK(!(strictness & EqualityStrictness::SubGeomOrdered));
  BOOST_CHECK(!(strictness & EqualityStrictness::SubPartOrdered));
  BOOST_CHECK(!(strictness & EqualityStrictness::InternalPointOrdered));
  BOOST_CHECK_EQUAL(strictness.toString(), "CheckPoint");
}

BOOST_AUTO_TEST_CASE(testFlagCombination)
{
  auto strictness = EqualityStrictness(EqualityStrictness::CheckCoverOrPoint) |
                    EqualityStrictness::SubGeomOrdered;
  BOOST_CHECK(strictness & EqualityStrictness::CheckCoverOrPoint);
  BOOST_CHECK(strictness & EqualityStrictness::SubGeomOrdered);
}

BOOST_AUTO_TEST_CASE(testConflictCheckCoverWithInternalPoint)
{
  auto strictness = EqualityStrictness(EqualityStrictness::CheckCoverOrPoint);
  BOOST_CHECK_THROW(strictness | EqualityStrictness::InternalPointOrdered,
                    SFCGAL::Exception);
}

BOOST_AUTO_TEST_CASE(testConflictMultiInternalPoint)
{
  auto strictness =
      EqualityStrictness(EqualityStrictness::InternalPointOrdered);
  BOOST_CHECK_THROW(strictness | EqualityStrictness::InternalPointShifted,
                    SFCGAL::Exception);
}

BOOST_AUTO_TEST_CASE(testConflictInternalPointThenCheckCover)
{
  auto strictness =
      EqualityStrictness(EqualityStrictness::InternalPointOrdered);
  BOOST_CHECK_THROW(strictness | EqualityStrictness::CheckCoverOrPoint,
                    SFCGAL::Exception);
}

BOOST_AUTO_TEST_CASE(testAssignFlag)
{
  EqualityStrictness strictness;
  strictness = EqualityStrictness::SubGeomOrdered;
  BOOST_CHECK(strictness & EqualityStrictness::SubGeomOrdered);
  BOOST_CHECK(!(strictness & EqualityStrictness::InternalPointOrdered));
}

BOOST_AUTO_TEST_CASE(testAssignEqualityStrictness)
{
  EqualityStrictness a(EqualityStrictness::SubPartOrdered);
  EqualityStrictness b;
  b = a;
  BOOST_CHECK(b & EqualityStrictness::SubPartOrdered);
}

BOOST_AUTO_TEST_CASE(testToString)
{
  EqualityStrictness strictness =
      EqualityStrictness() | EqualityStrictness::SubGeomOrdered;
  BOOST_CHECK_EQUAL(strictness.toString(), "CheckPoint | SubGeomOrdered");

  auto strictness2 = EqualityStrictness(EqualityStrictness::CheckCoverOrPoint) |
                     EqualityStrictness::SubGeomOrdered;
  BOOST_CHECK_EQUAL(strictness2.toString(), "CheckCover | SubGeomOrdered");

  auto strictness3 =
      EqualityStrictness(EqualityStrictness::InternalPointShifted);
  BOOST_CHECK_EQUAL(strictness3.toString(),
                    "CheckPoint | InternalPointShifted");
}

// --- almostEqual: basic point comparisons ---

BOOST_AUTO_TEST_CASE(testEqualPoints)
{
  auto gA = io::readWkt("POINT (1 2 3)");
  auto gB = io::readWkt("POINT (1 2 3)");
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0));
}

BOOST_AUTO_TEST_CASE(testUnequalPoints)
{
  auto gA = io::readWkt("POINT (1 2 3)");
  auto gB = io::readWkt("POINT (4 5 6)");
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, 0.0));
}

BOOST_AUTO_TEST_CASE(testEqualPointsWithTolerance)
{
  auto gA = io::readWkt("POINT (1.0 2.0)");
  auto gB = io::readWkt("POINT (1.0001 2.0001)");
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 1e-3));
}

BOOST_AUTO_TEST_CASE(testUnequalPointsTypeMismatch)
{
  auto gA = io::readWkt("POINT (1 2)");
  auto gB = io::readWkt("LINESTRING (1 2, 3 4)");
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, 0.0));
}

BOOST_AUTO_TEST_CASE(testEmptyPoints)
{
  auto gA = io::readWkt("POINT EMPTY");
  auto gB = io::readWkt("POINT EMPTY");
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0));
}

BOOST_AUTO_TEST_CASE(testOneEmptyOneNonEmpty)
{
  auto gA = io::readWkt("POINT EMPTY");
  auto gB = io::readWkt("POINT (1 2)");
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, 0.0));
}

BOOST_AUTO_TEST_CASE(testNonEmptyOneEmpty)
{
  auto gA = io::readWkt("POINT (1 2)");
  auto gB = io::readWkt("POINT EMPTY");
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, 0.0));
}

// --- almostEqual: CheckCoverOrPoint mode ---

BOOST_AUTO_TEST_CASE(testCheckCoverEqual)
{
  auto gA         = io::readWkt("POLYGON ((0 0, 0 1, 1 1, 1 0, 0 0))");
  auto gB         = io::readWkt("POLYGON ((0 0, 0 1, 1 1, 1 0, 0 0))");
  auto strictness = EqualityStrictness(EqualityStrictness::CheckCoverOrPoint);
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

// --- almostEqual: Multi* SubGeomOrdered ---

BOOST_AUTO_TEST_CASE(testMultiPointOrdered)
{
  auto gA         = io::readWkt("MULTIPOINT ((0 0), (1 1))");
  auto gB         = io::readWkt("MULTIPOINT ((0 0), (1 1))");
  auto strictness = EqualityStrictness(EqualityStrictness::SubGeomOrdered);
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testMultiPointNonOrdered)
{
  auto gA         = io::readWkt("MULTIPOINT ((0 0), (1 1))");
  auto gB         = io::readWkt("MULTIPOINT ((1 1), (0 0))");
  auto strictness = EqualityStrictness::pointNonOrdered();
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testMultiPointOrderedMismatch)
{
  auto gA         = io::readWkt("MULTIPOINT ((0 0), (1 1))");
  auto gB         = io::readWkt("MULTIPOINT ((1 1), (0 0))");
  auto strictness = EqualityStrictness(EqualityStrictness::SubGeomOrdered);
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testMultiLineStringNonOrdered)
{
  auto gA         = io::readWkt("MULTILINESTRING ((0 0, 1 1), (2 2, 3 3))");
  auto gB         = io::readWkt("MULTILINESTRING ((2 2, 3 3), (0 0, 1 1))");
  auto strictness = EqualityStrictness::pointNonOrdered();
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testMultiPolygonNonOrdered)
{
  auto gA = io::readWkt(
      "MULTIPOLYGON (((0 0, 0 1, 1 1, 1 0, 0 0)),((2 2, 2 3, 3 3, 3 2, 2 2)))");
  auto gB = io::readWkt(
      "MULTIPOLYGON (((2 2, 2 3, 3 3, 3 2, 2 2)),((0 0, 0 1, 1 1, 1 0, 0 0)))");
  auto strictness = EqualityStrictness::pointNonOrdered();
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testMultiPolygonOrderedMismatch)
{
  auto gA = io::readWkt(
      "MULTIPOLYGON (((0 0, 0 1, 1 1, 1 0, 0 0)),((2 2, 2 3, 3 3, 3 2, 2 2)))");
  auto gB = io::readWkt(
      "MULTIPOLYGON (((2 2, 2 3, 3 3, 3 2, 2 2)),((0 0, 0 1, 1 1, 1 0, 0 0)))");
  auto strictness = EqualityStrictness(EqualityStrictness::SubGeomOrdered);
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testMultiGeomCountMismatch)
{
  auto gA = io::readWkt("MULTIPOINT ((0 0), (1 1))");
  auto gB = io::readWkt("MULTIPOINT ((0 0))");
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, 0.0));
}

// --- almostEqual: SubPartOrdered ---

BOOST_AUTO_TEST_CASE(testPolygonSubPartOrdered)
{
  auto gA = io::readWkt(
      "POLYGON ((0 0, 0 10, 10 10, 10 0, 0 0), (2 2, 2 4, 4 4, 4 2, 2 2))");
  auto gB = io::readWkt(
      "POLYGON ((0 0, 0 10, 10 10, 10 0, 0 0), (2 2, 2 4, 4 4, 4 2, 2 2))");
  auto strictness = EqualityStrictness(EqualityStrictness::SubPartOrdered);
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testPolygonSubPartNonOrdered)
{
  auto gA = io::readWkt(
      "POLYGON ((0 0, 0 10, 10 10, 10 0, 0 0), (2 2, 2 4, 4 4, 4 2, 2 2))");
  auto gB = io::readWkt(
      "POLYGON ((0 0, 0 10, 10 10, 10 0, 0 0), (3 3, 3 5, 5 5, 5 3, 3 3))");
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, 0.0));
}

// --- almostEqual: InternalPointShifted ---

BOOST_AUTO_TEST_CASE(testLineStringShifted)
{
  auto gA = io::readWkt("LINESTRING (1 0, 2 0, 3 0, 4 0)");
  auto gB = io::readWkt("LINESTRING (3 0, 4 0, 1 0, 2 0)");
  auto strictness =
      EqualityStrictness(EqualityStrictness::InternalPointShifted);
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testLineStringShiftedClosed)
{
  auto gA = io::readWkt("LINESTRING (0 0, 1 0, 1 1, 0 1, 0 0)");
  auto gB = io::readWkt("LINESTRING (1 0, 1 1, 0 1, 0 0, 1 0)");
  auto strictness =
      EqualityStrictness(EqualityStrictness::InternalPointShifted);
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testLineStringShiftedMismatch)
{
  auto gA = io::readWkt("LINESTRING (1 0, 2 0, 3 0, 4 0)");
  auto gB = io::readWkt("LINESTRING (5 0, 6 0, 7 0, 8 0)");
  auto strictness =
      EqualityStrictness(EqualityStrictness::InternalPointShifted);
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testLineStringShiftedEmpty)
{
  auto gA = io::readWkt("LINESTRING EMPTY");
  auto gB = io::readWkt("LINESTRING EMPTY");
  auto strictness =
      EqualityStrictness(EqualityStrictness::InternalPointShifted);
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testLineStringShiftedSizeMismatch)
{
  auto gA = io::readWkt("LINESTRING (1 0, 2 0, 3 0)");
  auto gB = io::readWkt("LINESTRING (1 0, 2 0, 3 0, 4 0)");
  auto strictness =
      EqualityStrictness(EqualityStrictness::InternalPointShifted);
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

// --- almostEqual: InternalPointInverted ---

BOOST_AUTO_TEST_CASE(testLineStringInverted)
{
  auto gA = io::readWkt("LINESTRING (1 0, 2 0, 3 0, 4 0)");
  auto gB = io::readWkt("LINESTRING (4 0, 3 0, 2 0, 1 0)");
  auto strictness =
      EqualityStrictness(EqualityStrictness::InternalPointInverted);
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testLineStringInvertedAlsoMatchesOrdered)
{
  auto gA = io::readWkt("LINESTRING (1 0, 2 0, 3 0, 4 0)");
  auto gB = io::readWkt("LINESTRING (1 0, 2 0, 3 0, 4 0)");
  auto strictness =
      EqualityStrictness(EqualityStrictness::InternalPointInverted);
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testLineStringInvertedMismatch)
{
  auto gA = io::readWkt("LINESTRING (1 0, 2 0, 3 0, 4 0)");
  auto gB = io::readWkt("LINESTRING (5 0, 6 0, 7 0, 8 0)");
  auto strictness =
      EqualityStrictness(EqualityStrictness::InternalPointInverted);
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

// --- almostEqual: InternalPointOrdered ---

BOOST_AUTO_TEST_CASE(testLineStringOrdered)
{
  auto gA = io::readWkt("LINESTRING (1 0, 2 0, 3 0, 4 0)");
  auto gB = io::readWkt("LINESTRING (1 0, 2 0, 3 0, 4 0)");
  auto strictness =
      EqualityStrictness(EqualityStrictness::InternalPointOrdered);
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

BOOST_AUTO_TEST_CASE(testLineStringOrderedMismatch)
{
  auto gA = io::readWkt("LINESTRING (1 0, 2 0, 3 0, 4 0)");
  auto gB = io::readWkt("LINESTRING (4 0, 3 0, 2 0, 1 0)");
  auto strictness =
      EqualityStrictness(EqualityStrictness::InternalPointOrdered);
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, 0.0, strictness));
}

// --- almostEqual: TriangulatedSurface/PolyhedralSurface ---

BOOST_AUTO_TEST_CASE(testTriangulatedSurfaceEqual)
{
  auto gA = io::readWkt(
      "TIN (((0 0 0, 1 0 0, 0 1 0, 0 0 0)), ((0 0 0, 0 1 0, 0 0 1, 0 0 0)))");
  auto gB = io::readWkt(
      "TIN (((0 0 0, 1 0 0, 0 1 0, 0 0 0)), ((0 0 0, 0 1 0, 0 0 1, 0 0 0)))");
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0));
}

BOOST_AUTO_TEST_CASE(testPolyhedralSurfaceEqual)
{
  auto gA = io::readWkt(
      "POLYHEDRALSURFACE (((0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 0)), ((0 0 0, 0 0 "
      "1, 1 0 1, 1 0 0, 0 0 0)))");
  auto gB = io::readWkt(
      "POLYHEDRALSURFACE (((0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 0)), ((0 0 0, 0 0 "
      "1, 1 0 1, 1 0 0, 0 0 0)))");
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0));
}

// --- almostEqual: Solid ---

BOOST_AUTO_TEST_CASE(testSolidEqual)
{
  auto gA = io::readWkt(
      "SOLID ((((0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 0)), ((0 0 0, 0 0 1, 1 0 1, "
      "1 0 0, 0 0 0)), ((0 0 0, 0 1 0, 0 1 1, 0 0 1, 0 0 0)), ((1 1 1, 1 0 1, "
      "0 0 1, 0 1 1, 1 1 1)), ((1 1 1, 0 1 1, 0 1 0, 1 1 0, 1 1 1)), ((1 1 1, "
      "1 1 0, 1 0 0, 1 0 1, 1 1 1))))");
  auto gB = io::readWkt(
      "SOLID ((((0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 0)), ((0 0 0, 0 0 1, 1 0 1, "
      "1 0 0, 0 0 0)), ((0 0 0, 0 1 0, 0 1 1, 0 0 1, 0 0 0)), ((1 1 1, 1 0 1, "
      "0 0 1, 0 1 1, 1 1 1)), ((1 1 1, 0 1 1, 0 1 0, 1 1 0, 1 1 1)), ((1 1 1, "
      "1 1 0, 1 0 0, 1 0 1, 1 1 1))))");
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0));
}

// --- almostEqual: GeometryCollection ---

BOOST_AUTO_TEST_CASE(testGeometryCollectionEqual)
{
  auto gA =
      io::readWkt("GEOMETRYCOLLECTION (POINT (0 0), LINESTRING (0 0, 1 1))");
  auto gB =
      io::readWkt("GEOMETRYCOLLECTION (POINT (0 0), LINESTRING (0 0, 1 1))");
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, 0.0));
}

// --- almostEqual: negative tolerance (exact) ---

BOOST_AUTO_TEST_CASE(testNegativeToleranceExactEqual)
{
  auto gA = io::readWkt("POINT (1 2 3)");
  auto gB = io::readWkt("POINT (1 2 3)");
  BOOST_CHECK(algorithm::almostEqual(*gA, *gB, -1.0));
}

BOOST_AUTO_TEST_CASE(testNegativeToleranceExactUnequal)
{
  auto gA = io::readWkt("POINT (1 2 3)");
  auto gB = io::readWkt("POINT (1.001 2.001 3.001)");
  BOOST_CHECK(!algorithm::almostEqual(*gA, *gB, -1.0));
}

// --- allPointOrdered static ---

BOOST_AUTO_TEST_CASE(testAllPointOrdered)
{
  auto strictness = EqualityStrictness::allPointOrdered();
  BOOST_CHECK(strictness & EqualityStrictness::SubGeomOrdered);
  BOOST_CHECK(strictness & EqualityStrictness::SubPartOrdered);
  BOOST_CHECK(strictness & EqualityStrictness::InternalPointOrdered);
}

BOOST_AUTO_TEST_CASE(testPointNonOrdered)
{
  auto strictness = EqualityStrictness::pointNonOrdered();
  BOOST_CHECK(!(strictness & EqualityStrictness::SubGeomOrdered));
  BOOST_CHECK(!(strictness & EqualityStrictness::InternalPointOrdered));
}

BOOST_AUTO_TEST_CASE(testCoverSubGeomNonOrdered)
{
  auto strictness = EqualityStrictness::coverSubGeomNonOrdered();
  BOOST_CHECK(strictness & EqualityStrictness::CheckCoverOrPoint);
  BOOST_CHECK(!(strictness & EqualityStrictness::SubGeomOrdered));
}

BOOST_AUTO_TEST_CASE(testCoverSubGeomOrdered)
{
  auto strictness = EqualityStrictness::coverSubGeomOrdered();
  BOOST_CHECK(strictness & EqualityStrictness::CheckCoverOrPoint);
  BOOST_CHECK(strictness & EqualityStrictness::SubGeomOrdered);
}

BOOST_AUTO_TEST_SUITE_END()
