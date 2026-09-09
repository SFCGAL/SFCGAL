// Copyright (c) 2025-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

/**
 * @file test_minimal.cpp
 * @brief Minimal unit tests for SFCGALOP - one test per main functionality
 */

/// @brief Test module definition for Boost.Test framework
#define BOOST_TEST_MODULE SFCGALOP_Minimal_Tests
#include <boost/test/unit_test.hpp>

#include "../error_handler.hpp"
#include "../io.hpp"
#include "../operations/operations.hpp"

#include <SFCGAL/LineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>

#include <sstream>

/// @brief Test geometry loading from WKT format
BOOST_AUTO_TEST_CASE(test_load_wkt)
{
  std::string wkt  = "POINT(1.5 2.5)";
  auto        geom = load_geometry(wkt);
  BOOST_CHECK(geom != nullptr);
}

/// @brief Test geometry validation functionality
BOOST_AUTO_TEST_CASE(test_validate)
{
  SFCGAL::Polygon    polygon;
  SFCGAL::LineString ring;
  ring.addPoint(SFCGAL::Point(0, 0));
  ring.addPoint(SFCGAL::Point(0, 1));
  ring.addPoint(SFCGAL::Point(1, 1));
  ring.addPoint(SFCGAL::Point(1, 0));
  ring.addPoint(SFCGAL::Point(0, 0));
  polygon.setExteriorRing(ring);

  auto result = ErrorHandler::validate_geometry(polygon);
  BOOST_CHECK(result.valid);
}

/// @brief Test area calculation operation
BOOST_AUTO_TEST_CASE(test_area)
{
  SFCGAL::Polygon    polygon;
  SFCGAL::LineString ring;
  ring.addPoint(SFCGAL::Point(0, 0));
  ring.addPoint(SFCGAL::Point(0, 10));
  ring.addPoint(SFCGAL::Point(10, 10));
  ring.addPoint(SFCGAL::Point(10, 0));
  ring.addPoint(SFCGAL::Point(0, 0));
  polygon.setExteriorRing(ring);

  auto result = Operations::execute_operation("area", "", &polygon, nullptr);
  BOOST_CHECK(result.has_value());
}

/// @brief Test distance calculation operation
BOOST_AUTO_TEST_CASE(test_distance)
{
  SFCGAL::Point point1(0, 0);
  SFCGAL::Point point2(3, 4);

  auto result = Operations::execute_operation("distance", "", &point1, &point2);
  BOOST_CHECK(result.has_value());
}

/// @brief Test operations list retrieval
BOOST_AUTO_TEST_CASE(test_operations_list)
{
  auto ops = Operations::get_all_operations_info();
  BOOST_CHECK(!ops.empty());
}

/// @brief Test exception handling mechanism
BOOST_AUTO_TEST_CASE(test_exception)
{
  ErrorHandler::SfcgalopException ex("Test",
                                     ErrorHandler::ErrorCode::INVALID_GEOMETRY);
  BOOST_CHECK_EQUAL(std::string(ex.what()), "Test");
}

/// @brief Test WKT output formatting
BOOST_AUTO_TEST_CASE(test_output_wkt)
{
  SFCGAL::Point     point(1, 2);
  std::stringstream result;
  auto             *old_cout = std::cout.rdbuf(result.rdbuf());
  IO::print_result(std::make_optional(1.0), OutputFormat::WKT, 6);
  std::cout.rdbuf(old_cout);
  BOOST_CHECK(!result.str().empty());
}

/// @brief A precision of -1 asks for the exact form
BOOST_AUTO_TEST_CASE(test_output_precision_exact)
{
  auto geom = load_geometry("POINT (0.1 0.2)");
  BOOST_REQUIRE(geom != nullptr);

  std::stringstream exact;
  IO::print_result(std::make_optional(geom->clone()), OutputFormat::WKT, -1,
                   exact);
  BOOST_CHECK_EQUAL(exact.str(), "POINT (1/10 1/5)\n");

  std::stringstream rounded;
  IO::print_result(std::make_optional(geom->clone()), OutputFormat::WKT, 6,
                   rounded);
  BOOST_CHECK_EQUAL(rounded.str(), "POINT (0.100000 0.200000)\n");
}

/// @brief Test geometric intersection operation
BOOST_AUTO_TEST_CASE(test_intersection)
{
  std::string wkt1 = "POLYGON((0 0, 0 4, 4 4, 4 0, 0 0))";
  std::string wkt2 = "POLYGON((2 2, 2 6, 6 6, 6 2, 2 2))";

  auto geom1 = load_geometry(wkt1);
  auto geom2 = load_geometry(wkt2);
  BOOST_CHECK(geom1 != nullptr);
  BOOST_CHECK(geom2 != nullptr);

  auto result = Operations::execute_operation("intersection", "", geom1.get(),
                                              geom2.get());
  BOOST_CHECK(result.has_value());
}

/// @brief Test convex hull calculation
BOOST_AUTO_TEST_CASE(test_convexhull)
{
  std::string wkt  = "MULTIPOINT((0 0),(1 1),(1 0),(0 1))";
  auto        geom = load_geometry(wkt);
  BOOST_CHECK(geom != nullptr);

  auto result =
      Operations::execute_operation("convexhull", "", geom.get(), nullptr);
  BOOST_CHECK(result.has_value());
}

/// @brief Test centroid calculation
BOOST_AUTO_TEST_CASE(test_centroid)
{
  std::string wkt  = "MULTIPOINT((0 0),(1 1),(1 0),(0 1))";
  auto        geom = load_geometry(wkt);
  BOOST_CHECK(geom != nullptr);

  auto result =
      Operations::execute_operation("centroid", "", geom.get(), nullptr);
  BOOST_CHECK(result.has_value());
}

/// @brief Test centroid calculation
BOOST_AUTO_TEST_CASE(test_centroid_3d)
{
  std::string wkt  = "MULTIPOINT((0 0 0),(1 1 0),(1 0 1),(0 1 1))";
  auto        geom = load_geometry(wkt);
  BOOST_CHECK(geom != nullptr);

  auto result =
      Operations::execute_operation("centroid_3d", "", geom.get(), nullptr);
  BOOST_CHECK(result.has_value());
}

/// @brief Test geometry validity checking
BOOST_AUTO_TEST_CASE(test_is_valid)
{
  SFCGAL::Point point(1, 2);
  auto result = Operations::execute_operation("is_valid", "", &point, nullptr);
  BOOST_CHECK(result.has_value());
}

/// @brief A patch whose vertices leave the plane by more than the default
/// tolerance is only accepted once tolerance is raised above that distance.
BOOST_AUTO_TEST_CASE(test_is_valid_tolerance)
{
  auto geom = load_geometry(
      "POLYHEDRALSURFACE Z (((0 0 0,1 0 0,1 1 0.0000001,0 1 0,0 0 0)))");
  BOOST_REQUIRE(geom != nullptr);

  auto tooStrict = Operations::execute_operation("is_valid", "tolerance=1e-8",
                                                 geom.get(), nullptr);

  if (!tooStrict.has_value()) {
    BOOST_FAIL("is_valid returned an empty result for a strict tolerance");
    return;
  }

  BOOST_CHECK_EQUAL(std::get<bool>(tooStrict.value()), false);

  auto loose = Operations::execute_operation("is_valid", "tolerance=1e-6",
                                             geom.get(), nullptr);

  if (!loose.has_value()) {
    BOOST_FAIL("is_valid returned an empty result for a loose tolerance");
    return;
  }

  BOOST_CHECK_EQUAL(std::get<bool>(loose.value()), true);

  // An unparsable value is silently converted to 0.0 by parse_double, which
  // is stricter than the default and still rejects this geometry.
  auto garbage = Operations::execute_operation("is_valid", "tolerance=abc",
                                               geom.get(), nullptr);

  if (!garbage.has_value()) {
    BOOST_FAIL("is_valid returned an empty result for an invalid tolerance");
    return;
  }

  BOOST_CHECK_EQUAL(std::get<bool>(garbage.value()), false);
}

/// @brief is_simple applies the same tolerance to the planarity of its patches.
BOOST_AUTO_TEST_CASE(test_is_simple_tolerance)
{
  auto geom = load_geometry(
      "POLYHEDRALSURFACE Z (((0 0 0,1 0 0,1 1 0.0000001,0 1 0,0 0 0)))");
  BOOST_REQUIRE(geom != nullptr);

  auto tooStrict = Operations::execute_operation("is_simple", "tolerance=1e-8",
                                                 geom.get(), nullptr);

  if (!tooStrict.has_value()) {
    BOOST_FAIL("is_simple returned an empty result for a strict tolerance");
    return;
  }

  BOOST_CHECK_EQUAL(std::get<bool>(tooStrict.value()), false);

  auto loose = Operations::execute_operation("is_simple", "tolerance=1e-6",
                                             geom.get(), nullptr);

  if (!loose.has_value()) {
    BOOST_FAIL("is_simple returned an empty result for a loose tolerance");
    return;
  }

  BOOST_CHECK_EQUAL(std::get<bool>(loose.value()), true);
}

/// @brief straight_skeleton drops the segments shorter than the tolerance.
BOOST_AUTO_TEST_CASE(test_straight_skeleton_tolerance)
{
  auto geom = load_geometry("POLYGON ((0 0,10 0,10 4,0 4,0 0))");
  BOOST_REQUIRE(geom != nullptr);

  auto kept = Operations::execute_operation("straight_skeleton", "tolerance=2",
                                            geom.get(), nullptr);

  if (!kept.has_value()) {
    BOOST_FAIL("straight_skeleton returned an empty result for tolerance=2");
    return;
  }

  const auto &keptGeometry =
      std::get<std::unique_ptr<SFCGAL::Geometry>>(kept.value());

  if (keptGeometry == nullptr) {
    BOOST_FAIL("straight_skeleton returned a null geometry for tolerance=2");
    return;
  }

  BOOST_CHECK_EQUAL(keptGeometry->numGeometries(), 5U);

  auto filtered = Operations::execute_operation(
      "straight_skeleton", "tolerance=3", geom.get(), nullptr);

  if (!filtered.has_value()) {
    BOOST_FAIL("straight_skeleton returned an empty result for tolerance=3");
    return;
  }

  const auto &filteredGeometry =
      std::get<std::unique_ptr<SFCGAL::Geometry>>(filtered.value());

  if (filteredGeometry == nullptr) {
    BOOST_FAIL("straight_skeleton returned a null geometry for tolerance=3");
    return;
  }

  BOOST_CHECK_EQUAL(filteredGeometry->numGeometries(), 1U);
}
