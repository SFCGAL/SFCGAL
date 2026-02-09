// Copyright (c) 2025-2026, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/LineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/algorithm/sweep.h>
#include <SFCGAL/io/wkt.h>

#include <boost/test/unit_test.hpp>

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_SweepTest)

/**
 * @brief Test default anchor point behavior (backward compatibility)
 * Default anchor (0,0) should position profile origin on path
 */
BOOST_AUTO_TEST_CASE(testSweep_AnchorDefault)
{
  auto path = algorithm::create_circular_profile(10.0, 32); // Circle for path
  // Make it a helix
  for (size_t i = 0; i < path->numPoints(); ++i) {
    Point &pt = path->pointN(i);
    pt.z()    = i * 0.1;
  }

  auto profile = algorithm::create_rectangular_profile(2.0, 1.0);

  // Sweep with default options (anchor = 0,0)
  algorithm::SweepOptions options;
  auto                    result = algorithm::sweep(*path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_EQUAL(result->numGeometries(), path->numPoints() - 1);
}

/**
 * @brief Test custom anchor point (centered on profile)
 */
BOOST_AUTO_TEST_CASE(testSweep_AnchorCenter)
{
  auto path = algorithm::create_circular_profile(10.0, 32);
  for (size_t i = 0; i < path->numPoints(); ++i) {
    Point &pt = path->pointN(i);
    pt.z()    = i * 0.1;
  }

  // Rectangle from 0 to 2 in X, 0 to 1 in Y
  Polygon    rect;
  LineString ring;
  ring.addPoint(Point(0, 0));
  ring.addPoint(Point(2, 0));
  ring.addPoint(Point(2, 1));
  ring.addPoint(Point(0, 1));
  ring.addPoint(Point(0, 0));
  rect.addRing(ring);
  auto profile = std::make_unique<Polygon>(rect);

  // Sweep with anchor at rectangle center (1, 0.5)
  algorithm::SweepOptions options;
  options.anchor_x = 1.0;
  options.anchor_y = 0.5;

  auto result = algorithm::sweep(*path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
}

/**
 * @brief Test custom anchor point (offset)
 */
BOOST_AUTO_TEST_CASE(testSweep_AnchorOffset)
{
  // Simple straight path along X axis
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(10, 0, 0));
  auto path_ptr = std::make_unique<LineString>(path);

  // Simple profile: Point at (0,0) (represented as small triangle for sweep)
  // Let's use a small square centered at 0
  auto profile = algorithm::create_rectangular_profile(1.0, 1.0);

  // Sweep with anchor at (0.5, 0)
  algorithm::SweepOptions options;
  options.anchor_x = 0.5;
  options.anchor_y = 0.0;

  auto result = algorithm::sweep(*path_ptr, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
}

/**
 * @brief Test different anchor points on same path
 */
BOOST_AUTO_TEST_CASE(testSweep_MultipleAnchors)
{
  auto path    = algorithm::create_circular_profile(5.0, 16);
  auto profile = algorithm::create_rectangular_profile(1.0, 1.0);

  // Sweep with anchor at bottom-left corner (0, 0)
  algorithm::SweepOptions options1;
  options1.anchor_x = 0.0;
  options1.anchor_y = 0.0;
  auto result1      = algorithm::sweep(*path, *profile, options1);

  // Sweep with anchor at top-right corner (1, 1)
  algorithm::SweepOptions options2;
  options2.anchor_x = 1.0;
  options2.anchor_y = 1.0;
  auto result2      = algorithm::sweep(*path, *profile, options2);

  BOOST_CHECK(algorithm::isValid(*result1));
  BOOST_CHECK(algorithm::isValid(*result2));
}

/**
 * @brief Test negative anchor values
 */
BOOST_AUTO_TEST_CASE(testSweep_NegativeAnchor)
{
  auto path    = algorithm::create_circular_profile(5.0, 16);
  auto profile = algorithm::create_rectangular_profile(1.0, 1.0);

  // Sweep with negative anchor (-0.5, -0.5)
  algorithm::SweepOptions options;
  options.anchor_x = -0.5;
  options.anchor_y = -0.5;

  auto result = algorithm::sweep(*path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
}

/**
 * @brief Test closed path flag interaction with anchor
 */
BOOST_AUTO_TEST_CASE(testSweep_ClosedPathAnchor)
{
  // Create an open arc
  LineString arc;
  for (int i = 0; i <= 10; ++i) {
    double angle = i * M_PI / 10.0; // 0 to Pi
    arc.addPoint(Point(std::cos(angle) * 5, std::sin(angle) * 5, 0));
  }
  auto path    = std::make_unique<LineString>(arc);
  auto profile = algorithm::create_rectangular_profile(1.0, 1.0);

  // Sweep with closed path and custom anchor
  algorithm::SweepOptions options;
  options.closed_path = true; // Should close the semi-circle
  options.anchor_x    = 0.5;

  auto result = algorithm::sweep(*path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
}

/**
 * @brief Test Flat End Caps (Default)
 */
BOOST_AUTO_TEST_CASE(testSweep_EndCaps)
{
  // Straight path
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(10, 0, 0));

  auto profile      = algorithm::create_circular_profile(1.0, 8); // Polygon
  auto profile_poly = std::make_unique<Polygon>(*profile);

  // Sweep with flat caps (default) and custom anchor
  algorithm::SweepOptions options;
  options.start_cap = algorithm::SweepOptions::EndCapStyle::FLAT;
  options.end_cap   = algorithm::SweepOptions::EndCapStyle::FLAT;
  options.anchor_x  = 0.5; // Offset

  auto result = algorithm::sweep(path, *profile_poly, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));

  algorithm::SweepOptions options_no_caps;
  options_no_caps.start_cap = algorithm::SweepOptions::EndCapStyle::NONE;
  options_no_caps.end_cap   = algorithm::SweepOptions::EndCapStyle::NONE;

  auto result_no_caps = algorithm::sweep(path, *profile_poly, options_no_caps);
  BOOST_CHECK(result_no_caps != nullptr);
  BOOST_CHECK(algorithm::isValid(*result_no_caps));
}

/**
 * @brief Test caps on vertical line with profile starting at origin
 * Regression test for bug where caps failed to be added
 */
BOOST_AUTO_TEST_CASE(testSweep_VerticalLineWithCaps)
{
  // Vertical line from z=0 to z=1
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(0, 0, 1));

  // Circular profile (closed)
  auto profile = algorithm::create_circular_profile(0.1, 16);

  algorithm::SweepOptions options;
  options.frame_method = algorithm::SweepOptions::FrameMethod::SEGMENT_ALIGNED;
  options.start_cap    = algorithm::SweepOptions::EndCapStyle::FLAT;
  options.end_cap      = algorithm::SweepOptions::EndCapStyle::FLAT;

  auto result = algorithm::sweep(path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));

  // The result should have caps: 16 side faces + 2 cap polygons = 18 faces
  BOOST_CHECK_EQUAL(result->numGeometries(), 18);

  // Result should be closed (can be converted to SOLID)
  BOOST_CHECK(isClosed(*result));
}

/**
 * @brief Test independent cap control
 */
BOOST_AUTO_TEST_CASE(testSweep_IndependentCapControl)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(0, 0, 1));

  auto profile = algorithm::create_rectangular_profile(0.2, 0.2);

  // Test with both caps
  algorithm::SweepOptions options_both;
  options_both.start_cap = algorithm::SweepOptions::EndCapStyle::FLAT;
  options_both.end_cap   = algorithm::SweepOptions::EndCapStyle::FLAT;
  auto result_both       = algorithm::sweep(path, *profile, options_both);
  BOOST_CHECK_EQUAL(result_both->numGeometries(), 6); // 4 sides + 2 caps

  // Test with only start cap
  algorithm::SweepOptions options_start;
  options_start.start_cap = algorithm::SweepOptions::EndCapStyle::FLAT;
  options_start.end_cap   = algorithm::SweepOptions::EndCapStyle::NONE;
  auto result_start       = algorithm::sweep(path, *profile, options_start);
  BOOST_CHECK_EQUAL(result_start->numGeometries(), 5); // 4 sides + start cap

  // Test with only end cap
  algorithm::SweepOptions options_end;
  options_end.start_cap = algorithm::SweepOptions::EndCapStyle::NONE;
  options_end.end_cap   = algorithm::SweepOptions::EndCapStyle::FLAT;
  auto result_end       = algorithm::sweep(path, *profile, options_end);
  BOOST_CHECK_EQUAL(result_end->numGeometries(), 5); // 4 sides + end cap

  // Test with no caps
  algorithm::SweepOptions options_none;
  options_none.start_cap = algorithm::SweepOptions::EndCapStyle::NONE;
  options_none.end_cap   = algorithm::SweepOptions::EndCapStyle::NONE;
  auto result_none       = algorithm::sweep(path, *profile, options_none);
  BOOST_CHECK_EQUAL(result_none->numGeometries(), 4); // 4 sides only
}

BOOST_AUTO_TEST_SUITE_END()
