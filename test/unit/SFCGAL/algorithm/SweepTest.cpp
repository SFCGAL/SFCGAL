// Copyright (c) 2025-2026, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/LineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/algorithm/isClosed.h>
#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/algorithm/sweep.h>
#include <SFCGAL/io/wkt.h>

#include <boost/test/unit_test.hpp>

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_SweepTest)

/**
 * @brief Test default anchor point behavior
 */
BOOST_AUTO_TEST_CASE(testSweep_AnchorDefault)
{
  // Straight path along X
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(10, 0, 0));

  auto profile = algorithm::create_rectangular_profile(2.0, 1.0);

  algorithm::SweepOptions options;
  auto                    result = algorithm::sweep(path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
  auto &ps = result->as<PolyhedralSurface>();
  BOOST_CHECK_EQUAL(ps.numPatches(), 6);
  BOOST_CHECK(algorithm::isClosed(*result));
}

/**
 * @brief Test custom anchor point (centered on profile)
 */
BOOST_AUTO_TEST_CASE(testSweep_AnchorCenter)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(10, 0, 0));

  // Rectangle from 0 to 2 in X, 0 to 1 in Y — use Polygon(LineString)
  LineString ring;
  ring.addPoint(Point(0, 0));
  ring.addPoint(Point(2, 0));
  ring.addPoint(Point(2, 1));
  ring.addPoint(Point(0, 1));
  ring.addPoint(Point(0, 0));
  Polygon profile(ring);

  // Sweep with anchor at rectangle center (1, 0.5)
  algorithm::SweepOptions options;
  options.anchor_x = 1.0;
  options.anchor_y = 0.5;

  auto result = algorithm::sweep(path, profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
  auto &ps = result->as<PolyhedralSurface>();
  BOOST_CHECK_EQUAL(ps.numPatches(), 6);
  BOOST_CHECK(algorithm::isClosed(*result));
}

/**
 * @brief Test custom anchor point (offset)
 */
BOOST_AUTO_TEST_CASE(testSweep_AnchorOffset)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(10, 0, 0));

  auto profile = algorithm::create_rectangular_profile(1.0, 1.0);

  algorithm::SweepOptions options;
  options.anchor_x = 0.5;
  options.anchor_y = 0.0;

  auto result = algorithm::sweep(path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
  auto &ps = result->as<PolyhedralSurface>();
  BOOST_CHECK_EQUAL(ps.numPatches(), 6);
  BOOST_CHECK(algorithm::isClosed(*result));
}

/**
 * @brief Test different anchor points on same path
 */
BOOST_AUTO_TEST_CASE(testSweep_MultipleAnchors)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(5, 0, 0));
  path.addPoint(Point(5, 5, 0));

  auto profile = algorithm::create_rectangular_profile(1.0, 1.0);

  algorithm::SweepOptions options1;
  options1.anchor_x     = 0.0;
  options1.anchor_y     = 0.0;
  options1.frame_method = algorithm::SweepOptions::FrameMethod::SEGMENT_ALIGNED;
  auto result1          = algorithm::sweep(path, *profile, options1);

  algorithm::SweepOptions options2;
  options2.anchor_x     = 1.0;
  options2.anchor_y     = 1.0;
  options2.frame_method = algorithm::SweepOptions::FrameMethod::SEGMENT_ALIGNED;
  auto result2          = algorithm::sweep(path, *profile, options2);

  BOOST_CHECK(result1 != nullptr);
  BOOST_CHECK(result2 != nullptr);
  BOOST_CHECK_EQUAL(result1->as<PolyhedralSurface>().numPatches(), 10);
  BOOST_CHECK(algorithm::isValid(*result1));
}

/**
 * @brief Test closed path flag
 */
BOOST_AUTO_TEST_CASE(testSweep_ClosedPath)
{
  // Square path (closed)
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(5, 0, 0));
  path.addPoint(Point(5, 5, 0));
  path.addPoint(Point(0, 5, 0));
  path.addPoint(Point(0, 0, 0));

  auto profile = algorithm::create_rectangular_profile(0.5, 0.5);

  algorithm::SweepOptions options;
  options.frame_method = algorithm::SweepOptions::FrameMethod::SEGMENT_ALIGNED;

  auto result = algorithm::sweep(path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_EQUAL(result->as<PolyhedralSurface>().numPatches(), 16);
  BOOST_CHECK(algorithm::isClosed(*result));
}

/**
 * @brief Test Flat End Caps (Default)
 */
BOOST_AUTO_TEST_CASE(testSweep_EndCaps)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(10, 0, 0));

  auto profile = algorithm::create_circular_profile(1.0, 8);

  algorithm::SweepOptions options;
  options.start_cap = algorithm::SweepOptions::EndCapStyle::FLAT;
  options.end_cap   = algorithm::SweepOptions::EndCapStyle::FLAT;

  auto result = algorithm::sweep(path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_EQUAL(result->as<PolyhedralSurface>().numPatches(), 10);
  BOOST_CHECK(algorithm::isClosed(*result));

  algorithm::SweepOptions options_no_caps;
  options_no_caps.start_cap = algorithm::SweepOptions::EndCapStyle::NONE;
  options_no_caps.end_cap   = algorithm::SweepOptions::EndCapStyle::NONE;

  auto result_no_caps = algorithm::sweep(path, *profile, options_no_caps);
  BOOST_CHECK(result_no_caps != nullptr);
  BOOST_CHECK(algorithm::isValid(*result_no_caps));
  BOOST_CHECK_EQUAL(result_no_caps->as<PolyhedralSurface>().numPatches(), 8);
  BOOST_CHECK(!algorithm::isClosed(*result_no_caps));
}

/**
 * @brief Test caps on vertical line
 */
BOOST_AUTO_TEST_CASE(testSweep_VerticalLineWithCaps)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(0, 0, 1));

  auto profile = algorithm::create_circular_profile(0.1, 16);

  algorithm::SweepOptions options;
  options.frame_method = algorithm::SweepOptions::FrameMethod::SEGMENT_ALIGNED;
  options.start_cap    = algorithm::SweepOptions::EndCapStyle::FLAT;
  options.end_cap      = algorithm::SweepOptions::EndCapStyle::FLAT;

  auto result = algorithm::sweep(path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));

  // 16 side quads + 2 caps = 18 polygons
  auto &ps = result->as<PolyhedralSurface>();
  BOOST_CHECK_EQUAL(ps.numPatches(), 18);

  BOOST_CHECK(algorithm::isClosed(*result));
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

  // Both caps: 4 side quads + 2 caps = 6
  algorithm::SweepOptions options_both;
  options_both.start_cap = algorithm::SweepOptions::EndCapStyle::FLAT;
  options_both.end_cap   = algorithm::SweepOptions::EndCapStyle::FLAT;
  auto result_both       = algorithm::sweep(path, *profile, options_both);
  BOOST_CHECK_EQUAL(result_both->as<PolyhedralSurface>().numPatches(), 6);
  BOOST_CHECK(algorithm::isClosed(*result_both));

  // Start cap only: 4 + 1 = 5
  algorithm::SweepOptions options_start;
  options_start.start_cap = algorithm::SweepOptions::EndCapStyle::FLAT;
  options_start.end_cap   = algorithm::SweepOptions::EndCapStyle::NONE;
  auto result_start       = algorithm::sweep(path, *profile, options_start);
  BOOST_CHECK_EQUAL(result_start->as<PolyhedralSurface>().numPatches(), 5);
  BOOST_CHECK(!algorithm::isClosed(*result_start));

  // End cap only: 4 + 1 = 5
  algorithm::SweepOptions options_end;
  options_end.start_cap = algorithm::SweepOptions::EndCapStyle::NONE;
  options_end.end_cap   = algorithm::SweepOptions::EndCapStyle::FLAT;
  auto result_end       = algorithm::sweep(path, *profile, options_end);
  BOOST_CHECK_EQUAL(result_end->as<PolyhedralSurface>().numPatches(), 5);
  BOOST_CHECK(!algorithm::isClosed(*result_end));

  // No caps: 4 sides only
  algorithm::SweepOptions options_none;
  options_none.start_cap = algorithm::SweepOptions::EndCapStyle::NONE;
  options_none.end_cap   = algorithm::SweepOptions::EndCapStyle::NONE;
  auto result_none       = algorithm::sweep(path, *profile, options_none);
  BOOST_CHECK_EQUAL(result_none->as<PolyhedralSurface>().numPatches(), 4);
  BOOST_CHECK(!algorithm::isClosed(*result_none));
}

// ---------------------------------------------------------------------------
// Error cases
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testSweep_InvalidPath_TooFewPoints)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));

  auto profile = algorithm::create_rectangular_profile(1.0, 1.0);
  BOOST_CHECK_THROW(algorithm::sweep(path, *profile), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testSweep_InvalidProfile_Point)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(0, 0, 1));

  Point profile(1, 1, 0);
  BOOST_CHECK_THROW(algorithm::sweep(path, profile), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testSweep_InvalidCircularProfile_NegativeRadius)
{
  BOOST_CHECK_THROW(algorithm::create_circular_profile(-1.0),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testSweep_InvalidCircularProfile_TooFewSegments)
{
  BOOST_CHECK_THROW(algorithm::create_circular_profile(1.0, 2),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testSweep_InvalidRectProfile_NegativeDim)
{
  BOOST_CHECK_THROW(algorithm::create_rectangular_profile(-1.0, 1.0),
                    std::invalid_argument);
  BOOST_CHECK_THROW(algorithm::create_rectangular_profile(1.0, -1.0),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testSweep_InvalidChamferProfile)
{
  BOOST_CHECK_THROW(algorithm::create_chamfer_profile(-1.0),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testSweep_InvalidFilletProfile)
{
  BOOST_CHECK_THROW(algorithm::create_fillet_profile(-1.0),
                    std::invalid_argument);
  BOOST_CHECK_THROW(algorithm::create_fillet_profile(1.0, 0),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
