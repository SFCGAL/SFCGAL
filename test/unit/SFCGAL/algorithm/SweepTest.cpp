// Copyright (c) 2025-2026, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/LineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/algorithm/isClosed.h>
#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/algorithm/sweep.h>

#include <boost/test/unit_test.hpp>

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_SweepTest)

// ---------------------------------------------------------------------------
// Straight paths (continuous / RMF sweep)
// Profile points are quads (4 lateral + 2 flat caps = 6 patches).
// ---------------------------------------------------------------------------

/**
 * Sweep a 2×1 rectangle along a 10-unit straight X-axis path.
 * Produces a box: 4 lateral quad faces + start cap + end cap = 6 patches.
 */
BOOST_AUTO_TEST_CASE(testSweep_StraightPath_RectProfile)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(10, 0, 0));

  auto profile = algorithm::create_rectangular_profile(2.0, 1.0);

  auto result = algorithm::sweep(path, *profile);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::isClosed(*result));
  BOOST_CHECK_EQUAL(result->as<PolyhedralSurface>().numPatches(), 6);
}

/**
 * Sweep a regular 8-sided circular profile along a straight X-axis path
 * using the default Rotation Minimizing Frames method.
 * build_sweep_mesh creates quad faces: 8 lateral quads + 2 flat caps = 10 patches.
 * Without caps: 8 patches (open tube).
 */
BOOST_AUTO_TEST_CASE(testSweep_StraightPath_CircProfile_WithCaps)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(10, 0, 0));

  auto profile = algorithm::create_circular_profile(1.0, 8);

  // With flat caps
  {
    algorithm::SweepOptions opts;
    opts.start_cap = algorithm::SweepOptions::EndCapStyle::FLAT;
    opts.end_cap   = algorithm::SweepOptions::EndCapStyle::FLAT;

    auto result = algorithm::sweep(path, *profile, opts);
    BOOST_REQUIRE(result != nullptr);
    BOOST_CHECK(algorithm::isValid(*result));
    BOOST_CHECK(algorithm::isClosed(*result));
    BOOST_CHECK_EQUAL(result->as<PolyhedralSurface>().numPatches(), 10);
  }

  // Without caps: open tube, not closed
  {
    algorithm::SweepOptions opts;
    opts.start_cap = algorithm::SweepOptions::EndCapStyle::NONE;
    opts.end_cap   = algorithm::SweepOptions::EndCapStyle::NONE;

    auto result = algorithm::sweep(path, *profile, opts);
    BOOST_REQUIRE(result != nullptr);
    BOOST_CHECK(algorithm::isValid(*result));
    BOOST_CHECK(!algorithm::isClosed(*result));
    BOOST_CHECK_EQUAL(result->as<PolyhedralSurface>().numPatches(), 8);
  }
}

// ---------------------------------------------------------------------------
// Multi-segment paths (SEGMENT_ALIGNED, miter joins at corners)
// ---------------------------------------------------------------------------

/**
 * L-shaped 2-segment path (0,0,0)→(5,0,0)→(5,5,0) with a 1×1 square profile.
 *
 * SEGMENT_ALIGNED creates miter joins at the bend. Flat caps close the ends.
 * Per segment: 4 profile sides × 2 triangles = 8 lateral triangles.
 * 2 segments × 8 + 2 caps = 18 patches total.
 */
BOOST_AUTO_TEST_CASE(testSweep_LPath_MiterJoin)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(5, 0, 0));
  path.addPoint(Point(5, 5, 0));

  auto profile = algorithm::create_rectangular_profile(1.0, 1.0);

  algorithm::SweepOptions opts;
  opts.frame_method = algorithm::SweepOptions::FrameMethod::SEGMENT_ALIGNED;
  opts.start_cap    = algorithm::SweepOptions::EndCapStyle::FLAT;
  opts.end_cap      = algorithm::SweepOptions::EndCapStyle::FLAT;

  auto result = algorithm::sweep(path, *profile, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_EQUAL(result->as<PolyhedralSurface>().numPatches(), 18);
}

/**
 * Square closed path (0,0,0)→(5,0,0)→(5,5,0)→(0,5,0)→(0,0,0).
 * SEGMENT_ALIGNED + closed_path: miter joins at all 4 corners, no caps.
 * 4 segments × 8 lateral triangles = 32 patches.
 */
BOOST_AUTO_TEST_CASE(testSweep_ClosedPath_Square)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(5, 0, 0));
  path.addPoint(Point(5, 5, 0));
  path.addPoint(Point(0, 5, 0));
  path.addPoint(Point(0, 0, 0)); // geometrically closed

  auto profile = algorithm::create_rectangular_profile(0.5, 0.5);

  algorithm::SweepOptions opts;
  opts.frame_method = algorithm::SweepOptions::FrameMethod::SEGMENT_ALIGNED;

  auto result = algorithm::sweep(path, *profile, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::isClosed(*result));
  BOOST_CHECK_EQUAL(result->as<PolyhedralSurface>().numPatches(), 32);
}

// ---------------------------------------------------------------------------
// Vertical path (Z axis, challenges frame computation)
// ---------------------------------------------------------------------------

/**
 * Sweep a 16-sided circular profile along a vertical Z-axis path.
 * Verifies SEGMENT_ALIGNED frame computation for T ∥ Z_up (fallback to Y).
 * 16 sides × 2 triangles + 2 caps = 34 patches.
 */
BOOST_AUTO_TEST_CASE(testSweep_VerticalPath_CircProfile_WithCaps)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(0, 0, 1));

  auto profile = algorithm::create_circular_profile(0.1, 16);

  algorithm::SweepOptions opts;
  opts.frame_method = algorithm::SweepOptions::FrameMethod::SEGMENT_ALIGNED;
  opts.start_cap    = algorithm::SweepOptions::EndCapStyle::FLAT;
  opts.end_cap      = algorithm::SweepOptions::EndCapStyle::FLAT;

  auto result = algorithm::sweep(path, *profile, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::isClosed(*result));
  BOOST_CHECK_EQUAL(result->as<PolyhedralSurface>().numPatches(), 34);
}

// ---------------------------------------------------------------------------
// End cap control
// ---------------------------------------------------------------------------

/**
 * Independent control of start and end caps on a straight vertical path.
 * Rectangle profile = 4 lateral faces (quad); caps add 1 each.
 *
 * none/start/end/both → 4/5/5/6 patches.
 */
BOOST_AUTO_TEST_CASE(testSweep_CapControl_Independent)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(0, 0, 1));

  auto profile = algorithm::create_rectangular_profile(0.2, 0.2);

  auto make = [&](algorithm::SweepOptions::EndCapStyle s,
                  algorithm::SweepOptions::EndCapStyle e) {
    algorithm::SweepOptions opts;
    opts.start_cap = s;
    opts.end_cap   = e;
    return algorithm::sweep(path, *profile, opts);
  };
  using Cap = algorithm::SweepOptions::EndCapStyle;

  auto r_none  = make(Cap::NONE, Cap::NONE);
  auto r_start = make(Cap::FLAT, Cap::NONE);
  auto r_end   = make(Cap::NONE, Cap::FLAT);
  auto r_both  = make(Cap::FLAT, Cap::FLAT);

  BOOST_CHECK_EQUAL(r_none->as<PolyhedralSurface>().numPatches(), 4);
  BOOST_CHECK(!algorithm::isClosed(*r_none));

  BOOST_CHECK_EQUAL(r_start->as<PolyhedralSurface>().numPatches(), 5);
  BOOST_CHECK(!algorithm::isClosed(*r_start));

  BOOST_CHECK_EQUAL(r_end->as<PolyhedralSurface>().numPatches(), 5);
  BOOST_CHECK(!algorithm::isClosed(*r_end));

  BOOST_CHECK_EQUAL(r_both->as<PolyhedralSurface>().numPatches(), 6);
  BOOST_CHECK(algorithm::isClosed(*r_both));
}

// ---------------------------------------------------------------------------
// Error cases
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testSweep_Error_PathTooFewPoints)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));

  auto profile = algorithm::create_rectangular_profile(1.0, 1.0);
  BOOST_CHECK_THROW(algorithm::sweep(path, *profile), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testSweep_Error_ProfileIsPoint)
{
  LineString path;
  path.addPoint(Point(0, 0, 0));
  path.addPoint(Point(0, 0, 1));

  Point profile(1, 1, 0);
  BOOST_CHECK_THROW(algorithm::sweep(path, profile), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testSweep_Error_CircularProfile_NegativeRadius)
{
  BOOST_CHECK_THROW(algorithm::create_circular_profile(-1.0, 8),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testSweep_Error_CircularProfile_TooFewSegments)
{
  BOOST_CHECK_THROW(algorithm::create_circular_profile(1.0, 2),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testSweep_Error_RectProfile_NegativeDimension)
{
  BOOST_CHECK_THROW(algorithm::create_rectangular_profile(-1.0, 1.0),
                    std::invalid_argument);
  BOOST_CHECK_THROW(algorithm::create_rectangular_profile(1.0, -1.0),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
