// Copyright (c) 2025-2026, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/LineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
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
BOOST_AUTO_TEST_CASE(testDefaultAnchor)
{
  // Simple straight path
  std::unique_ptr<LineString> path(new LineString());
  path->addPoint(Point(0, 0, 0));
  path->addPoint(Point(10, 0, 0));

  // Square profile centered at origin
  std::unique_ptr<Polygon> profile(new Polygon());
  LineString               &ring = profile->exteriorRing();
  ring.addPoint(Point(-0.5, -0.5));
  ring.addPoint(Point(0.5, -0.5));
  ring.addPoint(Point(0.5, 0.5));
  ring.addPoint(Point(-0.5, 0.5));
  ring.addPoint(Point(-0.5, -0.5));

  // Sweep with default options (anchor = 0,0)
  algorithm::SweepOptions options;
  auto result = algorithm::sweep(*path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK(result->numPolygons() > 0);

  // Profile origin (0,0,0) should be on the path
  // First point of path is (0,0,0), so transformed profile should contain
  // points offset from (0,0,0)
  BOOST_CHECK(result->isValid());
}

/**
 * @brief Test custom anchor point with Polygon profile
 * Rectangle with anchor at center should center the profile on path
 */
BOOST_AUTO_TEST_CASE(testCustomAnchorPolygon)
{
  // Simple straight path
  std::unique_ptr<LineString> path(new LineString());
  path->addPoint(Point(0, 0, 0));
  path->addPoint(Point(10, 0, 0));

  // Rectangle profile: bottom-left at (0,0), top-right at (2,1)
  std::unique_ptr<Polygon> profile(new Polygon());
  LineString               &ring = profile->exteriorRing();
  ring.addPoint(Point(0, 0));
  ring.addPoint(Point(2, 0));
  ring.addPoint(Point(2, 1));
  ring.addPoint(Point(0, 1));
  ring.addPoint(Point(0, 0));

  // Sweep with anchor at rectangle center (1, 0.5)
  algorithm::SweepOptions options;
  options.anchor_x = 1.0;
  options.anchor_y = 0.5;

  auto result = algorithm::sweep(*path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK(result->isValid());

  // With anchor (1, 0.5), the point (1, 0.5) of profile should be on path
  // This centers the rectangle on the path
  BOOST_CHECK(result->numPolygons() > 0);
}

/**
 * @brief Test anchor with LineString profile (open profile)
 * Anchor should work the same way for open profiles
 */
BOOST_AUTO_TEST_CASE(testCustomAnchorLineString)
{
  // Simple straight path
  std::unique_ptr<LineString> path(new LineString());
  path->addPoint(Point(0, 0, 0));
  path->addPoint(Point(10, 0, 0));

  // LineString profile (horizontal line from -1 to 1)
  std::unique_ptr<LineString> profile(new LineString());
  profile->addPoint(Point(-1, 0));
  profile->addPoint(Point(1, 0));

  // Sweep with anchor at (0.5, 0)
  algorithm::SweepOptions options;
  options.anchor_x = 0.5;
  options.anchor_y = 0.0;

  auto result = algorithm::sweep(*path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK(result->isValid());

  // The point (0.5, 0) of profile should be positioned on path
  BOOST_CHECK(result->numPolygons() > 0);
}

/**
 * @brief Test anchor at corner of profile
 * Anchor at profile corner should position that corner on path
 */
BOOST_AUTO_TEST_CASE(testAnchorAtCorner)
{
  // Simple straight path
  std::unique_ptr<LineString> path(new LineString());
  path->addPoint(Point(0, 0, 0));
  path->addPoint(Point(5, 0, 0));

  // Square profile from (0,0) to (1,1)
  std::unique_ptr<Polygon> profile(new Polygon());
  LineString               &ring = profile->exteriorRing();
  ring.addPoint(Point(0, 0));
  ring.addPoint(Point(1, 0));
  ring.addPoint(Point(1, 1));
  ring.addPoint(Point(0, 1));
  ring.addPoint(Point(0, 0));

  // Sweep with anchor at bottom-left corner (0, 0)
  algorithm::SweepOptions options1;
  options1.anchor_x = 0.0;
  options1.anchor_y = 0.0;
  auto result1 = algorithm::sweep(*path, *profile, options1);
  BOOST_CHECK(result1 != nullptr);
  BOOST_CHECK(result1->isValid());

  // Sweep with anchor at top-right corner (1, 1)
  algorithm::SweepOptions options2;
  options2.anchor_x = 1.0;
  options2.anchor_y = 1.0;
  auto result2 = algorithm::sweep(*path, *profile, options2);
  BOOST_CHECK(result2 != nullptr);
  BOOST_CHECK(result2->isValid());

  // Both should produce valid geometries
  BOOST_CHECK(result1->numPolygons() > 0);
  BOOST_CHECK(result2->numPolygons() > 0);
}

/**
 * @brief Test negative anchor coordinates
 * Negative anchor values should create valid offset geometry
 */
BOOST_AUTO_TEST_CASE(testNegativeAnchor)
{
  // Simple straight path
  std::unique_ptr<LineString> path(new LineString());
  path->addPoint(Point(0, 0, 0));
  path->addPoint(Point(10, 0, 0));

  // Square profile from (0,0) to (1,1)
  std::unique_ptr<Polygon> profile(new Polygon());
  LineString               &ring = profile->exteriorRing();
  ring.addPoint(Point(0, 0));
  ring.addPoint(Point(1, 0));
  ring.addPoint(Point(1, 1));
  ring.addPoint(Point(0, 1));
  ring.addPoint(Point(0, 0));

  // Sweep with negative anchor (-0.5, -0.5)
  // This should offset the profile in positive direction
  algorithm::SweepOptions options;
  options.anchor_x = -0.5;
  options.anchor_y = -0.5;

  auto result = algorithm::sweep(*path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK(result->isValid());
  BOOST_CHECK(result->numPolygons() > 0);
}

/**
 * @brief Test anchor with closed path
 * Anchor should work correctly with closed loop paths
 */
BOOST_AUTO_TEST_CASE(testAnchorWithClosedPath)
{
  // Closed rectangular path
  std::unique_ptr<LineString> path(new LineString());
  path->addPoint(Point(0, 0, 0));
  path->addPoint(Point(5, 0, 0));
  path->addPoint(Point(5, 5, 0));
  path->addPoint(Point(0, 5, 0));
  path->addPoint(Point(0, 0, 0)); // Closes the loop

  // Small square profile
  std::unique_ptr<Polygon> profile(new Polygon());
  LineString               &ring = profile->exteriorRing();
  ring.addPoint(Point(-0.2, -0.2));
  ring.addPoint(Point(0.2, -0.2));
  ring.addPoint(Point(0.2, 0.2));
  ring.addPoint(Point(-0.2, 0.2));
  ring.addPoint(Point(-0.2, -0.2));

  // Sweep with closed path and custom anchor
  algorithm::SweepOptions options;
  options.closed_path = true;
  options.anchor_x    = 0.1;
  options.anchor_y    = 0.1;

  auto result = algorithm::sweep(*path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK(result->isValid());
  BOOST_CHECK(result->numPolygons() > 0);
}

/**
 * @brief Test anchor with end caps
 * Anchor should work correctly when end caps are enabled
 */
BOOST_AUTO_TEST_CASE(testAnchorWithCaps)
{
  // Simple straight path
  std::unique_ptr<LineString> path(new LineString());
  path->addPoint(Point(0, 0, 0));
  path->addPoint(Point(10, 0, 0));

  // Circular-ish profile (octagon)
  std::unique_ptr<Polygon> profile(new Polygon());
  LineString               &ring = profile->exteriorRing();
  const double              r    = 0.5;
  const int                 n    = 8;
  for (int i = 0; i <= n; ++i) {
    double angle = 2.0 * M_PI * i / n;
    ring.addPoint(Point(r * cos(angle), r * sin(angle)));
  }

  // Sweep with flat caps (default) and custom anchor
  algorithm::SweepOptions options;
  options.start_cap = algorithm::SweepOptions::EndCapStyle::FLAT;
  options.end_cap   = algorithm::SweepOptions::EndCapStyle::FLAT;
  options.anchor_x  = 0.25; // Offset anchor
  options.anchor_y  = 0.0;

  auto result = algorithm::sweep(*path, *profile, options);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK(result->isValid());

  // Should have end caps plus tube walls
  BOOST_CHECK(result->numPolygons() > n); // At least tube segments + 2 caps

  // Test without caps
  algorithm::SweepOptions options_no_caps;
  options_no_caps.start_cap = algorithm::SweepOptions::EndCapStyle::NONE;
  options_no_caps.end_cap   = algorithm::SweepOptions::EndCapStyle::NONE;
  options_no_caps.anchor_x  = 0.25;
  options_no_caps.anchor_y  = 0.0;

  auto result_no_caps = algorithm::sweep(*path, *profile, options_no_caps);

  BOOST_CHECK(result_no_caps != nullptr);
  BOOST_CHECK(result_no_caps->isValid());

  // Should have fewer polygons without caps
  BOOST_CHECK(result_no_caps->numPolygons() < result->numPolygons());
}

BOOST_AUTO_TEST_SUITE_END()
