// Copyright (c) 2025-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <memory>
#include <sstream>
#include <string>

#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/io/OBJ.h"
#include "SFCGAL/io/STL.h"
#include "SFCGAL/io/vtk.h"
#include "SFCGAL/io/wkt.h"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_RecursionLimitTest)

/**
 * @brief Create a deeply nested GeometryCollection
 */
auto
createDeeplyNestedCollection(int depth) -> std::unique_ptr<Geometry>
{
  auto  root    = std::make_unique<GeometryCollection>();
  auto *current = root.get();

  for (int i = 0; i < depth - 1; ++i) {
    auto  next    = std::make_unique<GeometryCollection>();
    auto *nextPtr = next.get();
    current->addGeometry(std::move(next));
    current = nextPtr;
  }

  current->addGeometry(std::make_unique<Point>(0, 0));
  return root;
}

BOOST_AUTO_TEST_CASE(test_vtk_recursion_limit)
{
  // Default limit is 32. Create 34 levels of nesting.
  auto               deeplyNested = createDeeplyNestedCollection(34);
  std::ostringstream oss;
  BOOST_CHECK_THROW(io::VTK::save(*deeplyNested, oss), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(test_stl_recursion_limit)
{
  // Default limit is 32. Create 34 levels of nesting.
  auto               deeplyNested = createDeeplyNestedCollection(34);
  std::ostringstream oss;
  BOOST_CHECK_THROW(io::STL::save(*deeplyNested, oss), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(test_obj_recursion_limit)
{
  // Default limit is 32. Create 34 levels of nesting.
  auto               deeplyNested = createDeeplyNestedCollection(34);
  std::ostringstream oss;
  BOOST_CHECK_THROW(io::OBJ::save(*deeplyNested, oss), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(test_wkt_recursion_limit)
{
  // Construct a deeply nested WKT string
  std::string wkt = "POINT(0 0)";
  for (int i = 0; i < 34; ++i) {
    wkt.insert(0, "GEOMETRYCOLLECTION(");
    wkt += ")";
  }

  BOOST_CHECK_THROW(io::readWkt(wkt), Exception);
}

BOOST_AUTO_TEST_SUITE_END()
