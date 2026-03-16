// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/tools/old/interface.hpp>
#include <vector>

#include "SFCGAL/detail/io/WktWriter_p.h"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;
using namespace SFCGAL::detail::io;

BOOST_AUTO_TEST_SUITE(SFCGAL_io_WktWriterTest)

//-- WKT POINT

BOOST_AUTO_TEST_CASE(testFixZeroNeg)
{
  std::vector<double> testValues      = {0.9999, 0.01, -0.0, 0.0,  0.9,
                                         4.9,    5.5,  5.05, 5.49, 5.51};
  std::vector<double> expectedResults = {
      1.0, 0.9999, 0.9999, -1.0, -0.9999, -0.9999, // 0.9999
      0.0, 0.0,    0.01,   0.0,  0.0,     -0.01,   // 0.01
      0.0, 0.0,    0.0,    0.0,  0.0,     0.0,     // -0.0
      0.0, 0.0,    0.0,    0.0,  0.0,     -0.0,    // 0.0
      1.0, 0.9,    0.9,    -1.0, -0.9,    -0.9,    // 0.9
      5.0, 4.9,    4.9,    -5.0, -4.9,    -4.9,    // 4.9
      6.0, 5.5,    5.5,    -6.0, -5.5,    -5.5,    // 5.5
      5.0, 5.05,   5.05,   -5,   -5.05,   -5.05,   // 5.05
      5.0, 5.49,   5.49,   -5.0, -5.49,   -5.49,   // 5.49
      6.0, 5.51,   5.51,   -6,   -5.51,   -5.51};  // 5.51

  const unsigned int nrTestValues     = testValues.size();
  const unsigned int nrExpectedValues = expectedResults.size();
  BOOST_CHECK_EQUAL(nrExpectedValues, nrTestValues * 2 * 3);

  int resultIncr = 0;
  for (const double val : testValues) {
    for (const auto valSign : {1, -1}) {
      for (int precision = 0; precision < 3; ++precision) {
        const double result = fixZeroNeg(valSign * val, precision);
        BOOST_CHECK_EQUAL(result, expectedResults[resultIncr]);
        resultIncr++;
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
