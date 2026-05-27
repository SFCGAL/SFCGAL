// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#define BOOST_TEST_MODULE UnitTestSFCGAL

#define BOOST_TEST_ALTERNATIVE_INIT_API

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

#include <CGAL/assertions.h>

#include "SFCGAL/detail/tools/Log.h"

// CGAL's default msg_assert_fail prints to stderr before calling the behaviour
// handler. In tests, this creates noise even when we catch the resulting
// exceptions. Install a silent handler so assertions are still thrown (via
// set_error_behaviour in sfcgal_init) but never go to stderr.
static void
silentCGALHandler(const char * /*expr*/, const char * /*file*/,
                  const char * /*func*/, int /*line*/, const char * /*msg*/)
{
}

/**
 * @brief RAII initializer that silences CGAL assertion stderr output.
 *
 * CGAL's default handler prints to stderr before calling the behaviour
 * handler (which is set to THROW_EXCEPTION by sfcgal_init). In tests,
 * this creates noise even when the resulting exception is caught.
 * This struct installs a silent handler during static initialization.
 */
struct CGALSilentInitializer {
  CGALSilentInitializer()
  {
    CGAL::set_error_handler(silentCGALHandler);
    CGAL::set_warning_handler(silentCGALHandler);
  }
};

static CGALSilentInitializer cgalSilentInitializer; // NOLINT(cert-err58-cpp)

auto
init_unit_test_suite(int /*unused*/, char **const /*unused*/) -> test_suite *
{
  //	std::cerr << "init test suite" << std::endl;
  SFCGAL::Logger::get()->setLogLevel(SFCGAL::Logger::Info);
  // CGAL must throw exceptions (not abort) so try/catch in algorithm code
  // can handle precondition violations gracefully.
  CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
  return nullptr;
}
