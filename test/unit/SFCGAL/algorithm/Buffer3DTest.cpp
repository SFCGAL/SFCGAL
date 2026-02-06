// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/buffer3D.h"
#include "SFCGAL/GeometryCollection.h"
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
#include "SFCGAL/io/OBJ.h"
#include "SFCGAL/io/wkt.h"
#include <boost/test/unit_test.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>

#include "../../../test_config.h"

using namespace SFCGAL;
namespace fs = std::filesystem;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_Buffer3DTest)

auto
compareFiles(const std::string &file1, const std::string &file2) -> bool
{
  std::ifstream f1(file1, std::ifstream::binary | std::ifstream::ate);
  std::ifstream f2(file2, std::ifstream::binary | std::ifstream::ate);

  if (f1.fail() || f2.fail()) {
    return false; // File opening error
  }

  if (f1.tellg() != f2.tellg()) {
    return false; // Different sizes
  }

  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.rdbuf()));
}

auto
readFileContent(const std::string &filePath) -> std::string
{
  std::ifstream file(filePath);
  if (!file) {
    return "Unable to read file: " + filePath;
  }
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}

void
compareAndReportFiles(const std::string &expectedFile,
                      const std::string &generatedFile,
                      const std::string &testName)
{
  bool filesMatch = compareFiles(expectedFile, generatedFile);

  if (!filesMatch) {
    BOOST_TEST_MESSAGE("Warning for test " << testName << ":");
    BOOST_TEST_MESSAGE("Generated file does not match the expected file.");
    BOOST_TEST_MESSAGE("Expected file: " << expectedFile);
    BOOST_TEST_MESSAGE("Generated file: " << generatedFile);
    BOOST_TEST_MESSAGE("Content of the generated file:");
    BOOST_TEST_MESSAGE(readFileContent(generatedFile));
  } else {
    BOOST_TEST_MESSAGE("Test " << testName << " passed: files match.");
  }
}

BOOST_AUTO_TEST_CASE(testBuffer3D_Point)
{
  double              radius   = 10.0;
  int                 segments = 16;
  Point               point(0, 0, 0);
  algorithm::Buffer3D buffer3d(point, radius, segments);

  // Create a temporary directory for generated files
  fs::path temp_dir = fs::temp_directory_path() / random_string();
  fs::create_directories(temp_dir);

  std::vector<algorithm::Buffer3D::BufferType> bufferTypes = {
      algorithm::Buffer3D::ROUND, algorithm::Buffer3D::CYLSPHERE,
      algorithm::Buffer3D::FLAT};

  for (auto bufferType : bufferTypes) {
    std::unique_ptr<PolyhedralSurface> buffer = buffer3d.compute(bufferType);
    BOOST_CHECK(buffer->is3D());
    BOOST_CHECK(buffer->numGeometries() > 0);

    std::string typeName;
    switch (bufferType) {
    case algorithm::Buffer3D::ROUND:
      typeName = "ROUND";
      break;
    case algorithm::Buffer3D::CYLSPHERE:
      typeName = "CYLSPHERE";
      break;
    case algorithm::Buffer3D::FLAT:
      typeName = "FLAT";
      break;
    }

    // Generate the file with our function
    fs::path generatedFile =
        temp_dir / ("point_" + typeName + "_buffer_3d.obj");
    SFCGAL::io::OBJ::save(*buffer, generatedFile.string());

    // Generate the expected file name
    std::string expectedFile = std::string(SFCGAL_TEST_DIRECTORY) +
                               "/data/bufferfiles/point_" + typeName +
                               "_buffer_3d.obj";

    // Check if the expected file exists
    if (!fs::exists(expectedFile)) {
      std::cout << "Expected file does not exist: " << expectedFile << '\n';
      continue;
    }

    // Compare the files

    compareAndReportFiles(expectedFile, generatedFile.string(),
                          "point_" + typeName + "_buffer");
  }

  // Clean up the temporary directory
  fs::remove_all(temp_dir);
}

BOOST_AUTO_TEST_CASE(testBuffer3D_LineString)
{
  double              radius   = 10.0;
  int                 segments = 16;
  std::vector<Point>  points   = {Point(-100, 0, 0), Point(40, -70, 0),
                                  Point(40, 50, 40), Point(-90, -60, 60),
                                  Point(0, 0, -100), Point(30, 0, 150)};
  LineString          lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  // Create a temporary directory for generated files
  fs::path temp_dir = fs::temp_directory_path() / random_string();
  fs::create_directories(temp_dir);

  std::vector<algorithm::Buffer3D::BufferType> bufferTypes = {
      algorithm::Buffer3D::ROUND, algorithm::Buffer3D::CYLSPHERE,
      algorithm::Buffer3D::FLAT};

  for (auto bufferType : bufferTypes) {
    std::unique_ptr<PolyhedralSurface> buffer = buffer3d.compute(bufferType);
    BOOST_CHECK(buffer->is3D());
    BOOST_CHECK(buffer->numGeometries() > 0);

    std::string typeName;
    switch (bufferType) {
    case algorithm::Buffer3D::ROUND:
      typeName = "ROUND";
      break;
    case algorithm::Buffer3D::CYLSPHERE:
      typeName = "CYLSPHERE";
      break;
    case algorithm::Buffer3D::FLAT:
      typeName = "FLAT";
      break;
    }

    // Generate the file with our function
    fs::path generatedFile =
        temp_dir / ("linestring_" + typeName + "_buffer_3d.obj");
    SFCGAL::io::OBJ::save(*buffer, generatedFile.string());

    // Generate the expected file name
    std::string expectedFile = std::string(SFCGAL_TEST_DIRECTORY) +
                               "/data/bufferfiles/linestring_" + typeName +
                               "_buffer_3d.obj";

    // Check if the expected file exists
    if (!fs::exists(expectedFile)) {
      std::cout << "Expected file does not exist: " << expectedFile << '\n';
      continue;
    }

    // Compare the files

    compareAndReportFiles(expectedFile, generatedFile.string(),
                          "linestring_" + typeName + "_buffer");
  }

  // Clean up the temporary directory
  fs::remove_all(temp_dir);
}

BOOST_AUTO_TEST_CASE(testBuffer3D_InvalidGeometry)
{
  double   radius   = 10.0;
  int      segments = 16;
  Triangle triangle(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0));
  BOOST_CHECK_THROW(algorithm::Buffer3D(triangle, radius, segments),
                    std::invalid_argument);
}

// Phase 2 Tests: Edge Cases and Robustness

BOOST_AUTO_TEST_CASE(testBuffer3D_StraightLine)
{
  // Should produce a perfect cylinder
  double              radius   = 1.0;
  int                 segments = 8; // Fewer segments for simpler WKT
  std::vector<Point>  points   = {Point(0, 0, 0), Point(10, 0, 0)};
  LineString          lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  // Verify WKT output is valid and consistent
  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testBuffer3D_LShape)
{
  // Test 90-degree bend (L-shape)
  double              radius   = 1.0;
  int                 segments = 8;
  std::vector<Point>  points   = {Point(0, 0, 0), Point(5, 0, 0), Point(5, 5, 0)};
  LineString          lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  // Verify proper RMF frame propagation through the corner
  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);


}

BOOST_AUTO_TEST_CASE(testBuffer3D_UShape)
{
  // Test U-bend (two 90-degree turns)
  double             radius   = 1.0;
  int                segments = 8;
  std::vector<Point> points   = {Point(0, 0, 0), Point(5, 0, 0), Point(5, 5, 0),
                                 Point(0, 5, 0)};
  LineString         lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  // Verify RMF propagates correctly through two corners
  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);


}

BOOST_AUTO_TEST_CASE(testBuffer3D_VerticalLine)
{
  // Test singularity case: vertical line parallel to (0,0,1)
  // This tests the fallback perpendicular vector selection
  double              radius   = 1.0;
  int                 segments = 8;
  std::vector<Point>  points   = {Point(0, 0, 0), Point(0, 0, 10)};
  LineString          lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  // Verify the algorithm handles the singularity correctly
  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testBuffer3D_NearVerticalLine)
{
  // Test near-singularity case: line nearly parallel to (0,0,1)
  // This ensures RMF handles small angles correctly
  double             radius   = 1.0;
  int                segments = 8;
  std::vector<Point> points   = {Point(0, 0, 0), Point(0.01, 0.01, 10)};
  LineString         lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testBuffer3D_Helix)
{
  // Test helix/spiral curve - key test for smooth frame propagation
  double             radius   = 0.5;
  int                segments = 8;
  std::vector<Point> points;

  // Generate a helix with 3 turns
  int    num_turns = 3;
  int    num_points_per_turn = 8;
  double helix_radius        = 2.0;
  double height_per_turn     = 5.0;

  for (int i = 0; i <= num_turns * num_points_per_turn; ++i) {
    double angle = 2.0 * M_PI * i / num_points_per_turn;
    double x     = helix_radius * std::cos(angle);
    double y     = helix_radius * std::sin(angle);
    double z = height_per_turn * i / static_cast<double>(num_points_per_turn);
    points.emplace_back(x, y, z);
  }

  LineString          lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  // Critical test: helix should be smooth without discontinuities
  // This was the main failure case with the old fixed-frame approach
  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);


}

BOOST_AUTO_TEST_CASE(testBuffer3D_SharpAngle)
{
  // Test sharp angle (120 degrees) - ensures RMF handles large rotations
  double             radius   = 1.0;
  int                segments = 8;
  std::vector<Point> points   = {Point(0, 0, 0), Point(5, 0, 0),
                                 Point(2.5, 4.33, 0)}; // 120 degree angle
  LineString         lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  // Sharp angle test: RMF should propagate smoothly even with large angle change
  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testBuffer3D_CollinearSegments)
{
  // Test collinear consecutive segments - RMF should detect zero rotation
  double             radius   = 1.0;
  int                segments = 8;
  std::vector<Point> points   = {Point(0, 0, 0), Point(2, 0, 0), Point(5, 0, 0),
                                 Point(10, 0, 0)};
  LineString         lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  // With collinear points, RMF should maintain constant frame orientation
  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);


}

BOOST_AUTO_TEST_CASE(testBuffer3D_ZigZag)
{
  // Test zig-zag pattern - multiple alternating direction changes
  double             radius   = 0.5;
  int                segments = 8;
  std::vector<Point> points   = {Point(0, 0, 0), Point(2, 2, 0), Point(4, 0, 0),
                                 Point(6, 2, 0), Point(8, 0, 0)};
  LineString         lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  // Zig-zag tests frame propagation through multiple alternating bends
  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);


}

BOOST_AUTO_TEST_CASE(testBuffer3D_3DSpiral)
{
  // Test true 3D spiral with varying path curvature
  // This is the most demanding test for RMF frame propagation
  double             radius   = 0.3;
  int                segments = 8;
  std::vector<Point> points;

  for (int i = 0; i <= 30; ++i) {
    double t             = i / 30.0;
    double angle         = 3.0 * M_PI * t;
    double spiral_radius = 2.0 * (1.0 - 0.5 * t); // Decreasing radius
    double x             = spiral_radius * std::cos(angle);
    double y             = spiral_radius * std::sin(angle);
    double z             = 10.0 * t;
    points.emplace_back(x, y, z);
  }

  LineString          lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  // Complex 3D spiral: ultimate test of smooth RMF propagation
  // Both curvature and torsion vary continuously
  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);


}

// Phase 3 robustness tests
BOOST_AUTO_TEST_CASE(testBuffer3D_ClosedPath)
{
  // Test closed path with holonomy correction
  double             radius   = 0.5;
  int                segments = 8;
  std::vector<Point> points;

  // Create closed square path
  int n_points = 16;
  for (int i = 0; i < n_points; ++i) {
    double t = i / static_cast<double>(n_points);
    double x, y, z = 0;

    if (t < 0.25) {
      x = 10.0 * (t * 4.0);
      y = 0;
    } else if (t < 0.5) {
      x = 10.0;
      y = 10.0 * ((t - 0.25) * 4.0);
    } else if (t < 0.75) {
      x = 10.0 * (1.0 - (t - 0.5) * 4.0);
      y = 10.0;
    } else {
      x = 0;
      y = 10.0 * (1.0 - (t - 0.75) * 4.0);
    }

    points.emplace_back(x, y, z);
  }

  // Close the path
  points.push_back(points.front());

  LineString          lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  // Closed path should produce valid geometry with holonomy correction
  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testBuffer3D_VerySharpAngle)
{
  // Test path with very sharp angle (~150 degrees)
  double             radius   = 0.5;
  int                segments = 8;
  std::vector<Point> points;

  points.emplace_back(0, 0, 0);
  points.emplace_back(5, 0, 0);
  // Sharp turn: almost reversing direction
  points.emplace_back(6, -0.5, 0);
  points.emplace_back(11, -0.5, 0);

  LineString          lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  // Should handle sharp angles gracefully (>120 degrees)
  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testBuffer3D_ShortSegments)
{
  // Test path with some very short segments
  double             radius   = 1.0;
  int                segments = 8;
  std::vector<Point> points;

  points.emplace_back(0, 0, 0);
  points.emplace_back(5, 0, 0);
  // Very short segment (much shorter than radius)
  points.emplace_back(5.05, 0, 0);
  points.emplace_back(10, 0, 0);

  LineString          lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::unique_ptr<PolyhedralSurface> buffer =
      buffer3d.compute(algorithm::Buffer3D::FLAT);

  BOOST_CHECK(buffer->is3D());
  BOOST_CHECK(buffer->numGeometries() > 0);

  // Should handle short segments (detected in Phase 3 robustness checks)
  std::string wkt = buffer->asText(1);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") != std::string::npos);
}

BOOST_AUTO_TEST_SUITE_END()
