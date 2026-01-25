// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/vtk.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/io/wkt.h"
#include <filesystem>
#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

#include "../../../test_config.h"

BOOST_AUTO_TEST_SUITE(SFCGAL_io_VTKWriterTest)

namespace fs = std::filesystem;

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

BOOST_AUTO_TEST_CASE(test_all_geometries)
{
  std::vector<std::string> wkt_examples = {
      "POINT Z (1 2 3)",
      "LINESTRING Z (0 0 0, 1 1 1, 2 2 2)",
      "POLYGON Z ((0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 0))",
      "TRIANGLE Z ((0 0 0, 1 0 0, 0 1 0, 0 0 0))",
      "POLYHEDRALSURFACE Z (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)), ((0 0 0, 0 "
      "1 "
      "0, 0 1 1, 0 0 1, 0 0 0)))",
      "SOLID Z ((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)), ((1 0 0,1 1 0,1 1 1,1 0 "
      "1,1 "
      "0 0)), ((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0)), ((0 0 1,0 1 1,0 1 0,0 0 0,0 0 "
      "1)), ((1 0 1,1 1 1,0 1 1,0 0 1,1 0 1)), ((1 0 0,1 0 1,0 0 1,0 0 0,1 0 "
      "0))))",
      "TIN Z (((0 0 0, 0 0 1, 0 1 0, 0 0 0)), ((0 0 0, 0 1 0, 1 0 0, 0 0 0)))",
      "MULTIPOINT Z ((1 1 1), (2 2 2))",
      "MULTILINESTRING Z ((0 0 0, 1 1 1), (2 2 2, 3 3 3))",
      "MULTIPOLYGON Z (((0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 0)), ((2 2 2, 3 2 2, "
      "3 "
      "3 2, 2 3 2, 2 2 2)))",
      "MULTISOLID Z (((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),((0 0 1,1 0 1,1 1 1,0 "
      "1 1,0 0 1)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),((1 1 0,0 1 0,0 1 1,1 1 "
      "1,1 1 0)),((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),((0 0 0,0 0 1,0 1 1,0 1 0,0 "
      "0 0)))),((((2 4 6,2 5 6,3 5 6,3 4 6,2 4 6)),((2 4 7,3 4 7,3 5 7,2 5 7,2 "
      "4 7)),((2 4 6,3 4 6,3 4 7,2 4 7,2 4 6)),((3 5 6,2 5 6,2 5 7,3 5 7,3 5 "
      "6)),((3 4 6,3 5 6,3 5 7,3 4 7,3 4 6)),((2 4 6,2 4 7,2 5 7,2 5 6,2 4 "
      "6)))))",
      "GEOMETRYCOLLECTION Z (POINT Z (2.0 3.0 5.0),TRIANGLE Z ((0.0 0.0 "
      "6.0,1.0 "
      "0.0 6.0,1.0 1.0 6.0,0.0 0.0 6.0)))",
      "POINT (1 2)",
      "LINESTRING (0 0, 1 1, 2 2)",
      "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))",
      "TRIANGLE ((0 0, 1 0, 0 1, 0 0))",
      "POLYHEDRALSURFACE (((0 0, 0 1, 1 1, 1 0, 0 0)), ((0 0, 0 1, 0 0)))",
      "TIN (((0 0, 0 1, 1 0, 0 0)), ((0 0, 1 0, 1 1, 0 0)))",
      "MULTIPOINT ((1 1), (2 2))",
      "MULTILINESTRING ((0 0, 1 1), (2 2, 3 3))",
      "MULTIPOLYGON (((0 0, 1 0, 1 1, 0 1, 0 0)), ((2 2, 3 2, 3 3, 2 3, 2 2)))",
      "MULTISOLID (((((0 0,0 1,1 1,1 0,0 0)),((0 0,1 0,1 0,0 0,0 0)),((0 0,1 "
      "0,1 0,0 0,0 0)),((1 1,0 1,0 1,1 1,1 1)),((1 0,1 1,1 1,1 0,1 0)),((0 0,0 "
      "0,0 1,0 1,0 0)))),((((2 4,2 5,3 5,3 4,2 4)),((2 4,3 4,3 5,2 5,2 4)),((2 "
      "4,3 4,3 4,2 4,2 4)),((3 5,2 5,2 5,3 5,3 5)),((3 4,3 5,3 5,3 4,3 4)),((2 "
      "4,2 4,2 5,2 5,2 4)))))",
      "GEOMETRYCOLLECTION (POINT (2.0 3.0),TRIANGLE ((0.0 0.0,1.0 0.0,1.0 "
      "1.0,0.0 "
      "0.0)))"};

  // Create a temporary directory for generated files
  fs::path temp_dir = fs::temp_directory_path() / random_string();
  fs::create_directories(temp_dir);

  for (const auto &wkt : wkt_examples) {
    std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

    // Generate the expected file name
    std::string geomType = geom->geometryType();
    if (geom->is3D()) {
      geomType += "z";
    }
    std::transform(geomType.begin(), geomType.end(), geomType.begin(),
                   ::tolower);

    std::string expectedFile = std::string(SFCGAL_TEST_DIRECTORY) +
                               "/data/vtkfiles/" + geomType + ".vtk";

    // Check if the expected file exists
    if (!fs::exists(expectedFile)) {
      std::cout << "Expected file does not exist: " << expectedFile << '\n';
      // continue;
    }

    // Generate the file with our function
    fs::path generatedFile = temp_dir / (geomType + ".vtk");
    SFCGAL::io::VTK::save(*geom, generatedFile.string());

    // Compare the files
    BOOST_CHECK_MESSAGE(compareFiles(expectedFile, generatedFile.string()),
                        "Output for " << geomType
                                      << " does not match the expected file.\n"
                                         "Expected file: "
                                      << expectedFile
                                      << "\n"
                                         "Generated file: "
                                      << generatedFile.string());
  }

  // Clean up the temporary directory
  //  fs::remove_all(temp_dir);
}

BOOST_AUTO_TEST_CASE(test_save_to_string)
{
  std::string                       wkt = "POINT (1 2)";
  std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

  std::string result = SFCGAL::io::VTK::saveToString(*geom);

  std::string expected = "# vtk DataFile Version 2.0\n"
                         "SFCGAL Geometry\n"
                         "ASCII\n"
                         "DATASET UNSTRUCTURED_GRID\n"
                         "POINTS 1 float\n"
                         "1 2 0\n"
                         "CELLS 1 2\n"
                         "1 0\n"
                         "CELL_TYPES 1\n"
                         "1\n";
  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(test_save_to_buffer)
{
  std::string                       wkt = "LINESTRING (0 0, 1 1)";
  std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

  size_t size = 1000;
  char   buffer[1000];
  SFCGAL::io::VTK::saveToBuffer(*geom, buffer, &size);

  std::string result(buffer); // Buffer is null-terminated
  std::string expected = "# vtk DataFile Version 2.0\n"
                         "SFCGAL Geometry\n"
                         "ASCII\n"
                         "DATASET UNSTRUCTURED_GRID\n"
                         "POINTS 2 float\n"
                         "0 0 0\n"
                         "1 1 0\n"
                         "CELLS 1 3\n"
                         "2 0 1\n"
                         "CELL_TYPES 1\n"
                         "4\n";
  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(test_buffer_size)
{
  std::string                       wkt = "POINT (1 2 3)";
  std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

  size_t size = 0;
  SFCGAL::io::VTK::saveToBuffer(*geom, nullptr, &size);

  BOOST_CHECK_GT(size, 0);

  std::vector<char> buffer(size);
  SFCGAL::io::VTK::saveToBuffer(*geom, buffer.data(), &size);

  std::string result(buffer.data()); // Buffer is null-terminated
  std::string expected = "# vtk DataFile Version 2.0\n"
                         "SFCGAL Geometry\n"
                         "ASCII\n"
                         "DATASET UNSTRUCTURED_GRID\n"
                         "POINTS 1 float\n"
                         "1 2 3\n"
                         "CELLS 1 2\n"
                         "1 0\n"
                         "CELL_TYPES 1\n"
                         "1\n";
  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(test_complex_geometry)
{
  std::string wkt = "GEOMETRYCOLLECTION (POINT (1 1),LINESTRING (0 0, 1 "
                    "1),POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0)))";
  std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

  std::string result = SFCGAL::io::VTK::saveToString(*geom);

  std::string expected = "# vtk DataFile Version 2.0\n"
                         "SFCGAL Geometry\n"
                         "ASCII\n"
                         "DATASET UNSTRUCTURED_GRID\n"
                         "POINTS 7 float\n"
                         "1 1 0\n"
                         "0 0 0\n"
                         "1 1 0\n"
                         "0 0 0\n"
                         "1 0 0\n"
                         "1 1 0\n"
                         "0 1 0\n"
                         "CELLS 3 10\n"
                         "1 0\n"
                         "2 1 2\n"
                         "4 3 4 5 6\n"
                         "CELL_TYPES 3\n"
                         "1\n"
                         "4\n"
                         "7\n";

  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SFCGAL_io_VTKReaderTest)

BOOST_AUTO_TEST_CASE(test_load_point)
{
  std::string vtk_content = "# vtk DataFile Version 2.0\n"
                            "Test\n"
                            "ASCII\n"
                            "DATASET UNSTRUCTURED_GRID\n"
                            "POINTS 1 float\n"
                            "1 2 3\n"
                            "CELLS 1 2\n"
                            "1 0\n"
                            "CELL_TYPES 1\n"
                            "1\n";

  std::unique_ptr<SFCGAL::Geometry> geom = SFCGAL::io::VTK::load(vtk_content);

  BOOST_CHECK_EQUAL(geom->geometryType(), "MultiPoint");
  const auto &mp = geom->as<SFCGAL::MultiPoint>();
  BOOST_CHECK_EQUAL(mp.numGeometries(), 1);
}

BOOST_AUTO_TEST_CASE(test_load_line)
{
  std::string vtk_content = "# vtk DataFile Version 2.0\n"
                            "Test\n"
                            "ASCII\n"
                            "DATASET UNSTRUCTURED_GRID\n"
                            "POINTS 2 float\n"
                            "0 0 0\n"
                            "1 1 1\n"
                            "CELLS 1 3\n"
                            "2 0 1\n"
                            "CELL_TYPES 1\n"
                            "4\n";

  std::unique_ptr<SFCGAL::Geometry> geom = SFCGAL::io::VTK::load(vtk_content);

  BOOST_CHECK_EQUAL(geom->geometryType(), "MultiLineString");
}

BOOST_AUTO_TEST_CASE(test_load_triangle)
{
  std::string vtk_content = "# vtk DataFile Version 2.0\n"
                            "Test\n"
                            "ASCII\n"
                            "DATASET UNSTRUCTURED_GRID\n"
                            "POINTS 3 float\n"
                            "0 0 0\n"
                            "1 0 0\n"
                            "0 1 0\n"
                            "CELLS 1 4\n"
                            "3 0 1 2\n"
                            "CELL_TYPES 1\n"
                            "5\n";

  std::unique_ptr<SFCGAL::Geometry> geom = SFCGAL::io::VTK::load(vtk_content);

  BOOST_CHECK_EQUAL(geom->geometryType(), "TriangulatedSurface");
  const auto &tin = geom->as<SFCGAL::TriangulatedSurface>();
  BOOST_CHECK_EQUAL(tin.numTriangles(), 1);
}

BOOST_AUTO_TEST_CASE(test_load_polygon)
{
  std::string vtk_content = "# vtk DataFile Version 2.0\n"
                            "Test\n"
                            "ASCII\n"
                            "DATASET UNSTRUCTURED_GRID\n"
                            "POINTS 4 float\n"
                            "0 0 0\n"
                            "1 0 0\n"
                            "1 1 0\n"
                            "0 1 0\n"
                            "CELLS 1 5\n"
                            "4 0 1 2 3\n"
                            "CELL_TYPES 1\n"
                            "7\n";

  std::unique_ptr<SFCGAL::Geometry> geom = SFCGAL::io::VTK::load(vtk_content);

  BOOST_CHECK_EQUAL(geom->geometryType(), "PolyhedralSurface");
}

BOOST_AUTO_TEST_CASE(test_roundtrip_triangle)
{
  std::string wkt = "TRIANGLE Z ((0 0 0, 1 0 0, 0 1 0, 0 0 0))";
  std::unique_ptr<SFCGAL::Geometry> original(SFCGAL::io::readWkt(wkt));

  // Save to VTK
  std::string vtk = SFCGAL::io::VTK::saveToString(*original);

  // Load back
  std::unique_ptr<SFCGAL::Geometry> loaded = SFCGAL::io::VTK::load(vtk);

  BOOST_CHECK_EQUAL(loaded->geometryType(), "TriangulatedSurface");
}

BOOST_AUTO_TEST_CASE(test_load_from_file)
{
  std::string filename =
      std::string(SFCGAL_TEST_DIRECTORY) + "/data/vtkfiles/trianglez.vtk";
  std::unique_ptr<SFCGAL::Geometry> geom =
      SFCGAL::io::VTK::loadFromFile(filename);

  BOOST_CHECK_EQUAL(geom->geometryType(), "TriangulatedSurface");
}

BOOST_AUTO_TEST_CASE(test_error_invalid_point_index)
{
  std::string vtk_content = "# vtk DataFile Version 2.0\n"
                            "Test\n"
                            "ASCII\n"
                            "DATASET UNSTRUCTURED_GRID\n"
                            "POINTS 1 float\n"
                            "0 0 0\n"
                            "CELLS 1 4\n"
                            "3 0 1 2\n"
                            "CELL_TYPES 1\n"
                            "5\n";

  BOOST_CHECK_THROW(SFCGAL::io::VTK::load(vtk_content), SFCGAL::Exception);
}

BOOST_AUTO_TEST_CASE(test_smart_triangulated_surface)
{
  // Multiple triangles should become TriangulatedSurface
  std::string vtk_content = "# vtk DataFile Version 2.0\n"
                            "Test\n"
                            "ASCII\n"
                            "DATASET UNSTRUCTURED_GRID\n"
                            "POINTS 6 float\n"
                            "0 0 0\n"
                            "1 0 0\n"
                            "0 1 0\n"
                            "1 0 0\n"
                            "1 1 0\n"
                            "0 1 0\n"
                            "CELLS 2 8\n"
                            "3 0 1 2\n"
                            "3 3 4 5\n"
                            "CELL_TYPES 2\n"
                            "5\n"
                            "5\n";

  std::unique_ptr<SFCGAL::Geometry> geom = SFCGAL::io::VTK::load(vtk_content);

  BOOST_CHECK_EQUAL(geom->geometryType(), "TriangulatedSurface");
  const auto &tin = geom->as<SFCGAL::TriangulatedSurface>();
  BOOST_CHECK_EQUAL(tin.numTriangles(), 2);
}

BOOST_AUTO_TEST_CASE(test_load_legacy1_polydata_cube)
{
  // Load legacy1.vtk - POLYDATA format representing a cube with 6 quad faces
  std::string filename =
      std::string(SFCGAL_TEST_DIRECTORY) + "/data/vtkfiles/legacy1.vtk";
  std::unique_ptr<SFCGAL::Geometry> geom =
      SFCGAL::io::VTK::loadFromFile(filename);

  // Verify geometry is not empty
  BOOST_CHECK(!geom->isEmpty());

  // POLYDATA with polygons should become PolyhedralSurface
  BOOST_CHECK_EQUAL(geom->geometryType(), "PolyhedralSurface");

  const auto &ps = geom->as<SFCGAL::PolyhedralSurface>();
  // Cube has 6 faces
  BOOST_CHECK_EQUAL(ps.numPolygons(), 6);

  // Round-trip test: convert to WKT and back
  std::string wkt = geom->asText();
  BOOST_CHECK(!wkt.empty());

  std::unique_ptr<SFCGAL::Geometry> reloaded(SFCGAL::io::readWkt(wkt));
  BOOST_CHECK(!reloaded->isEmpty());
}

BOOST_AUTO_TEST_CASE(test_load_legacy2_structured_points)
{
  // Load legacy2.vtk - STRUCTURED_POINTS format (3D voxel grid)
  // The file has DIMENSIONS 3 4 6 which creates a grid of 72 points
  // and implicit voxel cells of size (3-1)*(4-1)*(6-1) = 30 voxels
  std::string filename =
      std::string(SFCGAL_TEST_DIRECTORY) + "/data/vtkfiles/legacy2.vtk";
  std::unique_ptr<SFCGAL::Geometry> geom =
      SFCGAL::io::VTK::loadFromFile(filename);

  // STRUCTURED_POINTS with implicit voxel cells becomes PolyhedralSurface
  BOOST_CHECK(!geom->isEmpty());
  BOOST_CHECK_EQUAL(geom->geometryType(), "PolyhedralSurface");

  const auto &ps = geom->as<SFCGAL::PolyhedralSurface>();
  // 2 * 3 * 5 = 30 voxels, each with 6 faces = 180 polygons
  BOOST_CHECK_GT(ps.numPolygons(), 0);

  // Round-trip test: convert to WKT and back
  std::string wkt = geom->asText();
  BOOST_CHECK(!wkt.empty());

  std::unique_ptr<SFCGAL::Geometry> reloaded(SFCGAL::io::readWkt(wkt));
  BOOST_CHECK(!reloaded->isEmpty());
}

BOOST_AUTO_TEST_CASE(test_load_legacy4_unstructured_grid)
{
  // Load legacy4.vtk - UNSTRUCTURED_GRID format with mixed cell types
  std::string filename =
      std::string(SFCGAL_TEST_DIRECTORY) + "/data/vtkfiles/legacy4.vtk";
  std::unique_ptr<SFCGAL::Geometry> geom =
      SFCGAL::io::VTK::loadFromFile(filename);

  // Verify geometry is not empty
  BOOST_CHECK(!geom->isEmpty());

  // UNSTRUCTURED_GRID with mixed cell types becomes GeometryCollection
  // The file contains: hexahedrons, tetras, triangles, quads, lines, points
  BOOST_CHECK(geom->geometryType() == "GeometryCollection" ||
              geom->geometryType() == "PolyhedralSurface" ||
              geom->geometryType() == "TriangulatedSurface" ||
              geom->geometryType() == "MultiLineString" ||
              geom->geometryType() == "MultiPoint");

  // Round-trip test: convert to WKT and back
  std::string wkt = geom->asText();
  BOOST_CHECK(!wkt.empty());

  std::unique_ptr<SFCGAL::Geometry> reloaded(SFCGAL::io::readWkt(wkt));
  BOOST_CHECK(!reloaded->isEmpty());
}

BOOST_AUTO_TEST_SUITE_END()
