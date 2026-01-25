// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/STL.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/io/io_utils.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <vector>

namespace SFCGAL::io::STL {

// Utility function to normalize "-0" strings to "0" for consistent output
auto
normalizeZeroString(const std::string &str) -> std::string
{
  return (str == "-0") ? "0" : str;
}

// Convert value using kernel and normalize -0 to 0
template <typename T>
auto
normalizeKernelValue(const T &value) -> std::string
{
  std::ostringstream oss;
  oss << value;
  return normalizeZeroString(oss.str());
}

auto
save(const Geometry &geom, std::ostream &out) -> void
{
  std::vector<Triangle> all_triangles;

  std::function<void(const Geometry &)> process_geometry =
      [&](const Geometry &geom) {
        switch (geom.geometryTypeId()) {
        case TYPE_TRIANGLE: {
          all_triangles.push_back(geom.as<Triangle>());
          break;
        }
        case TYPE_POLYGON: {
          const auto         &poly = geom.as<Polygon>();
          TriangulatedSurface tin;
          triangulate::triangulatePolygon3D(poly, tin);
          for (size_t i = 0; i < tin.numTriangles(); ++i) {
            all_triangles.push_back(tin.triangleN(i));
          }
          break;
        }
        case TYPE_TRIANGULATEDSURFACE: {
          const auto &tin = geom.as<TriangulatedSurface>();
          for (size_t i = 0; i < tin.numTriangles(); ++i) {
            all_triangles.push_back(tin.triangleN(i));
          }
          break;
        }
        case TYPE_POLYHEDRALSURFACE: {
          const auto &phs = geom.as<PolyhedralSurface>();
          for (size_t i = 0; i < phs.numPolygons(); ++i) {
            process_geometry(phs.polygonN(i));
          }
          break;
        }
        case TYPE_SOLID: {
          const auto &solid = geom.as<Solid>();
          if (!solid.isEmpty()) {
            process_geometry(solid.exteriorShell());
          }
          break;
        }
        case TYPE_MULTIPOLYGON:
        case TYPE_MULTISOLID:
        case TYPE_GEOMETRYCOLLECTION: {
          const auto &geomcoll = geom.as<GeometryCollection>();
          for (size_t i = 0; i < geomcoll.numGeometries(); ++i) {
            process_geometry(geomcoll.geometryN(i));
          }
          break;
        }
        default:
          // Ignore other geometry types as they can't be represented in STL
          break;
        }
      };

  process_geometry(geom);

  // Write STL header
  out << "solid SFCGAL_export\n";

  // Write triangles
  for (const auto &triangle : all_triangles) {
    CGAL::Vector_3<Kernel> normal = CGAL::normal(
        triangle.vertex(0).toPoint_3(), triangle.vertex(1).toPoint_3(),
        triangle.vertex(2).toPoint_3());

    // Get normalized string values to avoid -0 output
    std::string normalizedX = normalizeKernelValue(normal.x());
    std::string normalizedY = normalizeKernelValue(normal.y());
    std::string normalizedZ = normalizeKernelValue(normal.z());

    out << "  facet normal " << normalizedX << " " << normalizedY << " "
        << normalizedZ << "\n";
    out << "    outer loop\n";
    for (int i = 0; i < 3; ++i) {
      const auto &vertex  = triangle.vertex(i);
      std::string vertexX = normalizeKernelValue(vertex.x());
      std::string vertexY = normalizeKernelValue(vertex.y());
      std::string vertexZ =
          vertex.is3D() ? normalizeKernelValue(vertex.z()) : "0";
      out << "      vertex " << vertexX << " " << vertexY << " " << vertexZ
          << "\n";
    }
    out << "    endloop\n";
    out << "  endfacet\n";
  }

  // Write STL footer
  out << "endsolid SFCGAL_export\n";
}

auto
save(const Geometry &geom, const std::string &filename) -> void
{
  std::ofstream out(filename);
  if (!out) {
    BOOST_THROW_EXCEPTION(
        Exception("Unable to open file " + filename + " for writing."));
  }
  save(geom, out);
}

auto
saveToString(const Geometry &geom) -> std::string
{
  std::ostringstream oss;
  save(geom, oss);
  return oss.str();
}

auto
saveToBuffer(const Geometry &geom, char *buffer, size_t *size) -> void
{
  std::string result = saveToString(geom);
  // Need space for content + null terminator
  size_t requiredSize = result.size() + 1;

  if ((buffer != nullptr) && *size >= requiredSize) {
    std::copy(result.begin(), result.end(), buffer);
    buffer[result.size()] = '\0'; // Null terminate
    *size                 = requiredSize;
  } else {
    *size = requiredSize;
  }
}

namespace {

/**
 * @brief Read and validate a keyword (case-insensitive)
 */
auto
expectKeyword(std::istream &in, const std::string &expected) -> bool
{
  std::string word;
  if (!(in >> word)) {
    return false;
  }
  return detail::caseInsensitiveEqual(word, expected);
}

} // anonymous namespace

auto
load(std::istream &in) -> std::unique_ptr<Geometry>
{
  std::vector<Triangle> triangles;

  // Skip whitespace and check for empty file
  detail::skipWhitespace(in);
  if (!in || in.peek() == EOF) {
    BOOST_THROW_EXCEPTION(Exception("Empty STL file"));
  }

  // Read "solid" keyword
  if (!expectKeyword(in, "solid")) {
    BOOST_THROW_EXCEPTION(Exception("STL file must start with 'solid'"));
  }

  // Read optional solid name (rest of line)
  std::string solidName;
  std::getline(in, solidName);

  // Parse facets
  while (in) {
    detail::skipWhitespace(in);
    if (!in || in.peek() == EOF) {
      break;
    }

    std::string keyword;
    if (!(in >> keyword)) {
      break;
    }

    // Convert to lowercase for comparison
    std::string keywordLower = keyword;
    std::transform(keywordLower.begin(), keywordLower.end(),
                   keywordLower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (keywordLower == "endsolid") {
      // End of solid - success
      break;
    }

    if (keywordLower != "facet") {
      BOOST_THROW_EXCEPTION(
          Exception("Expected 'facet' or 'endsolid', got: " + keyword));
    }

    // Parse "normal nx ny nz"
    if (!expectKeyword(in, "normal")) {
      BOOST_THROW_EXCEPTION(Exception("Expected 'normal' after 'facet'"));
    }

    // Read normal (we don't use it for geometry, but validate it)
    detail::parseDouble(in, "facet normal X");
    detail::parseDouble(in, "facet normal Y");
    detail::parseDouble(in, "facet normal Z");

    // Parse "outer loop"
    if (!expectKeyword(in, "outer")) {
      BOOST_THROW_EXCEPTION(Exception("Expected 'outer loop'"));
    }
    if (!expectKeyword(in, "loop")) {
      BOOST_THROW_EXCEPTION(Exception("Expected 'loop' after 'outer'"));
    }

    // Parse 3 vertices
    std::vector<Point> vertices;
    vertices.reserve(3);
    for (int i = 0; i < 3; ++i) {
      if (!expectKeyword(in, "vertex")) {
        BOOST_THROW_EXCEPTION(
            Exception("Expected 'vertex', triangle must have 3 vertices"));
      }
      double x = detail::parseDouble(in, "vertex X coordinate");
      double y = detail::parseDouble(in, "vertex Y coordinate");
      double z = detail::parseDouble(in, "vertex Z coordinate");
      vertices.emplace_back(x, y, z);
    }

    // Parse "endloop"
    if (!expectKeyword(in, "endloop")) {
      BOOST_THROW_EXCEPTION(Exception("Expected 'endloop'"));
    }

    // Parse "endfacet"
    if (!expectKeyword(in, "endfacet")) {
      BOOST_THROW_EXCEPTION(Exception("Expected 'endfacet'"));
    }

    // Create triangle
    triangles.emplace_back(vertices[0], vertices[1], vertices[2]);
  }

  if (triangles.empty()) {
    BOOST_THROW_EXCEPTION(Exception("No triangles found in STL file"));
  }

  // Create TriangulatedSurface
  auto tin = std::make_unique<TriangulatedSurface>();
  for (auto &triangle : triangles) {
    tin->addTriangle(triangle);
  }

  return tin;
}

auto
load(const std::string &stl) -> std::unique_ptr<Geometry>
{
  std::istringstream iss(stl);
  return load(iss);
}

auto
loadFromFile(const std::string &filename) -> std::unique_ptr<Geometry>
{
  std::ifstream in(filename);
  if (!in) {
    BOOST_THROW_EXCEPTION(
        Exception("Unable to open file " + filename + " for reading."));
  }
  return load(in);
}

} // namespace SFCGAL::io::STL
