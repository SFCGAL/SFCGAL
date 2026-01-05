// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/vtk.h"
#include "SFCGAL/io/io_utils.h"
#include "SFCGAL/Exception.h"
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
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <vector>

namespace SFCGAL::io::VTK {

// VTK cell type constants
constexpr int VTK_VERTEX    = 1;
constexpr int VTK_POLY_LINE = 4;
constexpr int VTK_TRIANGLE  = 5;
constexpr int VTK_POLYGON   = 7;

// NOLINTBEGIN(readability-function-cognitive-complexity)
void
save(const Geometry &geom, std::ostream &out)
{
  std::vector<Point>               all_points;
  std::vector<std::vector<size_t>> all_cells;
  std::vector<int>                 cell_types;

  std::function<void(const Geometry &)> process_geometry =
      [&](const Geometry &geometry) {
        switch (geometry.geometryTypeId()) {
        case TYPE_POINT: {
          const auto &point = geometry.as<Point>();
          all_points.push_back(point);
          all_cells.push_back({all_points.size() - 1});
          cell_types.push_back(VTK_VERTEX);
          break;
        }
        case TYPE_LINESTRING: {
          const auto         &linestring = geometry.as<LineString>();
          std::vector<size_t> line;
          for (size_t i = 0; i < linestring.numPoints(); ++i) {
            all_points.push_back(linestring.pointN(i));
            line.push_back(all_points.size() - 1);
          }
          all_cells.push_back(line);
          cell_types.push_back(VTK_POLY_LINE);
          break;
        }
        case TYPE_TRIANGLE: {
          const auto         &tri = geometry.as<Triangle>();
          std::vector<size_t> face;
          for (int i = 0; i < 3; ++i) {
            all_points.push_back(tri.vertex(i));
            face.push_back(all_points.size() - 1);
          }
          all_cells.push_back(face);
          cell_types.push_back(VTK_TRIANGLE);
          break;
        }
        case TYPE_POLYGON: {
          const auto         &poly = geometry.as<Polygon>();
          std::vector<size_t> face;
          for (size_t i = 0; i < poly.exteriorRing().numPoints() - 1; ++i) {
            all_points.push_back(poly.exteriorRing().pointN(i));
            face.push_back(all_points.size() - 1);
          }
          all_cells.push_back(face);
          cell_types.push_back(VTK_POLYGON);
          break;
        }
        case TYPE_TRIANGULATEDSURFACE: {
          const auto &triangulatedsurface = geometry.as<TriangulatedSurface>();
          for (size_t i = 0; i < triangulatedsurface.numPatches(); ++i) {
            process_geometry(triangulatedsurface.patchN(i));
          }
          break;
        }
        case TYPE_POLYHEDRALSURFACE: {
          const auto &polyhedralsurface = geometry.as<PolyhedralSurface>();
          for (size_t i = 0; i < polyhedralsurface.numPatches(); ++i) {
            process_geometry(polyhedralsurface.patchN(i));
          }
          break;
        }
        case TYPE_SOLID: {
          const auto &solid = geometry.as<Solid>();
          if (!solid.isEmpty()) {
            process_geometry(solid.exteriorShell());
          }
          break;
        }
        case TYPE_MULTIPOINT:
        case TYPE_MULTILINESTRING:
        case TYPE_MULTIPOLYGON:
        case TYPE_MULTISOLID:
        case TYPE_GEOMETRYCOLLECTION: {
          const auto &geometrycollection = geometry.as<GeometryCollection>();
          for (size_t i = 0; i < geometrycollection.numGeometries(); ++i) {
            process_geometry(geometrycollection.geometryN(i));
          }
          break;
        }
        default:
          BOOST_THROW_EXCEPTION(InappropriateGeometryException(
              "Unsupported geometry type: " + geometry.geometryType()));
        }
      };

  process_geometry(geom);

  // Write VTK header
  out << "# vtk DataFile Version 2.0\n";
  out << "SFCGAL Geometry\n";
  out << "ASCII\n";
  out << "DATASET UNSTRUCTURED_GRID\n";

  // Write points
  out << "POINTS " << all_points.size() << " float\n";
  for (const auto &point : all_points) {
    out << point.x() << " " << point.y() << " "
        << (point.is3D() ? point.z() : 0.0) << "\n";
  }

  // Write cells
  size_t total_cell_size = 0;
  for (const auto &cell : all_cells) {
    total_cell_size += cell.size() + 1; // +1 for the size prefix
  }
  out << "CELLS " << all_cells.size() << " " << total_cell_size << "\n";
  for (const auto &cell : all_cells) {
    out << cell.size();
    for (size_t idx : cell) {
      out << " " << idx;
    }
    out << "\n";
  }

  // Write cell types
  out << "CELL_TYPES " << cell_types.size() << "\n";
  for (int type : cell_types) {
    out << type << "\n";
  }
}
// NOLINTEND(readability-function-cognitive-complexity)

void
save(const Geometry &geom, const std::string &filename)
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

void
saveToBuffer(const Geometry &geom, char *buffer, size_t *size)
{
  std::string result = saveToString(geom);
  if ((buffer != nullptr) && *size >= result.size()) {
    std::copy(result.begin(), result.end(), buffer);
    *size = result.size();
  } else {
    *size = result.size();
  }
}

namespace {

/**
 * @brief Parsed VTK data structure
 */
struct VtkData {
  std::vector<Point>               points;
  std::vector<std::vector<size_t>> cells;
  std::vector<int>                 cellTypes;
};

/**
 * @brief Parse VTK data from input stream
 */
auto
parseVtkData(std::istream &in) -> VtkData
{
  VtkData data;

  // Read VTK header line (# vtk DataFile Version X.X)
  std::string line;
  std::getline(in, line);
  // Convert to lowercase for case-insensitive comparison
  std::string lineLower = line;
  std::transform(lineLower.begin(), lineLower.end(), lineLower.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (lineLower.find("vtk") == std::string::npos) {
    BOOST_THROW_EXCEPTION(Exception("Invalid VTK header: expected '# vtk'"));
  }

  // Read title line
  std::getline(in, line);

  // Read format (must be ASCII)
  std::string format;
  in >> format;
  std::transform(format.begin(), format.end(), format.begin(),
                 [](unsigned char c) { return std::toupper(c); });
  if (format != "ASCII") {
    BOOST_THROW_EXCEPTION(Exception("Only ASCII VTK files are supported"));
  }

  // Read dataset type
  std::string dataset, datasetType;
  in >> dataset >> datasetType;
  std::transform(dataset.begin(), dataset.end(), dataset.begin(),
                 [](unsigned char c) { return std::toupper(c); });
  std::transform(datasetType.begin(), datasetType.end(), datasetType.begin(),
                 [](unsigned char c) { return std::toupper(c); });

  if (dataset != "DATASET" || datasetType != "UNSTRUCTURED_GRID") {
    BOOST_THROW_EXCEPTION(
        Exception("Only UNSTRUCTURED_GRID dataset type is supported"));
  }

  // Parse sections
  while (in) {
    std::string section;
    if (!(in >> section)) {
      break;
    }

    std::transform(section.begin(), section.end(), section.begin(),
                   [](unsigned char c) { return std::toupper(c); });

    if (section == "POINTS") {
      size_t      numPoints = 0;
      std::string dataType;
      in >> numPoints >> dataType;

      data.points.reserve(numPoints);
      for (size_t i = 0; i < numPoints; ++i) {
        double x = detail::parseDouble(in, "POINTS X coordinate");
        double y = detail::parseDouble(in, "POINTS Y coordinate");
        double z = detail::parseDouble(in, "POINTS Z coordinate");
        data.points.emplace_back(x, y, z);
      }
    } else if (section == "CELLS") {
      size_t numCells  = 0;
      size_t totalSize = 0;
      in >> numCells >> totalSize;

      data.cells.reserve(numCells);
      for (size_t i = 0; i < numCells; ++i) {
        size_t cellSize = 0;
        in >> cellSize;

        std::vector<size_t> cell;
        cell.reserve(cellSize);
        for (size_t j = 0; j < cellSize; ++j) {
          size_t idx = 0;
          in >> idx;

          if (idx >= data.points.size()) {
            BOOST_THROW_EXCEPTION(Exception(
                "Cell references invalid point index: " + std::to_string(idx) +
                " (only " + std::to_string(data.points.size()) +
                " points available)"));
          }
          cell.push_back(idx);
        }
        data.cells.push_back(std::move(cell));
      }
    } else if (section == "CELL_TYPES") {
      size_t numTypes = 0;
      in >> numTypes;

      if (numTypes != data.cells.size()) {
        BOOST_THROW_EXCEPTION(
            Exception("CELL_TYPES count (" + std::to_string(numTypes) +
                      ") doesn't match CELLS count (" +
                      std::to_string(data.cells.size()) + ")"));
      }

      data.cellTypes.reserve(numTypes);
      for (size_t i = 0; i < numTypes; ++i) {
        int cellType = 0;
        in >> cellType;
        data.cellTypes.push_back(cellType);
      }
    }
    // Ignore other sections (POINT_DATA, CELL_DATA, etc.)
  }

  return data;
}

/**
 * @brief Create geometry from parsed VTK data
 */
auto
createGeometryFromVtkData(const VtkData &data) -> std::unique_ptr<Geometry>
{
  if (data.cells.empty()) {
    BOOST_THROW_EXCEPTION(Exception("No cells found in VTK file"));
  }

  if (data.cellTypes.size() != data.cells.size()) {
    BOOST_THROW_EXCEPTION(Exception("Cell types count mismatch"));
  }

  // Categorize cells by type
  std::vector<Point>                       vertexPoints;
  std::vector<std::unique_ptr<LineString>> lines;
  std::vector<std::unique_ptr<Triangle>>   triangles;
  std::vector<std::unique_ptr<Polygon>>    polygons;

  for (size_t i = 0; i < data.cells.size(); ++i) {
    const auto &cell     = data.cells[i];
    int         cellType = data.cellTypes[i];

    switch (cellType) {
    case VTK_VERTEX:
      if (!cell.empty()) {
        vertexPoints.push_back(data.points[cell[0]]);
      }
      break;

    case VTK_POLY_LINE: {
      auto ls = std::make_unique<LineString>();
      for (size_t idx : cell) {
        ls->addPoint(data.points[idx]);
      }
      lines.push_back(std::move(ls));
      break;
    }

    case VTK_TRIANGLE: {
      if (cell.size() != 3) {
        BOOST_THROW_EXCEPTION(
            Exception("VTK_TRIANGLE must have exactly 3 points"));
      }
      triangles.push_back(std::make_unique<Triangle>(
          data.points[cell[0]], data.points[cell[1]], data.points[cell[2]]));
      break;
    }

    case VTK_POLYGON: {
      auto ring = std::make_unique<LineString>();
      for (size_t idx : cell) {
        ring->addPoint(data.points[idx]);
      }
      // Close the ring
      if (!cell.empty()) {
        ring->addPoint(data.points[cell[0]]);
      }
      polygons.push_back(std::make_unique<Polygon>(ring.release()));
      break;
    }

    default:
      BOOST_THROW_EXCEPTION(
          Exception("Unsupported VTK cell type: " + std::to_string(cellType)));
    }
  }

  // Determine output geometry type based on what we found
  // Priority: surfaces > lines > points

  // If we have triangles and/or polygons, create a surface
  if (!triangles.empty() || !polygons.empty()) {
    // Check if all faces are triangles
    if (polygons.empty()) {
      // All VTK_TRIANGLE cells - create TriangulatedSurface
      auto tin = std::make_unique<TriangulatedSurface>();
      for (auto &tri : triangles) {
        tin->addTriangle(*tri);
      }
      return tin;
    }
    if (triangles.empty()) {
      // All polygons (no VTK_TRIANGLE) - check if they're all triangular
      bool allTriangles =
          std::all_of(polygons.begin(), polygons.end(),
                      [](const std::unique_ptr<Polygon> &p) {
                        return p->exteriorRing().numPoints() ==
                               4; // 3 vertices + close
                      });

      if (allTriangles) {
        auto tin = std::make_unique<TriangulatedSurface>();
        for (auto &poly : polygons) {
          const auto &ring = poly->exteriorRing();
          tin->addTriangle(
              Triangle(ring.pointN(0), ring.pointN(1), ring.pointN(2)));
        }
        return tin;
      }
      // Mixed polygon sizes - create PolyhedralSurface
      auto phs = std::make_unique<PolyhedralSurface>();
      for (auto &poly : polygons) {
        phs->addPolygon(*poly);
      }
      return phs;
    }
    // Mix of VTK_TRIANGLE and VTK_POLYGON - create PolyhedralSurface
    auto phs = std::make_unique<PolyhedralSurface>();
    for (auto &tri : triangles) {
      // Convert triangle to polygon
      auto ring = std::make_unique<LineString>();
      ring->addPoint(tri->vertex(0));
      ring->addPoint(tri->vertex(1));
      ring->addPoint(tri->vertex(2));
      ring->addPoint(tri->vertex(0)); // close
      phs->addPolygon(Polygon(ring.release()));
    }
    for (auto &poly : polygons) {
      phs->addPolygon(*poly);
    }
    return phs;
  }

  // If we have lines, create MultiLineString
  if (!lines.empty()) {
    auto mls = std::make_unique<MultiLineString>();
    for (auto &ls : lines) {
      mls->addGeometry(ls.release());
    }
    return mls;
  }

  // If we have points, create MultiPoint
  if (!vertexPoints.empty()) {
    auto mp = std::make_unique<MultiPoint>();
    for (const auto &pt : vertexPoints) {
      mp->addGeometry(new Point(pt));
    }
    return mp;
  }

  BOOST_THROW_EXCEPTION(Exception("No valid geometry found in VTK file"));
}

} // anonymous namespace

auto
load(std::istream &in) -> std::unique_ptr<Geometry>
{
  VtkData data = parseVtkData(in);
  return createGeometryFromVtkData(data);
}

auto
load(const std::string &vtk) -> std::unique_ptr<Geometry>
{
  std::istringstream iss(vtk);
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

} // namespace SFCGAL::io::VTK
