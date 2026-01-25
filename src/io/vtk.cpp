// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/vtk.h"
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
#include "SFCGAL/io/io_utils.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <vector>

namespace SFCGAL::io::VTK {

// VTK cell type constants
// Linear cells
constexpr int VTK_VERTEX         = 1;
constexpr int VTK_POLY_VERTEX    = 2;
constexpr int VTK_LINE           = 3;
constexpr int VTK_POLY_LINE      = 4;
constexpr int VTK_TRIANGLE       = 5;
constexpr int VTK_TRIANGLE_STRIP = 6;
constexpr int VTK_POLYGON        = 7;
constexpr int VTK_PIXEL          = 8;
constexpr int VTK_QUAD           = 9;
constexpr int VTK_TETRA          = 10;
constexpr int VTK_VOXEL          = 11;
constexpr int VTK_HEXAHEDRON     = 12;
constexpr int VTK_WEDGE          = 13;
constexpr int VTK_PYRAMID        = 14;

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

  if (dataset != "DATASET") {
    BOOST_THROW_EXCEPTION(Exception("Expected DATASET keyword"));
  }

  bool isPolyData         = (datasetType == "POLYDATA");
  bool isUnstructuredGrid = (datasetType == "UNSTRUCTURED_GRID");
  bool isStructuredPoints = (datasetType == "STRUCTURED_POINTS");

  if (!isPolyData && !isUnstructuredGrid && !isStructuredPoints) {
    BOOST_THROW_EXCEPTION(
        Exception("Only POLYDATA, UNSTRUCTURED_GRID, and STRUCTURED_POINTS "
                  "dataset types are supported, got: " +
                  datasetType));
  }

  // For STRUCTURED_POINTS: parse grid parameters and generate points
  if (isStructuredPoints) {
    size_t dimX = 0, dimY = 0, dimZ = 0;
    double originX = 0, originY = 0, originZ = 0;
    double spacingX = 1, spacingY = 1, spacingZ = 1;

    while (in) {
      std::string keyword;
      if (!(in >> keyword)) {
        break;
      }
      std::transform(keyword.begin(), keyword.end(), keyword.begin(),
                     [](unsigned char c) { return std::toupper(c); });

      if (keyword == "DIMENSIONS") {
        in >> dimX >> dimY >> dimZ;
      } else if (keyword == "ORIGIN") {
        in >> originX >> originY >> originZ;
      } else if (keyword == "SPACING" || keyword == "ASPECT_RATIO") {
        in >> spacingX >> spacingY >> spacingZ;
      } else if (keyword == "POINT_DATA" || keyword == "CELL_DATA") {
        // Skip data sections - we only need geometry
        break;
      }
    }

    if (dimX == 0 || dimY == 0 || dimZ == 0) {
      BOOST_THROW_EXCEPTION(
          Exception("STRUCTURED_POINTS requires valid DIMENSIONS"));
    }

    // Security: Validate dimensions to prevent overflow/OOM
    constexpr size_t MAX_DIM_PER_AXIS = 10000;
    constexpr size_t MAX_GRID_POINTS  = 100'000'000; // 100M points max

    if (dimX > MAX_DIM_PER_AXIS || dimY > MAX_DIM_PER_AXIS ||
        dimZ > MAX_DIM_PER_AXIS) {
      BOOST_THROW_EXCEPTION(Exception(
          "STRUCTURED_POINTS dimensions too large (max 10000 per axis)"));
    }

    size_t totalPoints = dimX * dimY * dimZ;
    if (totalPoints > MAX_GRID_POINTS) {
      BOOST_THROW_EXCEPTION(
          Exception("STRUCTURED_POINTS total points exceeds maximum (100M)"));
    }

    // Generate grid points
    data.points.reserve(totalPoints);
    for (size_t k = 0; k < dimZ; ++k) {
      for (size_t j = 0; j < dimY; ++j) {
        for (size_t i = 0; i < dimX; ++i) {
          double x = originX + static_cast<double>(i) * spacingX;
          double y = originY + static_cast<double>(j) * spacingY;
          double z = originZ + static_cast<double>(k) * spacingZ;
          data.points.emplace_back(x, y, z);
        }
      }
    }

    // Create implicit voxel cells for the grid (dimX-1) * (dimY-1) * (dimZ-1)
    if (dimX > 1 && dimY > 1 && dimZ > 1) {
      for (size_t k = 0; k < dimZ - 1; ++k) {
        for (size_t j = 0; j < dimY - 1; ++j) {
          for (size_t i = 0; i < dimX - 1; ++i) {
            // Voxel vertices in VTK order
            size_t idx0 = i + j * dimX + k * dimX * dimY;
            size_t idx1 = (i + 1) + j * dimX + k * dimX * dimY;
            size_t idx2 = i + (j + 1) * dimX + k * dimX * dimY;
            size_t idx3 = (i + 1) + (j + 1) * dimX + k * dimX * dimY;
            size_t idx4 = i + j * dimX + (k + 1) * dimX * dimY;
            size_t idx5 = (i + 1) + j * dimX + (k + 1) * dimX * dimY;
            size_t idx6 = i + (j + 1) * dimX + (k + 1) * dimX * dimY;
            size_t idx7 = (i + 1) + (j + 1) * dimX + (k + 1) * dimX * dimY;

            data.cells.push_back(
                {idx0, idx1, idx2, idx3, idx4, idx5, idx6, idx7});
            data.cellTypes.push_back(VTK_VOXEL);
          }
        }
      }
    }

    return data;
  }

  // Security: Maximum limits to prevent OOM attacks
  constexpr size_t MAX_CELLS     = 10'000'000;  // 10M cells max
  constexpr size_t MAX_CELL_SIZE = 1'000'000;   // 1M vertices per cell max
  constexpr size_t MAX_POINTS    = 100'000'000; // 100M points max

  // Helper lambda to read cells with size prefix
  auto readCellsWithSize = [&data, &in](size_t numCells) {
    if (numCells > MAX_CELLS) {
      BOOST_THROW_EXCEPTION(
          Exception("Number of cells exceeds maximum (10M): " +
                    std::to_string(numCells)));
    }

    for (size_t i = 0; i < numCells; ++i) {
      size_t cellSize = 0;
      in >> cellSize;

      if (cellSize > MAX_CELL_SIZE) {
        BOOST_THROW_EXCEPTION(Exception("Cell size exceeds maximum (1M): " +
                                        std::to_string(cellSize)));
      }

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
  };

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

      if (numPoints > MAX_POINTS) {
        BOOST_THROW_EXCEPTION(
            Exception("Number of points exceeds maximum (100M): " +
                      std::to_string(numPoints)));
      }

      data.points.reserve(numPoints);
      for (size_t i = 0; i < numPoints; ++i) {
        double x = detail::parseDouble(in, "POINTS X coordinate");
        double y = detail::parseDouble(in, "POINTS Y coordinate");
        double z = detail::parseDouble(in, "POINTS Z coordinate");
        data.points.emplace_back(x, y, z);
      }
    } else if (section == "CELLS") {
      // UNSTRUCTURED_GRID format
      size_t numCells  = 0;
      size_t totalSize = 0;
      in >> numCells >> totalSize;
      readCellsWithSize(numCells);
    } else if (section == "VERTICES") {
      // POLYDATA format - vertices
      size_t numCells  = 0;
      size_t totalSize = 0;
      in >> numCells >> totalSize;
      size_t startIdx = data.cells.size();
      readCellsWithSize(numCells);
      // Assign VTK_VERTEX type
      for (size_t i = startIdx; i < data.cells.size(); ++i) {
        data.cellTypes.push_back(VTK_VERTEX);
      }
    } else if (section == "LINES") {
      // POLYDATA format - lines
      size_t numCells  = 0;
      size_t totalSize = 0;
      in >> numCells >> totalSize;
      size_t startIdx = data.cells.size();
      readCellsWithSize(numCells);
      // Assign VTK_POLY_LINE type
      for (size_t i = startIdx; i < data.cells.size(); ++i) {
        data.cellTypes.push_back(VTK_POLY_LINE);
      }
    } else if (section == "POLYGONS") {
      // POLYDATA format - polygons
      size_t numCells  = 0;
      size_t totalSize = 0;
      in >> numCells >> totalSize;
      size_t startIdx = data.cells.size();
      readCellsWithSize(numCells);
      // Assign VTK_POLYGON type
      for (size_t i = startIdx; i < data.cells.size(); ++i) {
        data.cellTypes.push_back(VTK_POLYGON);
      }
    } else if (section == "TRIANGLE_STRIPS") {
      // POLYDATA format - triangle strips
      size_t numCells  = 0;
      size_t totalSize = 0;
      in >> numCells >> totalSize;
      size_t startIdx = data.cells.size();
      readCellsWithSize(numCells);
      // Assign VTK_TRIANGLE_STRIP type
      for (size_t i = startIdx; i < data.cells.size(); ++i) {
        data.cellTypes.push_back(VTK_TRIANGLE_STRIP);
      }
    } else if (section == "CELL_TYPES") {
      // UNSTRUCTURED_GRID format - cell types
      size_t numTypes = 0;
      in >> numTypes;

      if (numTypes != data.cells.size()) {
        BOOST_THROW_EXCEPTION(
            Exception("CELL_TYPES count (" + std::to_string(numTypes) +
                      ") doesn't match CELLS count (" +
                      std::to_string(data.cells.size()) + ")"));
      }

      data.cellTypes.clear();
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
// NOLINTBEGIN(readability-function-cognitive-complexity)
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

    case VTK_POLY_VERTEX:
      // Multiple vertices in one cell
      for (size_t idx : cell) {
        vertexPoints.push_back(data.points[idx]);
      }
      break;

    case VTK_LINE: {
      // Simple line with 2 points
      if (cell.size() >= 2) {
        auto ls = std::make_unique<LineString>();
        ls->addPoint(data.points[cell[0]]);
        ls->addPoint(data.points[cell[1]]);
        lines.push_back(std::move(ls));
      }
      break;
    }

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

    case VTK_TRIANGLE_STRIP: {
      // Triangle strip: vertices v0, v1, v2, v3, ... form triangles
      // (v0,v1,v2), (v2,v1,v3), (v2,v3,v4), ...
      if (cell.size() >= 3) {
        for (size_t j = 0; j + 2 < cell.size(); ++j) {
          if (j % 2 == 0) {
            triangles.push_back(std::make_unique<Triangle>(
                data.points[cell[j]], data.points[cell[j + 1]],
                data.points[cell[j + 2]]));
          } else {
            // Reverse winding for odd triangles
            triangles.push_back(std::make_unique<Triangle>(
                data.points[cell[j + 1]], data.points[cell[j]],
                data.points[cell[j + 2]]));
          }
        }
      }
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

    case VTK_QUAD: {
      // Quadrilateral - 4 points
      if (cell.size() != 4) {
        BOOST_THROW_EXCEPTION(Exception("VTK_QUAD must have exactly 4 points"));
      }
      auto ring = std::make_unique<LineString>();
      for (size_t idx : cell) {
        ring->addPoint(data.points[idx]);
      }
      ring->addPoint(data.points[cell[0]]); // close
      polygons.push_back(std::make_unique<Polygon>(ring.release()));
      break;
    }

    case VTK_PIXEL: {
      // Pixel - 4 points with different ordering than QUAD
      // VTK_PIXEL order: (i,j), (i+1,j), (i,j+1), (i+1,j+1)
      // Need to reorder to: 0, 1, 3, 2 for proper polygon winding
      if (cell.size() != 4) {
        BOOST_THROW_EXCEPTION(
            Exception("VTK_PIXEL must have exactly 4 points"));
      }
      auto ring = std::make_unique<LineString>();
      ring->addPoint(data.points[cell[0]]);
      ring->addPoint(data.points[cell[1]]);
      ring->addPoint(data.points[cell[3]]);
      ring->addPoint(data.points[cell[2]]);
      ring->addPoint(data.points[cell[0]]); // close
      polygons.push_back(std::make_unique<Polygon>(ring.release()));
      break;
    }

    case VTK_TETRA: {
      // Tetrahedron - 4 triangular faces
      if (cell.size() != 4) {
        BOOST_THROW_EXCEPTION(
            Exception("VTK_TETRA must have exactly 4 points"));
      }
      // Create 4 triangular faces with outward normals
      triangles.push_back(std::make_unique<Triangle>(
          data.points[cell[0]], data.points[cell[2]], data.points[cell[1]]));
      triangles.push_back(std::make_unique<Triangle>(
          data.points[cell[0]], data.points[cell[1]], data.points[cell[3]]));
      triangles.push_back(std::make_unique<Triangle>(
          data.points[cell[0]], data.points[cell[3]], data.points[cell[2]]));
      triangles.push_back(std::make_unique<Triangle>(
          data.points[cell[1]], data.points[cell[2]], data.points[cell[3]]));
      break;
    }

    case VTK_VOXEL: {
      // Voxel - 8 points with specific ordering
      // VTK order: (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k),
      //            (i,j,k+1), (i+1,j,k+1), (i,j+1,k+1), (i+1,j+1,k+1)
      if (cell.size() != 8) {
        BOOST_THROW_EXCEPTION(
            Exception("VTK_VOXEL must have exactly 8 points"));
      }
      // 6 quadrilateral faces (reordered for proper winding)
      // Bottom face (z=0): 0,1,3,2
      auto addQuadFace = [&](size_t a, size_t b, size_t c, size_t d) {
        auto ring = std::make_unique<LineString>();
        ring->addPoint(data.points[cell[a]]);
        ring->addPoint(data.points[cell[b]]);
        ring->addPoint(data.points[cell[c]]);
        ring->addPoint(data.points[cell[d]]);
        ring->addPoint(data.points[cell[a]]);
        polygons.push_back(std::make_unique<Polygon>(ring.release()));
      };
      addQuadFace(0, 2, 3, 1); // bottom
      addQuadFace(4, 5, 7, 6); // top
      addQuadFace(0, 1, 5, 4); // front
      addQuadFace(2, 6, 7, 3); // back
      addQuadFace(0, 4, 6, 2); // left
      addQuadFace(1, 3, 7, 5); // right
      break;
    }

    case VTK_HEXAHEDRON: {
      // Hexahedron - 8 points
      // VTK order: 0-3 bottom face, 4-7 top face (aligned with bottom)
      if (cell.size() != 8) {
        BOOST_THROW_EXCEPTION(
            Exception("VTK_HEXAHEDRON must have exactly 8 points"));
      }
      // 6 quadrilateral faces
      auto addQuadFace = [&](size_t a, size_t b, size_t c, size_t d) {
        auto ring = std::make_unique<LineString>();
        ring->addPoint(data.points[cell[a]]);
        ring->addPoint(data.points[cell[b]]);
        ring->addPoint(data.points[cell[c]]);
        ring->addPoint(data.points[cell[d]]);
        ring->addPoint(data.points[cell[a]]);
        polygons.push_back(std::make_unique<Polygon>(ring.release()));
      };
      addQuadFace(0, 3, 2, 1); // bottom
      addQuadFace(4, 5, 6, 7); // top
      addQuadFace(0, 1, 5, 4); // front
      addQuadFace(2, 3, 7, 6); // back
      addQuadFace(0, 4, 7, 3); // left
      addQuadFace(1, 2, 6, 5); // right
      break;
    }

    case VTK_WEDGE: {
      // Wedge (triangular prism) - 6 points
      // VTK order: 0-2 bottom triangle, 3-5 top triangle
      if (cell.size() != 6) {
        BOOST_THROW_EXCEPTION(
            Exception("VTK_WEDGE must have exactly 6 points"));
      }
      // 2 triangular faces + 3 quadrilateral faces
      triangles.push_back(std::make_unique<Triangle>(
          data.points[cell[0]], data.points[cell[2]], data.points[cell[1]]));
      triangles.push_back(std::make_unique<Triangle>(
          data.points[cell[3]], data.points[cell[4]], data.points[cell[5]]));
      auto addQuadFace = [&](size_t a, size_t b, size_t c, size_t d) {
        auto ring = std::make_unique<LineString>();
        ring->addPoint(data.points[cell[a]]);
        ring->addPoint(data.points[cell[b]]);
        ring->addPoint(data.points[cell[c]]);
        ring->addPoint(data.points[cell[d]]);
        ring->addPoint(data.points[cell[a]]);
        polygons.push_back(std::make_unique<Polygon>(ring.release()));
      };
      addQuadFace(0, 1, 4, 3);
      addQuadFace(1, 2, 5, 4);
      addQuadFace(2, 0, 3, 5);
      break;
    }

    case VTK_PYRAMID: {
      // Pyramid - 5 points
      // VTK order: 0-3 base quad, 4 apex
      if (cell.size() != 5) {
        BOOST_THROW_EXCEPTION(
            Exception("VTK_PYRAMID must have exactly 5 points"));
      }
      // 1 quadrilateral face (base) + 4 triangular faces
      auto ring = std::make_unique<LineString>();
      ring->addPoint(data.points[cell[0]]);
      ring->addPoint(data.points[cell[3]]);
      ring->addPoint(data.points[cell[2]]);
      ring->addPoint(data.points[cell[1]]);
      ring->addPoint(data.points[cell[0]]);
      polygons.push_back(std::make_unique<Polygon>(ring.release()));
      triangles.push_back(std::make_unique<Triangle>(
          data.points[cell[0]], data.points[cell[1]], data.points[cell[4]]));
      triangles.push_back(std::make_unique<Triangle>(
          data.points[cell[1]], data.points[cell[2]], data.points[cell[4]]));
      triangles.push_back(std::make_unique<Triangle>(
          data.points[cell[2]], data.points[cell[3]], data.points[cell[4]]));
      triangles.push_back(std::make_unique<Triangle>(
          data.points[cell[3]], data.points[cell[0]], data.points[cell[4]]));
      break;
    }

    default:
      // Skip unknown cell types
      break;
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
      bool allTriangles = std::all_of(polygons.begin(), polygons.end(),
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
// NOLINTEND(readability-function-cognitive-complexity)

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
