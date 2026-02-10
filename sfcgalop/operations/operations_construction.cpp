// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations_construction.hpp"

#include "SFCGAL/Envelope.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPolygon.h"
#ifndef _MSC_VER
  #include "SFCGAL/algorithm/alphaShapes.h"
#endif
#include "SFCGAL/algorithm/Chamfer.h"
#include "SFCGAL/algorithm/alphaWrapping3D.h"
#include "SFCGAL/algorithm/buffer3D.h"
#include "SFCGAL/algorithm/centroid.h"
#include "SFCGAL/algorithm/convexHull.h"
#include "SFCGAL/algorithm/insertPointsWithinTolerance.h"
#include "SFCGAL/algorithm/lineSubstring.h"
#include "SFCGAL/algorithm/minkowskiSum.h"
#include "SFCGAL/algorithm/minkowskiSum3D.h"
#include "SFCGAL/algorithm/offset.h"
#include "SFCGAL/algorithm/roofGeneration.h"
#include "SFCGAL/algorithm/straightSkeleton.h"
#include "SFCGAL/algorithm/sweep.h"
#include "SFCGAL/algorithm/tessellate.h"
#include "SFCGAL/triangulate/triangulate2DZ.h"

namespace Operations {

// NOLINTBEGIN(readability-function-cognitive-complexity)
/**
 * @brief Extract the angles parameter from the original string.
 *
 * Parses angles in JSON-like format: [[45,45,45,45],[30,30,30]]
 * Returns a vector of vectors of Kernel::FT values.
 *
 * @param param_str Original parameter string
 * @return Vector of vectors of angles (empty if not found or invalid)
 */
auto
extract_angles_param(const std::string &param_str)
    -> std::vector<std::vector<SFCGAL::Kernel::FT>>
{
  std::vector<std::vector<SFCGAL::Kernel::FT>> result;

  // Find "angles=" in the string
  auto pos = param_str.find("angles=");
  if (pos == std::string::npos) {
    return result; // No angles parameter
  }

  // Skip "angles="
  pos += 7;

  // Find the opening '[['
  if (pos >= param_str.size() || param_str[pos] != '[') {
    return result;
  }
  pos++; // Skip first '['

  // Parse each inner array
  while (pos < param_str.size()) {
    // Skip whitespace
    while (pos < param_str.size() && std::isspace(param_str[pos]) != 0) {
      pos++;
    }

    if (pos >= param_str.size()) {
      break;
    }

    // Check for end of outer array
    if (param_str[pos] == ']') {
      break;
    }

    // Skip comma between inner arrays
    if (param_str[pos] == ',') {
      pos++;
      continue;
    }

    // Expect start of inner array
    if (param_str[pos] != '[') {
      break;
    }
    pos++; // Skip '['

    // Parse values in inner array
    std::vector<SFCGAL::Kernel::FT> inner;
    std::string                     num_str;

    while (pos < param_str.size() && param_str[pos] != ']') {
      char chr = param_str[pos];
      if (chr == ',' || chr == ']') {
        if (!num_str.empty()) {
          inner.emplace_back(parse_double(trim(num_str)));
          num_str.clear();
        }
        if (chr == ',') {
          pos++;
        }
      } else if (std::isdigit(chr) != 0 || chr == '.' || chr == '-' ||
                 chr == '+') {
        num_str += chr;
        pos++;
      } else {
        pos++;
      }
    }

    // Don't forget the last number before ']'
    if (!num_str.empty()) {
      inner.emplace_back(parse_double(trim(num_str)));
    }

    if (pos < param_str.size() && param_str[pos] == ']') {
      pos++; // Skip closing ']'
    }

    if (!inner.empty()) {
      result.push_back(std::move(inner));
    }
  }

  return result;
}

/**
 * @brief Extract the weights parameter from the original string.
 *
 * Parses weights in JSON-like format: [[1.0,0,1.0,0],[1.0,1.0]]
 * Returns a vector of vectors of Kernel::FT values.
 *
 * @param param_str Original parameter string
 * @return Vector of vectors of weights (empty if not found or invalid)
 */
auto
extract_weights_param(const std::string &param_str)
    -> std::vector<std::vector<SFCGAL::Kernel::FT>>
{
  std::vector<std::vector<SFCGAL::Kernel::FT>> result;

  // Find "weights=" in the string
  auto pos = param_str.find("weights=");
  if (pos == std::string::npos) {
    return result; // No weights parameter
  }

  // Skip "weights="
  pos += 8;

  // Find the opening '[['
  if (pos >= param_str.size() || param_str[pos] != '[') {
    return result;
  }
  pos++; // Skip first '['

  // Parse each inner array (same logic as extract_angles_param)
  while (pos < param_str.size()) {
    // Skip whitespace
    while (pos < param_str.size() && std::isspace(param_str[pos]) != 0) {
      pos++;
    }

    if (pos >= param_str.size()) {
      break;
    }

    // Check for end of outer array
    if (param_str[pos] == ']') {
      break;
    }

    // Skip comma between inner arrays
    if (param_str[pos] == ',') {
      pos++;
      continue;
    }

    // Expect start of inner array
    if (param_str[pos] != '[') {
      break;
    }
    pos++; // Skip '['

    // Parse values in inner array
    std::vector<SFCGAL::Kernel::FT> inner;
    std::string                     num_str;

    while (pos < param_str.size() && param_str[pos] != ']') {
      char chr = param_str[pos];
      if (chr == ',' || chr == ']') {
        if (!num_str.empty()) {
          inner.emplace_back(parse_double(trim(num_str)));
          num_str.clear();
        }
        if (chr == ',') {
          pos++;
        }
      } else if (std::isdigit(chr) != 0 || chr == '.' || chr == '-' ||
                 chr == '+') {
        num_str += chr;
        pos++;
      } else {
        pos++;
      }
    }

    // Don't forget the last number before ']'
    if (!num_str.empty()) {
      inner.emplace_back(parse_double(trim(num_str)));
    }

    if (pos < param_str.size() && param_str[pos] == ']') {
      pos++; // Skip closing ']'
    }

    if (!inner.empty()) {
      result.push_back(std::move(inner));
    }
  }

  return result;
}
// NOLINTEND(readability-function-cognitive-complexity)

const std::vector<Operation> operations_construction = {
#ifndef _MSC_VER
    {"alpha_shapes", "Construction", "Compute alpha shapes from point cloud",
     false,
     "Parameters:\n  alpha=VALUE: Alpha parameter controlling shape detail "
     "(default: 1.0)\n  Smaller values create more detailed "
     "shapes\n\nExample:\n  sfcgalop -a \"MULTIPOINT((0 0),(1 0),(0.5 1),(2 "
     "0.5))\" alpha_shapes \"alpha=0.5\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params     = parse_params(args);
       double alpha      = params.count("alpha") ? params["alpha"] : 1.0;
       bool   allowHoles = true;
       return SFCGAL::algorithm::alphaShapes(*geom_a, alpha, allowHoles);
     }},
#endif

    {"alpha_wrapping_3d", "Construction",
     "Create 3D alpha wrapping surface from points", false,
     "Parameters:\n  alpha=VALUE: Alpha parameter (default: 1.0)\n  "
     "offset=VALUE: Offset parameter (default: 0)\n\nExample:\n  sfcgalop -a "
     "\"MULTIPOINT Z((0 0 0),(1 0 0),(0 1 0),(0 0 1))\" alpha_wrapping_3d "
     "\"alpha=1.5,offset=0\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double alpha  = params.count("alpha") ? params["alpha"] : 1.0;
       size_t offset =
           static_cast<size_t>(params.count("offset") ? params["offset"] : 0.0);
       auto alphaInt = static_cast<size_t>(
           alpha * 100); // Convert to integer representation
       return SFCGAL::algorithm::alphaWrapping3D(*geom_a, alphaInt, offset);
     }},

    {"boundary", "Construction",
     "Compute the topological boundary of a geometry", false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return geom_a->boundary();
     }},

    {"buffer_3d", "Construction", "Create a 3D buffer around points and lines",
     false,
     "Parameters:\n"
     "  radius=VALUE: Buffer radius (default: 1.0)\n"
     "  segments=N: Number of segments for curved surfaces (default: 16)\n"
     "  type=TYPE: Buffer type - round, cylsphere, or flat (default: round)\n\n"
     "Buffer types:\n"
     "  round     - Minkowski sum with a sphere (smooth result)\n"
     "  cylsphere - Union of cylinders and spheres (faster)\n"
     "  flat      - Construction using disk on bisector plane\n\n"
     "Only works with Point and LineString geometries.\n\n"
     "Examples:\n"
     "  sfcgalop -a \"POINT(0 0 0)\" buffer3d \"radius=2.5\"\n"
     "  sfcgalop -a \"LINESTRING Z(0 0 0,1 1 1)\" buffer3d "
     "\"radius=0.5,type=cylsphere\"\n"
     "  sfcgalop -a \"LINESTRING Z(0 0 0,1 0 0,1 1 0)\" buffer3d "
     "\"radius=0.3,segments=32,type=flat\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double radius = params.count("radius") ? params["radius"] : 1.0;
       int    segments =
           static_cast<int>(params.count("segments") ? params["segments"] : 16);

       // Parse buffer type from string
       SFCGAL::algorithm::Buffer3D::BufferType bufferType =
           SFCGAL::algorithm::Buffer3D::ROUND;

       // Check for type parameter in the args string
       if (args.find("type=") != std::string::npos) {
         size_t type_pos = args.find("type=");
         size_t start    = type_pos + 5; // length of "type="
         size_t end      = args.find(',', start);
         if (end == std::string::npos) {
           end = args.length();
         }
         std::string type_str = args.substr(start, end - start);
         // Trim and convert to lowercase
         type_str = trim(type_str);
         std::transform(
             type_str.begin(), type_str.end(), type_str.begin(),
             [](unsigned char chr) -> int { return std::tolower(chr); });

         if (type_str == "cylsphere" || type_str == "cyl" ||
             type_str == "cylinder") {
           bufferType = SFCGAL::algorithm::Buffer3D::CYLSPHERE;
         } else if (type_str == "flat" || type_str == "disk") {
           bufferType = SFCGAL::algorithm::Buffer3D::FLAT;
         }
         // else keep ROUND as default
       }

       try {
         SFCGAL::algorithm::Buffer3D buffer(*geom_a, radius, segments);
         return buffer.compute(bufferType);
       } catch (const std::exception &e) {
         std::cerr << "buffer3d error: " << e.what() << "\n";
         return std::nullopt;
       }
     }},

    {"centroid", "Construction", "Compute the geometric centroid of a geometry",
     false,
     "The computed centroid relies on 2D areas, even if geometries are "
     "3D-defined.\n",
     "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::centroid(*geom_a);
     }},

    {"chamfer", "Processing",
     "Apply a chamfer or fillet operation to a solid along an edge", true,
     "Parameters:\n"
     "  radius=VALUE: Primary radius/width (default: 0.1)\n"
     "  radius_y=VALUE: Secondary radius for asymmetric chamfer (default: -1, "
     "symmetric)\n"
     "  type=0|1: Operation type (0=flat, 1=round) (default: 0)\n"
     "  segments=INT: Number of segments for round fillet (default: 8)\n"
     "  epsilon=VALUE: Tolerance for edge matching (default: 1e-6)\n\n"
     "Examples:\n"
     "  sfcgalop -a \"SOLID((...))\" -b \"LINESTRING(...)\" chamfer "
     "\"radius=0.5\"\n"
     "  echo \"SOLID(...)\" | sfcgalop -b \"LINESTRING(...)\" chamfer "
     "\"type=1,radius=0.2,segments=16\"",
     "B", "",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_a || !geom_b) {
         std::cerr
             << "chamfer error: requires two geometries (solid and edge)\n";
         return std::nullopt;
       }

       if (geom_b->geometryTypeId() != SFCGAL::TYPE_LINESTRING &&
           geom_b->geometryTypeId() != SFCGAL::TYPE_MULTILINESTRING) {
         std::cerr << "chamfer error: second geometry must be a LineString or "
                      "MultiLineString (edge)\n";
         return std::nullopt;
       }

       auto                              params = parse_params(args);
       SFCGAL::algorithm::ChamferOptions options;

       if (params.count("radius"))
         options.radius = params["radius"];
       if (params.count("radius_y"))
         options.radius_y = params["radius_y"];
       if (params.count("segments"))
         options.segments = static_cast<int>(params["segments"]);
       if (params.count("epsilon"))
         options.epsilon = params["epsilon"];
       if (params.count("type")) {
         int type_int = static_cast<int>(params["type"]);
         options.type = (type_int == 1) ? SFCGAL::algorithm::ChamferType::ROUND
                                        : SFCGAL::algorithm::ChamferType::FLAT;
       }

       try {
         return SFCGAL::algorithm::chamfer(*geom_a, *geom_b, options);
       } catch (const std::exception &e) {
         std::cerr << "chamfer error: " << e.what() << "\n";
         return std::nullopt;
       }
     }},

    {"centroid_3d", "Construction",
     "Compute the geometric centroid of a 3D geometry", false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::centroid3D(*geom_a);
     }},

    {"convex_hull", "Construction", "Compute the 2D convex hull of a geometry",
     false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::convexHull(*geom_a);
     }},

    {"constrained_delaunay_triangulation", "Construction",
     "Compute Constrained Delaunay Triangulation of a geometry", true,
     "No parameters required.\n\nExample:\n  sfcgalop -a "
     "\"POLYGON((0 0,3 0,3 3,0 3,0 0))\" -b \"LINESTRING(0 1.5, 3 1.5)\" "
     "constrained_delaunay_triangulation",
     "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       auto triangulation = SFCGAL::triangulate::triangulate2DZ(*geom_a);
       SFCGAL::triangulate::triangulate2DZ(*geom_b, triangulation);

       return triangulation.getTriangulatedSurface();
     }},

    {"convex_hull_3d", "Construction",
     "Compute the 3D convex hull of a geometry", false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::convexHull3D(*geom_a);
     }},

    {"delaunay_triangulation", "Construction",
     "Compute Delaunay Triangulation of a geometry", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,3 0,3 "
     "3,0 3,0 0))\" delaunay_triangulation",
     "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::triangulate::triangulate2DZ(*geom_a)
           .getTriangulatedSurface();
     }},

    {"envelope", "Construction", "Compute the minimum bounding rectangle",
     false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       SFCGAL::Envelope env = geom_a->envelope();
       return env.toPolygon();
     }},

    {"extrude", "Construction", "Extrude a 2D geometry to create a 3D solid",
     false,
     "Parameters:\n  dx=X: X-axis extrusion distance\n  dy=Y: Y-axis extrusion "
     "distance\n  dz=Z: Z-axis extrusion distance (default: 1.0)\n\nExample:\n "
     " sfcgalop -a \"POLYGON((0 0,1 0,1 1,0 1,0 0))\" extrude "
     "\"dx=0,dy=0,dz=2\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double dx     = params["dx"];
       double dy     = params["dy"];
       double dz     = params.count("dz") ? params["dz"] : 1.0;
       return SFCGAL::algorithm::extrude(*geom_a, dx, dy, dz);
     }},

    {"extrude_straight_skeleton", "Construction",
     "Extrude a polygon using straight skeleton", false,
     "Creates a 3D roof-like surface by extruding a polygon along its straight "
     "skeleton.\n"
     "If building_height is specified, creates a building with vertical walls "
     "and roof.\n\n"
     "Parameters:\n"
     "  height=VALUE: Roof extrusion height (required)\n"
     "  building_height=VALUE: Height of vertical walls (optional, default: "
     "0)\n"
     "    - If > 0, creates a building with walls + roof\n"
     "    - If 0 or omitted, creates only the roof surface\n"
     "  angles=[[a1,a2,...],[b1,b2,...]]: Per-edge angles in degrees "
     "(optional)\n"
     "    - Each inner array corresponds to a ring (exterior, then holes)\n"
     "    - Each value is the angle for one edge segment\n"
     "    - Valid range: 0 < angle < 180 (90 = vertical)\n"
     "  weights=[[w1,w2,...],[v1,v2,...]]: Per-edge weights (tan of angles) "
     "(optional)\n"
     "    - Each inner array corresponds to a ring (exterior, then holes)\n"
     "    - Each value is the weight (tan) for one edge segment\n"
     "    - Cannot be used with angles parameter\n\n"
     "Examples:\n"
     "  # Roof only\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 10,0 10,0 0))\" "
     "extrude_straight_skeleton \"height=5\"\n"
     "  # Building with walls and roof\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 10,0 10,0 0))\" "
     "extrude_straight_skeleton \"height=3,building_height=9\"\n"
     "  # With custom angles\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 10,0 10,0 0))\" "
     "extrude_straight_skeleton \"height=5,angles=[[45,45,45,45]]\"\n"
     "  # With custom weights (gable pattern: vertical on edges 1,3)\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 5,0 5,0 0))\" "
     "extrude_straight_skeleton \"height=3,weights=[[1.0,0,1.0,0]]\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double height = params.count("height") != 0 ? params["height"] : 0.0;
       double building_height = params.count("building_height") != 0
                                    ? params["building_height"]
                                    : 0.0;

       auto angles  = extract_angles_param(args);
       auto weights = extract_weights_param(args);

       // Check for conflict
       if (!angles.empty() && !weights.empty()) {
         throw SFCGAL::Exception(
             "Cannot specify both angles and weights parameters");
       }

       if (angles.empty()) {
         angles = {{}};
       }
       if (weights.empty()) {
         weights = {{}};
       }

       if (building_height > 0.0) {
         // Building with walls + roof
         return SFCGAL::algorithm::extrudeStraightSkeleton(
             *geom_a, building_height, height, weights, angles);
       }
       // Roof only
       return SFCGAL::algorithm::extrudeStraightSkeleton(*geom_a, height,
                                                         weights, angles);
     }},

    {"generate_flat_roof", "Construction",
     "Generate a flat roof from a polygon", false,
     "Creates a flat roof by simple extrusion.\n\n"
     "Parameters:\n"
     "  height=VALUE: Roof height (required)\n\n"
     "Examples:\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 5,0 5,0 0))\" "
     "generate_flat_roof \"height=3\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double height = params.count("height") != 0 ? params["height"] : 3.0;

       SFCGAL::algorithm::RoofParameters roofParameters;
       roofParameters.type       = SFCGAL::algorithm::RoofType::FLAT;
       roofParameters.roofHeight = height;
       return SFCGAL::algorithm::generateRoof(geom_a->as<SFCGAL::Polygon>(),
                                              roofParameters);
     }},

    {"generate_gable_roof", "Construction",
     "Generate a gable roof from a polygon", false,
     "Creates a gable roof by automatically detecting the gable ends (shortest "
     "edges).\n"
     "Gable edges become vertical (90°), other edges use the specified slope "
     "angle.\n"
     "The roof height is determined by the slope angle and the geometry.\n\n"
     "Parameters:\n"
     "  slope_angle=VALUE: Slope angle for non-gable edges in degrees "
     "(default: "
     "45)\n"
     "  height=VALUE: Maximum clipping height (optional, default: 0 = no "
     "clipping)\n\n"
     "Examples:\n"
     "  # Rectangle with auto-detected gables\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 5,0 5,0 0))\" "
     "generate_gable_roof \"slope_angle=30\"\n"
     "  # With max clipping height\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 5,0 5,0 0))\" "
     "generate_gable_roof \"slope_angle=30,height=3\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double height = params.count("height") != 0 ? params["height"] : 0.0;
       double slope_angle =
           params.count("slope_angle") != 0 ? params["slope_angle"] : 45.0;

       SFCGAL::algorithm::RoofParameters roofParameters;
       roofParameters.type       = SFCGAL::algorithm::RoofType::GABLE;
       roofParameters.roofHeight = height;
       roofParameters.slopeAngle = slope_angle;
       return SFCGAL::algorithm::generateRoof(geom_a->as<SFCGAL::Polygon>(),
                                              roofParameters);
     }},

    {"generate_hipped_roof", "Construction",
     "Generate a hipped roof from a polygon", false,
     "Creates a hipped roof using straight skeleton extrusion.\n\n"
     "Parameters:\n"
     "  height=VALUE: Roof height (default: 3.0)\n\n"
     "Examples:\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 5,0 5,0 0))\" "
     "generate_hipped_roof \"height=3\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double height = params.count("height") != 0 ? params["height"] : 3.0;

       SFCGAL::algorithm::RoofParameters roofParameters;
       roofParameters.type       = SFCGAL::algorithm::RoofType::HIPPED;
       roofParameters.roofHeight = height;
       return SFCGAL::algorithm::generateRoof(geom_a->as<SFCGAL::Polygon>(),
                                              roofParameters);
     }},

    {"generate_roof", "Construction",
     "Generate a roof from a polygon using a specified type", false,
     "Unified roof generation dispatching to flat, hipped, gable, or "
     "skillion.\n\n"
     "Parameters:\n"
     "  type=VALUE: Roof type: flat, hipped, gable, skillion (default: "
     "gable)\n"
     "  height=VALUE: Roof height (default: 3.0)\n"
     "  slope_angle=VALUE: Slope angle in degrees for gable/skillion "
     "(default: 30)\n"
     "  primary_edge=VALUE: Index of the sloped edge for skillion "
     "(default: 0)\n\n"
     "Examples:\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 5,0 5,0 0))\" "
     "generate_roof \"type=hipped,height=5\"\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 5,0 5,0 0))\" "
     "generate_roof \"type=gable,height=3,slope_angle=45\"\n"
     "  sfcgalop -a \"POLYGON((0 0,4 0,4 4,0 4,0 0))\" "
     "generate_roof \"type=skillion,height=3,primary_edge=2\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double height = params.count("height") != 0 ? params["height"] : 3.0;
       double slope_angle =
           params.count("slope_angle") != 0 ? params["slope_angle"] : 30.0;
       size_t primary_edge = params.count("primary_edge") != 0
                                 ? static_cast<size_t>(params["primary_edge"])
                                 : 0;

       SFCGAL::algorithm::RoofParameters roofParameters;
       roofParameters.roofHeight       = height;
       roofParameters.slopeAngle       = slope_angle;
       roofParameters.primaryEdgeIndex = primary_edge;

       // Parse roof type from string
       std::string type_str;
       auto        pos = args.find("type=");
       if (pos != std::string::npos) {
         auto start = pos + 5;
         auto end   = args.find_first_of(",\"", start);
         type_str   = args.substr(start, end - start);
       }

       if (type_str == "flat") {
         roofParameters.type = SFCGAL::algorithm::RoofType::FLAT;
       } else if (type_str == "hipped") {
         roofParameters.type = SFCGAL::algorithm::RoofType::HIPPED;
       } else if (type_str == "skillion") {
         roofParameters.type = SFCGAL::algorithm::RoofType::SKILLION;
       } else {
         roofParameters.type = SFCGAL::algorithm::RoofType::GABLE;
       }

       return SFCGAL::algorithm::generateRoof(geom_a->as<SFCGAL::Polygon>(),
                                              roofParameters);
     }},

    {"generate_skillion_roof", "Construction",
     "Generate a skillion (mono-pitch) roof from a polygon", false,
     "Creates a skillion roof with one sloped edge and other edges "
     "vertical.\n\n"
     "Parameters:\n"
     "  height=VALUE: Roof height (default: 3.0)\n"
     "  slope_angle=VALUE: Slope angle in degrees (default: 30)\n"
     "  primary_edge=VALUE: Index of the sloped edge (default: 0)\n\n"
     "Examples:\n"
     "  # Default slope at 30°\n"
     "  sfcgalop -a \"POLYGON((0 0,4 0,4 4,0 4,0 0))\" "
     "generate_skillion_roof \"height=3\"\n"
     "  # Custom edge and angle\n"
     "  sfcgalop -a \"POLYGON((0 0,4 0,4 4,0 4,0 0))\" "
     "generate_skillion_roof \"height=3,primary_edge=1,slope_angle=25\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double height = params.count("height") != 0 ? params["height"] : 3.0;
       double slope_angle =
           params.count("slope_angle") != 0 ? params["slope_angle"] : 30.0;
       size_t primary_edge = params.count("primary_edge") != 0
                                 ? static_cast<size_t>(params["primary_edge"])
                                 : 0;

       SFCGAL::algorithm::RoofParameters roofParameters;
       roofParameters.type             = SFCGAL::algorithm::RoofType::SKILLION;
       roofParameters.roofHeight       = height;
       roofParameters.slopeAngle       = slope_angle;
       roofParameters.primaryEdgeIndex = primary_edge;
       return SFCGAL::algorithm::generateRoof(geom_a->as<SFCGAL::Polygon>(),
                                              roofParameters);
     }},

    {"insert_points_within_tolerance", "Construction",
     "Insert points from geometry B into geometry A within tolerance", true,
     "Parameters:\n  tolerance=VALUE: Maximum distance for point insertion "
     "(default: 1e-9)\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,10 0,10 10,0 "
     "10,0 0))\" -b \"MULTILINESTRING((5 0,5 10),(0 5,10 5))\" "
     "insert_points_within_tolerance \"tolerance=0.01\"",
     "A, B, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       auto   params = parse_params(args);
       double tolerance =
           params.count("tolerance") ? params["tolerance"] : 1e-9;
       return SFCGAL::algorithm::insertPointsWithinTolerance(*geom_a, *geom_b,
                                                             tolerance);
     }},

    {"line_substring", "Construction",
     "Extract substring from linestring by fraction", false,
     "Parameters:\n  start=VALUE: Start fraction (0.0 to 1.0)\n  end=VALUE: "
     "End fraction (default: 1.0)\n\nExample:\n  sfcgalop -a \"LINESTRING(0 "
     "0,10 0,10 10)\" line_substring \"start=0.25,end=0.75\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double start  = params["start"];
       double end    = params.count("end") ? params["end"] : 1.0;
       if (geom_a->geometryTypeId() == SFCGAL::TYPE_LINESTRING) {
         const auto &lineString = geom_a->as<SFCGAL::LineString>();
         return SFCGAL::algorithm::lineSubstring(lineString, start, end);
       }
       return std::nullopt;
     }},

    {"medial_axis", "Construction",
     "Compute the approximate medial axis of a polygon", false,
     "Computes the approximate medial axis for a polygon as a "
     "MultiLineString.\n"
     "The medial axis represents the 'skeleton' of the polygon, where each "
     "point\n"
     "is equidistant from the nearest polygon edges.\n\n"
     "Input A: Polygon geometry\n\n"
     "Parameters:\n"
     "  project_to_edges=BOOL: If true, extend free endpoints to polygon "
     "boundary\n"
     "                         using edge midpoint method (default: false)\n\n"
     "Examples:\n"
     "  sfcgalop -a \"POLYGON((0 0,6 0,6 3,0 3,0 0))\" medial_axis\n"
     "  sfcgalop -a \"POLYGON((0 0,6 0,6 3,0 3,0 0))\" medial_axis "
     "\"project_to_edges=true\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto params = parse_params(args);
       bool projectToEdges =
           parse_boolean_param(params, "project_to_edges", args, false);
       return SFCGAL::algorithm::approximateMedialAxis(*geom_a, projectToEdges);
     }},

    {"minkowski_sum", "Construction", "Compute Minkowski sum of two geometries",
     true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       if (geom_b->geometryTypeId() == SFCGAL::TYPE_POLYGON) {
         return SFCGAL::algorithm::minkowskiSum(*geom_a,
                                                geom_b->as<SFCGAL::Polygon>());
       }
       return std::nullopt;
     }},

    {"minkowski_sum_3d", "Construction",
     "Compute 3D Minkowski sum of two geometries", true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::minkowskiSum3D(*geom_a, *geom_b);
     }},

    {"offset", "Construction", "Create an offset polygon at specified distance",
     false,
     "Parameters:\n  distance=VALUE: Offset distance (default: 1.0)\n  "
     "Positive values create outward offset\n  Negative values create inward "
     "offset\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,3 0,3 3,0 3,0 0))\" "
     "offset \"0.5\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params   = parse_params(args);
       double distance = params.count("distance") ? params["distance"] : 1.0;
       return SFCGAL::algorithm::offset(*geom_a, distance);
     }},

    {"straight_skeleton", "Construction",
     "Compute the straight skeleton of a polygon", false,
     "Parameters:\n  auto_orientation=BOOL: Enable automatic orientation "
     "correction (default: false)\n"
     "                         Accepts: true/false, t/f, 1/0, TRUE/FALSE "
     "(case-insensitive)\n\n"
     "Example:\n  sfcgalop -a \"POLYGON((0 0,4 0,4 "
     "4,0 4,0 0))\" straight_skeleton \"auto_orientation=true\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto params = parse_params(args);
       bool autoOrientation =
           parse_boolean_param(params, "auto_orientation", args, false);
       return SFCGAL::algorithm::straightSkeleton(*geom_a, autoOrientation);
     }},

    {"tessellate", "Construction",
     "Tessellate a geometry into triangular faces", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,3 0,3 "
     "3,0 3,0 0))\" tessellate",
     "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::tessellate(*geom_a);
     }},

    {"triangulate", "Construction",
     "Triangulate a geometry (alias for tessellate)", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,3 0,3 "
     "3,0 3,0 0))\" triangulate",
     "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::tessellate(*geom_a);
     }},

    {"sweep", "Construction",
     "Sweep a 2D profile along a 3D path to create a 3D surface", true,
     "Parameters:\n"
     "  frame_method=0|1|3: Frame computation method (default: 0)\n"
     "    0 = ROTATION_MINIMIZING (RMF) - minimal twist for general paths\n"
     "    1 = FRENET - Frenet-Serret frames\n"
     "    2 = SEGMENT_ALIGNED - segment aligned frames (ideal for "
     "architecture)\n"
     "  start_cap=0|1|2: Start cap style (default: 1)\n"
     "    0 = NONE - no cap\n"
     "    1 = FLAT - flat planar cap\n"
     "  end_cap=0|1: End cap style (default: 1)\n"
     "  closed_path=BOOL: Whether path forms a closed loop (default: false)\n"
     "  anchor_x=VALUE: X coordinate of anchor point in profile space "
     "(default: 0.0)\n"
     "  anchor_y=VALUE: Y coordinate of anchor point in profile space "
     "(default: 0.0)\n\n"
     "The profile must be a 2D geometry (Polygon or LineString).\n"
     "The path must be a 3D LineString with at least 2 points.\n"
     "The anchor point specifies which profile point is positioned on the "
     "path.\n\n"
     "Examples:\n"
     "  sfcgalop -a \"LINESTRING Z(0 0 0,1 0 1,2 0 2)\" -b \"POLYGON((0 0,0.1 "
     "0,0.1 0.1,0 0.1,0 0))\" sweep\n"
     "  sfcgalop -a \"LINESTRING Z(0 0 0,1 1 0,1 1 1,0 1 1,0 0 1,0 0 0)\" -b "
     "\"POLYGON((-0.1 -0.1,0.1 -0.1,0.1 0.1,-0.1 0.1,-0.1 -0.1))\" sweep "
     "\"closed_path=true\"\n"
     "  sfcgalop -a \"LINESTRING Z(0 0 0,10 0 0)\" -b \"POLYGON((0 0,2 0,2 1,0 "
     "1,0 0))\" sweep \"anchor_x=1,anchor_y=0.5\"",
     "A, B, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }

       // Check that geom_a is a LineString
       if (geom_a->geometryTypeId() != SFCGAL::TYPE_LINESTRING) {
         std::cerr
             << "sweep error: first geometry must be a LineString (path)\n";
         return std::nullopt;
       }

       const auto &path = geom_a->as<SFCGAL::LineString>();

       // Parse parameters
       auto params = parse_params(args);

       SFCGAL::algorithm::SweepOptions options;

       // Frame method
       int frame_method = static_cast<int>(
           params.count("frame_method") ? params["frame_method"] : 0);
       switch (frame_method) {
       case 1:
         options.frame_method =
             SFCGAL::algorithm::SweepOptions::FrameMethod::FRENET;
         break;
       case 2:
         options.frame_method =
             SFCGAL::algorithm::SweepOptions::FrameMethod::SEGMENT_ALIGNED;
         break;
       default:
         options.frame_method =
             SFCGAL::algorithm::SweepOptions::FrameMethod::ROTATION_MINIMIZING;
       }

       // Start cap style
       int start_cap = static_cast<int>(
           params.count("start_cap") ? params["start_cap"] : 1);
       switch (start_cap) {
       case 0:
         options.start_cap = SFCGAL::algorithm::SweepOptions::EndCapStyle::NONE;
         break;
       default:
         options.start_cap = SFCGAL::algorithm::SweepOptions::EndCapStyle::FLAT;
       }

       // End cap style
       int end_cap =
           static_cast<int>(params.count("end_cap") ? params["end_cap"] : 1);
       switch (end_cap) {
       case 0:
         options.end_cap = SFCGAL::algorithm::SweepOptions::EndCapStyle::NONE;
         break;
       default:
         options.end_cap = SFCGAL::algorithm::SweepOptions::EndCapStyle::FLAT;
       }

       // Closed path
       options.closed_path =
           parse_boolean_param(params, "closed_path", args, false);

       // Anchor point parameters
       options.anchor_x = params.count("anchor_x") ? params["anchor_x"] : 0.0;
       options.anchor_y = params.count("anchor_y") ? params["anchor_y"] : 0.0;

       try {
         return SFCGAL::algorithm::sweep(path, *geom_b, options);
       } catch (const std::exception &e) {
         std::cerr << "sweep error: " << e.what() << "\n";
         return std::nullopt;
       }
     }}};

} // namespace Operations
