// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations_transformations.hpp"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/algorithm/force2D.h"
#include "SFCGAL/algorithm/force3D.h"
#include "SFCGAL/algorithm/forceMeasured.h"
#if SFCGAL_CGAL_VERSION_MAJOR >= 6
  #include "SFCGAL/algorithm/polygonRepair.h"
#endif
#include "SFCGAL/algorithm/rotate.h"
#include "SFCGAL/algorithm/scale.h"
#include "SFCGAL/algorithm/simplification.h"
#include "SFCGAL/algorithm/split3D.h"
#include "SFCGAL/algorithm/surfaceSimplification.h"
#include "SFCGAL/algorithm/translate.h"
#include "SFCGAL/detail/transform/ForceOrderPoints.h"

namespace Operations {

auto
parseStopPredicate(const std::map<std::string, double> &params)
    -> std::optional<SFCGAL::algorithm::SimplificationStopPredicate>
{
  if (params.count("count") != 0) {
    double count_value = params.at("count");
    if (count_value <= 0 || count_value != std::floor(count_value)) {
      return std::nullopt;
    }

    return SFCGAL::algorithm::SimplificationStopPredicate::edgeCount(
        static_cast<size_t>(count_value));
  }

  double ratio = params.count("ratio") != 0 ? params.at("ratio") : 0.5;
  if (ratio <= 0.0 || ratio >= 1.0) {
    return std::nullopt;
  }

  return SFCGAL::algorithm::SimplificationStopPredicate::edgeCountRatio(ratio);
}

auto
parseSimplificationStrategy(const std::string                   &args,
                            const std::map<std::string, double> &params)
    -> std::optional<SFCGAL::algorithm::SimplificationStrategy>
{
  using Strategy = SFCGAL::algorithm::SimplificationStrategy;

  if (params.count("strategy") == 0) {
    return Strategy::EDGE_LENGTH;
  }

  auto pos = args.find("strategy=");
  if (pos == std::string::npos) {
    return Strategy::EDGE_LENGTH;
  }

  auto start = pos + 9;
  auto end   = args.find(',', start);
  if (end == std::string::npos) {
    end = args.length();
  }

  std::string value = trim(args.substr(start, end - start));

  if (value == "edge_length") {
    return Strategy::EDGE_LENGTH;
  }

#ifdef SFCGAL_WITH_EIGEN
  if (value == "garland_heckbert") {
    return Strategy::GARLAND_HECKBERT;
  }
  if (value == "lindstrom_turk") {
    return Strategy::LINDSTROM_TURK;
  }
#endif

  return std::nullopt;
}

const std::vector<Operation> operations_transformations = {
    {"force_2d", "Transformations",
     "Remove Z coordinates to create 2D geometry", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"POINT Z(1 2 3)\" "
     "force_2d",
     "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto result = geom_a->clone();
       SFCGAL::algorithm::force2D(*result);
       return result;
     }},

    {"force_3d", "Transformations", "Add Z coordinates to create 3D geometry",
     false,
     "Parameters:\n  z=VALUE: Z coordinate value to assign (default: "
     "0.0)\n\nExample:\n  sfcgalop -a \"POINT(1 2)\" force_3d \"z=5.0\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double z      = params.count("z") ? params["z"] : 1.0;
       auto   result = geom_a->clone();
       SFCGAL::algorithm::force3D(*result, z);
       return result;
     }},

    {"force_measured", "Transformations", "Add measure coordinates to geometry",
     false,
     "Parameters:\n  m=VALUE: Measure coordinate value (default: "
     "0.0)\n\nExample:\n  sfcgalop -a \"POINT(1 2)\" force_measured \"m=10.0\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double m      = params.count("m") ? params["m"] : 1.0;
       auto   result = geom_a->clone();
       SFCGAL::algorithm::forceMeasured(*result, m);
       return result;
     }},

    {"forceLHR", "Transformations",
     "Force a Left Handed Rule on the given Geometry", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"MULTIPOLYGON (((9 "
     "9, 9 1, 1 1, 2 4, 7 7, 9 9)))\" "
     "forceLHR",
     "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto                                result = geom_a->clone();
       SFCGAL::transform::ForceOrderPoints force(/* ccw */ true);
       result->accept(force);
       return result;
     }},

    {"forceRHR", "Transformations",
     "Force a Right Handed Rule on the given Geometry", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"MULTIPOLYGON (((9 "
     "9,7 7,2 4,1 1,9 1,9 9)))\" "
     "forceRHR",
     "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto                                result = geom_a->clone();
       SFCGAL::transform::ForceOrderPoints force(/* ccw */ false);
       result->accept(force);
       return result;
     }},

#if SFCGAL_CGAL_VERSION_MAJOR >= 6
    {"polygon_repair", "Transformations", "Repair invalid polygons with rules",
     false,
     "Parameters:\n  method=0|1|2|3 (default: 0)\n\nMethods:\n  0 = "
     "Even-odd\n  1 = Non-zero winding\n  2 = Union of all polygons\n  3 = "
     "Intersection of all polygons\n\n"
     "Example:\n  sfcgalop -a \"POLYGON((0 0, 2 2, 2 0, 0 2, 0 0)) "
     "polygon_repair \"method=1\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto params = parse_params(args);
       int  method =
           static_cast<int>(params.count("method") ? params["method"] : 0);

       SFCGAL::algorithm::PolygonRepairRule rule;
       switch (method) {
       case 1:
         rule = SFCGAL::algorithm::PolygonRepairRule::NON_ZERO_RULE;
         break;
       case 2:
         rule = SFCGAL::algorithm::PolygonRepairRule::UNION_RULE;
         break;
       case 3:
         rule = SFCGAL::algorithm::PolygonRepairRule::INTERSECTION_RULE;
         break;
       default:
         rule = SFCGAL::algorithm::PolygonRepairRule::EVEN_ODD_RULE;
       }

       return SFCGAL::algorithm::polygonRepair(*geom_a, rule);
     }},
#endif // SFCGAL_CGAL_VERSION_MAJOR >= 6

    {"rotate", "Transformations", "Rotate geometry around specified axis",
     false,
     "Parameters:\n  angle=DEGREES: Rotation angle in degrees\n  axis=x|y|z: "
     "Rotation axis (default: z)\n\nExamples:\n  sfcgalop -a \"POINT(1 0)\" "
     "rotate \"angle=90\"           # Rotate 90° around Z-axis\n  sfcgalop -a "
     "\"POINT(1 0 0)\" rotate \"angle=90,axis=y\"  # Rotate 90° around Y-axis",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto        params    = parse_params(args);
       double      angle_deg = params.count("angle") ? params["angle"] : 0.0;
       std::string axis      = "z"; // default to Z-axis

       // Check if axis parameter is provided (axis won't be in parse_params
       // because it's not a double) Parse manually for non-numeric parameters
       if (args.find("axis=") != std::string::npos) {
         size_t axis_pos = args.find("axis=");
         size_t start    = axis_pos + 5; // length of "axis="
         size_t end      = args.find(',', start);
         if (end == std::string::npos) {
           end = args.length();
         }
         axis = args.substr(start, end - start);
       }

       // Convert degrees to radians
       double angle_rad = (angle_deg * M_PI) / 180.0;

       auto result = geom_a->clone();

       try {
         if (axis == "x" || axis == "X") {
           SFCGAL::algorithm::rotateX(*result, angle_rad);
         } else if (axis == "y" || axis == "Y") {
           SFCGAL::algorithm::rotateY(*result, angle_rad);
         } else {
           // Default to Z-axis rotation
           SFCGAL::algorithm::rotateZ(*result, angle_rad);
         }
         return result;
       } catch (const std::exception &e) {
         std::cerr << "Rotation failed: " << e.what() << "\n";
         return std::nullopt;
       }
     }},

    {"scale", "Transformations", "Scale geometry by specified factors", false,
     "Parameters:\n  s=VALUE: Uniform scaling factor\n  OR\n  sx=X: X-axis "
     "scaling factor\n  sy=Y: Y-axis scaling factor\n  sz=Z: Z-axis scaling "
     "factor\n\nExamples:\n  sfcgalop -a \"POINT(1 1)\" scale \"s=2\"\n  "
     "sfcgalop -a \"POINT(1 1)\" scale \"sx=2,sy=1,sz=1\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double scale  = params.count("s") ? params["s"] : 1.0;
       double sx     = params.count("sx") ? params["sx"] : scale;
       double sy     = params.count("sy") ? params["sy"] : scale;
       double sz     = params.count("sz") ? params["sz"] : scale;
       auto   result = geom_a->clone();
       SFCGAL::algorithm::scale(*result, sx, sy, sz);
       return result;
     }},

    {"simplify", "Transformations",
     "Simplify geometry by removing vertices within tolerance", false,
     "Parameters:\n  tolerance=VALUE: Distance tolerance for vertex removal "
     "(default: 0.01)\n\nExample:\n  sfcgalop -a \"LINESTRING(0 0,0.01 0.01,1 "
     "1,2 2)\" simplify \"tolerance=0.1\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params    = parse_params(args);
       double tolerance = params.count("tolerance") ? params["tolerance"] : 1.0;
       bool   preserveTopology = true;
       return SFCGAL::algorithm::simplify(*geom_a, tolerance, preserveTopology);
     }},

    {"split_3d", "Transformations", "Split geometry with a plane", false,
     "Split the given geometry with a plane defined by a point and a normal "
     "vector\n\n"
     "Input A: A PolyhedralSurface, aSolid or a TIN geometry\n\n"
     "Parameters:\n"
     "  ptx=X: X-coordinate of a point belonging to the splitting plane\n"
     "  pty=Y: Y-coordinate of a point belonging to the splitting plane\n"
     "  ptz=Z: Z-coordinate of a point belonging to the splitting plane\n"
     "  nx=X: X-coordinate of the normal vector of the splitting plane\n"
     "  ny=Y: Y-coordinate of the normal vector of the splitting plane\n"
     "  nz=Z: Z-coordinate of the normal vector of the splitting plane\n"
     "  close_geometries=BOOL: If true, ensures resulting geometries are "
     "closed\n\n"
     "Example:\n "
     "  sfcgalop -a \"SOLID(...)\" split "
     "\"ptx=1,pty=0,ptz=0,nx=1,ny=1,nz=0,close_geometries=true\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double ptx    = params["ptx"];
       double pty    = params["pty"];
       double ptz    = params["ptz"];
       double nx     = params["nx"];
       double ny     = params["ny"];
       double nz     = params["nz"];
       bool   closeGeometries =
           parse_boolean_param(params, "close_geometries", args, false);
       const SFCGAL::Point            plane_pt(ptx, pty, ptz);
       const SFCGAL::Kernel::Vector_3 plane_normal(nx, ny, nz);
       auto result = SFCGAL::algorithm::split3D(*geom_a, plane_pt, plane_normal,
                                                closeGeometries);
       return result;
     }},

    {"surfacesimplification", "Transformations",
     "Simplify a 3D surface mesh using edge collapse", false,
     "Parameters:\n"
     "  ratio=VALUE: Edge count ratio to keep (0.0 to 1.0, default: 0.5)\n"
     "  count=VALUE: Target edge count (alternative to ratio)\n"
     "  strategy=VALUE: Simplification strategy (edge_length, "
     "garland_heckbert, lindstrom_turk, default: edge_length)\n\n"
     "Strategies:\n"
     "  edge_length: Uses edge length cost with midpoint placement (default)\n"
#ifdef SFCGAL_WITH_EIGEN
     "  garland_heckbert: Uses quadric error metrics (requires Eigen)\n"
     "  lindstrom_turk: Uses Lindstrom-Turk cost/placement (requires Eigen)\n"
#endif
     "\nExamples:\n"
     "  sfcgalop -a mesh.obj surfacesimplification \"ratio=0.5\"\n"
     "  sfcgalop -a mesh.obj surfacesimplification \"count=1000\"\n"
     "  sfcgalop -a mesh.obj surfacesimplification "
     "\"ratio=0.3,strategy=edge_length\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       if (!geom_a) {
         return std::nullopt;
       }

       // Check if geometry is supported
       auto geom_type = geom_a->geometryTypeId();
       if (geom_type != SFCGAL::TYPE_TRIANGULATEDSURFACE &&
           geom_type != SFCGAL::TYPE_POLYHEDRALSURFACE &&
           geom_type != SFCGAL::TYPE_SOLID &&
           geom_type != SFCGAL::TYPE_MULTISOLID) {
         return std::nullopt;
       }

       auto params = parse_params(args);

       auto stopPredicate = parseStopPredicate(params);
       if (!stopPredicate) {
         return std::nullopt;
       }

       auto strategy = parseSimplificationStrategy(args, params);
       if (!strategy) {
         return std::nullopt;
       }

       try {
         return SFCGAL::algorithm::surfaceSimplification(
             *geom_a, *stopPredicate, *strategy);
       } catch (const std::invalid_argument &e) {
         // Invalid parameters (ratio out of range, etc.)
         return std::nullopt;
       } catch (const std::bad_alloc &e) {
         // Memory allocation failure during simplification
         return std::nullopt;
       } catch (const std::runtime_error &e) {
         // CGAL or geometric algorithm errors
         return std::nullopt;
       } catch (const std::exception &e) {
         // Other standard exceptions
         return std::nullopt;
       }
     }},

    {"translate", "Transformations", "Translate geometry by specified offset",
     false,
     "Parameters:\n  dx=X: X-axis translation\n  dy=Y: Y-axis translation\n  "
     "dz=Z: Z-axis translation\n\nExample:\n  sfcgalop -a \"POINT(0 0)\" "
     "translate \"dx=5,dy=3,dz=1\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double dx     = params["dx"];
       double dy     = params["dy"];
       double dz     = params["dz"];
       auto   result = geom_a->clone();
       SFCGAL::algorithm::translate(*result, dx, dy, dz);
       return result;
     }}};

} // namespace Operations
