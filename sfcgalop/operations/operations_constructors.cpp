// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations_constructors.hpp"
#include "../constructors.hpp"

namespace Operations {

const std::vector<Operation> operations_constructors = {
    {"make_box", "Constructors", "Create a 3D box primitive", false,
     "Parameters:\n  x_extent=VALUE: Length in X direction (default: 1.0)\n  "
     "y_extent=VALUE: Length in Y direction (default: 1.0)\n  z_extent=VALUE: "
     "Length in Z direction (default: 1.0)\n\nExample:\n  sfcgalop make_box "
     "\"x_extent=2,y_extent=3,z_extent=1\"",
     "params", "G",
     [](const std::string &args, const SFCGAL::Geometry *,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params   = parse_params(args);
       double x_extent = params.count("x_extent") ? params["x_extent"] : 1.0;
       double y_extent = params.count("y_extent") ? params["y_extent"] : 1.0;
       double z_extent = params.count("z_extent") ? params["z_extent"] : 1.0;
       return Constructors::make_box(x_extent, y_extent, z_extent);
     }},

    {"make_cone", "Constructors",
     "Create a 3D cone primitive (supports truncated cones)", false,
     "Parameters:\n  base_x=VALUE: X coordinate of base center (default: "
     "0.0)\n  base_y=VALUE: Y coordinate of base center (default: 0.0)\n  "
     "base_z=VALUE: Z coordinate of base center (default: 0.0)\n  "
     "axis_x=VALUE: X component of cone axis (default: 0.0)\n  axis_y=VALUE: Y "
     "component of cone axis (default: 0.0)\n  axis_z=VALUE: Z component of "
     "cone axis (default: 1.0)\n  bottom_radius=VALUE: Cone bottom radius "
     "(default: "
     "1.0)\n  top_radius=VALUE: Cone top radius - 0.0 for regular cone "
     "(default: 0.0)\n  "
     "height=VALUE: Cone height (default: 1.0)\n  num_radial=N: Number "
     "of radial divisions (default: 32)\n\nExamples:\n  sfcgalop make_cone "
     "\"bottom_radius=2,height=4,num_radial=24\"\n  sfcgalop make_cone "
     "\"bottom_radius=3,top_radius=1,height=5\" # Truncated cone",
     "params", "G",
     [](const std::string &args, const SFCGAL::Geometry *,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double base_x = params.count("base_x") ? params["base_x"] : 0.0;
       double base_y = params.count("base_y") ? params["base_y"] : 0.0;
       double base_z = params.count("base_z") ? params["base_z"] : 0.0;
       double axis_x = params.count("axis_x") ? params["axis_x"] : 0.0;
       double axis_y = params.count("axis_y") ? params["axis_y"] : 0.0;
       double axis_z = params.count("axis_z") ? params["axis_z"] : 1.0;
       double bottom_radius =
           params.count("bottom_radius") ? params["bottom_radius"] : 1.0;
       double top_radius =
           params.count("top_radius") ? params["top_radius"] : 0.0;
       double height     = params.count("height") ? params["height"] : 1.0;
       auto   num_radial = static_cast<unsigned int>(
           params.count("num_radial") ? params["num_radial"] : 32);
       return Constructors::make_cone(base_x, base_y, base_z, axis_x, axis_y,
                                      axis_z, bottom_radius, top_radius, height,
                                      num_radial);
     }},

    {"make_cube", "Constructors", "Create a 3D cube primitive", false,
     "Parameters:\n  size=VALUE: Edge length of the cube (default: "
     "1.0)\n\nExample:\n  sfcgalop make_cube \"size=2.0\"",
     "params", "G",
     [](const std::string &args, const SFCGAL::Geometry *,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double size   = params.count("size") ? params["size"] : 1.0;
       return Constructors::make_cube(size);
     }},

    {"make_cylinder", "Constructors", "Create a 3D cylinder primitive", false,
     "Parameters:\n  base_x=VALUE: X coordinate of base center (default: "
     "0.0)\n  base_y=VALUE: Y coordinate of base center (default: 0.0)\n  "
     "base_z=VALUE: Z coordinate of base center (default: 0.0)\n  "
     "axis_x=VALUE: X component of cylinder axis (default: 0.0)\n  "
     "axis_y=VALUE: Y component of cylinder axis (default: 0.0)\n  "
     "axis_z=VALUE: Z component of cylinder axis (default: 1.0)\n  "
     "radius=VALUE: Cylinder radius (default: 1.0)\n  height=VALUE: Cylinder "
     "height (default: 1.0)\n  num_radial=N: Number of radial divisions "
     "(default: 32)\n\nExample:\n  sfcgalop make_cylinder "
     "\"radius=1.5,height=3,num_radial=16\"",
     "params", "G",
     [](const std::string &args, const SFCGAL::Geometry *,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params     = parse_params(args);
       double base_x     = params.count("base_x") ? params["base_x"] : 0.0;
       double base_y     = params.count("base_y") ? params["base_y"] : 0.0;
       double base_z     = params.count("base_z") ? params["base_z"] : 0.0;
       double axis_x     = params.count("axis_x") ? params["axis_x"] : 0.0;
       double axis_y     = params.count("axis_y") ? params["axis_y"] : 0.0;
       double axis_z     = params.count("axis_z") ? params["axis_z"] : 1.0;
       double radius     = params.count("radius") ? params["radius"] : 1.0;
       double height     = params.count("height") ? params["height"] : 1.0;
       auto   num_radial = static_cast<unsigned int>(
           params.count("num_radial") ? params["num_radial"] : 32);
       return Constructors::make_cylinder(base_x, base_y, base_z, axis_x,
                                          axis_y, axis_z, radius, height,
                                          num_radial);
     }},

    {"make_sphere", "Constructors",
     "Create a 3D sphere primitive using icosahedron subdivision", false,
     "Parameters:\n  x=X_COORD: X coordinate of center (default: 0.0)\n  "
     "y=Y_COORD: Y coordinate of center (default: 0.0)\n  z=Z_COORD: Z "
     "coordinate of center (default: 0.0)\n  radius=VALUE: Sphere radius "
     "(default: 1.0)\n  num_subdivisions=N: Number of icosahedron subdivisions "
     "(default: 2)\n\n"
     "Example:\n  sfcgalop make_sphere "
     "\"x=0,y=0,z=0,radius=2.5,num_subdivisions=3\"",
     "params", "G",
     [](const std::string &args, const SFCGAL::Geometry *,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double x      = params.count("x") ? params["x"] : 0.0;
       double y      = params.count("y") ? params["y"] : 0.0;
       double z      = params.count("z") ? params["z"] : 0.0;
       double radius = params.count("radius") ? params["radius"] : 1.0;

       unsigned int num_subdivisions = 2; // default
       if (params.count("num_subdivisions")) {
         num_subdivisions =
             static_cast<unsigned int>(params["num_subdivisions"]);
       }

       return Constructors::make_sphere(x, y, z, radius, num_subdivisions);
     }},

    {"make_torus", "Constructors", "Create a 3D torus primitive", false,
     "Parameters:\n  center_x=VALUE: X coordinate of center (default: 0.0)\n  "
     "center_y=VALUE: Y coordinate of center (default: 0.0)\n  center_z=VALUE: "
     "Z coordinate of center (default: 0.0)\n  axis_x=VALUE: X component of "
     "torus axis (default: 0.0)\n  axis_y=VALUE: Y component of torus axis "
     "(default: 0.0)\n  axis_z=VALUE: Z component of torus axis (default: "
     "1.0)\n  major_radius=VALUE: Major radius of torus (default: 2.0)\n  "
     "minor_radius=VALUE: Minor radius of torus (default: 0.5)\n  num_major=N: "
     "Number of major divisions (default: 32)\n  num_minor=N: Number of minor "
     "divisions (default: 16)\n\nExample:\n  sfcgalop make_torus "
     "\"major_radius=3,minor_radius=1,num_major=24,num_minor=12\"",
     "params", "G",
     [](const std::string &args, const SFCGAL::Geometry *,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params   = parse_params(args);
       double center_x = params.count("center_x") ? params["center_x"] : 0.0;
       double center_y = params.count("center_y") ? params["center_y"] : 0.0;
       double center_z = params.count("center_z") ? params["center_z"] : 0.0;
       double axis_x   = params.count("axis_x") ? params["axis_x"] : 0.0;
       double axis_y   = params.count("axis_y") ? params["axis_y"] : 0.0;
       double axis_z   = params.count("axis_z") ? params["axis_z"] : 1.0;
       double major_radius =
           params.count("major_radius") ? params["major_radius"] : 2.0;
       double minor_radius =
           params.count("minor_radius") ? params["minor_radius"] : 0.5;
       auto num_major = static_cast<unsigned int>(
           params.count("num_major") ? params["num_major"] : 32);
       auto num_minor = static_cast<unsigned int>(
           params.count("num_minor") ? params["num_minor"] : 16);
       return Constructors::make_torus(center_x, center_y, center_z, axis_x,
                                       axis_y, axis_z, major_radius,
                                       minor_radius, num_major, num_minor);
     }}};

} // namespace Operations
