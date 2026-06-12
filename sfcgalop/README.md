#SFCGALOP - SFCGAL Geometry Operations CLI

A command-line interface for performing geometric operations using the SFCGAL library largely inspired by geosop.

## Features

## Building

### Build Instructions

```bash
cmake ... -DSFCGAL_BUILD_CLI=ON
```

## Usage

### Basic Syntax

```bash
sfcgalop [options] -a <WKT/WKB> [-b <WKT/WKB>] [operation] [params]
```

**Note**: If no operation is specified, the input geometry will be displayed in the requested format (useful for format conversion).

### Options

```
   Option                       Description
   -a, --geom-a=ARG             Source for geometry A (WKT, WKB, file, stdin)
   -b, --geom-b=ARG             Source for geometry B (WKT, WKB, file, stdin)
   -f, --format=ARG             Output format: wkt, wkb, txt/ewkt, obj, geojson/json, stl, vtk (default: wkt)
   -p, --precision=N            Decimal precision for output (default: 6)

   Output Control:
   -q, --quiet                  Disable result output
   -t, --time                   Print execution time
   -v, --verbose                Verbose output
   -l, --list                   List all available operations
   --validate                   Validate geometries before operations

   Information:
   -h, --help                   Show this help message
   -V, --version                Show version information
   --help-op=OPERATION          Show help for a specific opera
```

### Examples

#### Calculate area of a polygon

```bash
sfcgalop -a "POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))" area
#Output : 100
```

#### Calculate distance between two points

```bash
sfcgalop -a "POINT(0 0)" -b "POINT(3 4)" distance
#Output : 5
```

#### Compute intersection with validation

```bash
sfcgalop --validate \
  -a "POLYGON((0 0, 4 0, 4 4, 0 4, 0 0))" \
  -b "POLYGON((2 2, 6 2, 6 6, 2 6, 2 2))" \
  intersection -p 0
#Output : POLYGON((2 2, 2 4, 4 4, 4 2, 2 2))
```

#### Read from stdin

```bash
echo "POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))" | \
  sfcgalop -a stdin convexhull
```

#### Convert geometry format

```bash
sfcgalop -a "POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))" -f wkb
#Output : Binary WKB representation
```

#### Convert geometry to OBJ format

```bash
sfcgalop -a "TRIANGLE((0 0 0, 1 0 0, 0 1 0, 0 0 0))" -f obj
```

#### Convert geometry to STL format

```bash
sfcgalop -a "SOLID((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)), ((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)), ((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0)), ((0 0 1,0 1 1,0 1 0,0 0 0,0 0 1)), ((1 0 1,1 1 1,0 1 1,0 0 1,1 0 1)), ((1 0 0,1 0 1,0 0 1,0 0 0,1 0 0))))" -f stl
```

#### Convert geometry to VTK format

```bash
sfcgalop -a "TIN Z (((0 0 0, 0 0 1, 0 1 0, 0 0 0)), ((0 0 0, 0 1 0, 1 0 0, 0 0 0)))" -f vtk
```

#### Display geometry from file

```bash
echo "POINT(0 0)" > /tmp/geometry.wkt
sfcgalop -a /tmp/geometry.wkt
#Output : Geometry displayed in WKT format
```

#### 3D operations

```bash
sfcgalop -a "SOLID((((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),
  ((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),
  ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),
  ((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),
  ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),
  ((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1))))" volume
#Output : 1
```

## Available Operations

Use `sfcgalop --list` to list all operations

```
╔══════════════════════════════════════════════════════════════╗
║                 SFCGAL Available Operations                  ║
╚══════════════════════════════════════════════════════════════╝

╭───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│                                                   Complete Operations Reference                                                   │
├──────────────────────────────────────┬──────────────┬──────────┬──────────────────────────────────────────────────────────────────┤
│ Operation                            │ Input        │ Output   │ Description                                                      │
├──────────────────────────────────────┼──────────────┼──────────┼──────────────────────────────────────────────────────────────────┤
│ ▶ Analysis                           │              │          │                                                                  │
│   normal                             │ A            │ G        │ Compute surface normal vector for polygon/triangle               │
│   orientation                        │ A            │ D        │ Determine polygon ring orientation (clockwise/counter-clockwise) │
│   partition                          │ A, params    │ G        │ Partition polygon into simpler pieces                            │
│   visibility                         │ A, params    │ G        │ Compute visibility polygon from a point in polygon               │
│                                      │              │          │                                                                  │
│ ▶ Collections                        │              │          │                                                                  │
│   collect                            │ A, B         │ G        │ Combine two geometries into a collection                         │
│   collection_extract                 │ A            │ G        │ Extract polygons from geometry collection                        │
│   collection_homogenize              │ A            │ G        │ Convert collection to appropriate multi-type                     │
│   collection_to_multi                │ A            │ G        │ Convert collection to multi-geometry type                        │
│                                      │              │          │                                                                  │
│ ▶ Construction                       │              │          │                                                                  │
│   alpha_shapes                       │ A, params    │ G        │ Compute alpha shapes from point cloud                            │
│   alpha_wrapping_3d                  │ A, params    │ G        │ Create 3D alpha wrapping surface from points                     │
│   boundary                           │ A            │ G        │ Compute the topological boundary of a geometry                   │
│   buffer_3d                          │ A, params    │ G        │ Create a 3D buffer around points and lines                       │
│   centroid                           │ A            │ G        │ Compute the geometric centroid of a geometry                     │
│   centroid_3d                        │ A            │ G        │ Compute the geometric centroid of a 3D geometry                  │
│   convex_hull                        │ A            │ G        │ Compute the 2D convex hull of a geometry                         │
│   constrained_delaunay_triangulation │ A, B         │ G        │ Compute Constrained Delaunay Triangulation of a geometry         │
│   convex_hull_3d                     │ A            │ G        │ Compute the 3D convex hull of a geometry                         │
│   delaunay_triangulation             │ A            │ G        │ Compute Delaunay Triangulation of a geometry                     │
│   envelope                           │ A            │ G        │ Compute the minimum bounding rectangle                           │
│   extrude                            │ A, params    │ G        │ Extrude a 2D geometry to create a 3D solid                       │
│   extrude_straight_skeleton          │ A, params    │ G        │ Extrude a polygon using straight skeleton                        │
│   generate_flat_roof                 │ A, params    │ G        │ Generate a flat roof from a polygon                              │
│   generate_gable_roof                │ A, params    │ G        │ Generate a gable roof from a polygon                             │
│   generate_hipped_roof               │ A, params    │ G        │ Generate a hipped roof from a polygon                            │
│   generate_roof                      │ A, params    │ G        │ Generate a roof from a polygon using a specified type            │
│   generate_skillion_roof             │ A, params    │ G        │ Generate a skillion (mono-pitch) roof from a polygon             │
│   insert_points_within_tolerance     │ A, B, params │ G        │ Insert points from geometry B into geometry A within tolerance   │
│   line_substring                     │ A, params    │ G        │ Extract substring from linestring by fraction                    │
│   medial_axis                        │ A, params    │ G        │ Compute the approximate medial axis of a polygon                 │
│   minkowski_sum                      │ A, B         │ G        │ Compute Minkowski sum of two geometries                          │
│   minkowski_sum_3d                   │ A, B         │ G        │ Compute 3D Minkowski sum of two geometries                       │
│   offset                             │ A, params    │ G        │ Create an offset polygon at specified distance                   │
│   straight_skeleton                  │ A, params    │ G        │ Compute the straight skeleton of a polygon                       │
│   tessellate                         │ A            │ G        │ Tessellate a geometry into triangular faces                      │
│   triangulate                        │ A            │ G        │ Triangulate a geometry (alias for tessellate)                    │
│   sweep                              │ A, B, params │ G        │ Sweep a 2D profile along a 3D path to create a 3D surface        │
│                                      │              │          │                                                                  │
│ ▶ Constructors                       │              │          │                                                                  │
│   make_box                           │ params       │ G        │ Create a 3D box primitive                                        │
│   make_cone                          │ params       │ G        │ Create a 3D cone primitive (supports truncated cones)            │
│   make_cube                          │ params       │ G        │ Create a 3D cube primitive                                       │
│   make_cylinder                      │ params       │ G        │ Create a 3D cylinder primitive                                   │
│   make_sphere                        │ params       │ G        │ Create a 3D sphere primitive using icosahedron subdivision       │
│   make_torus                         │ params       │ G        │ Create a 3D torus primitive                                      │
│                                      │              │          │                                                                  │
│ ▶ Conversions                        │              │          │                                                                  │
│   to_solid                           │ A            │          │ Convert a PolyhedralSurface to a Solid                           │
│                                      │              │          │                                                                  │
│ ▶ Metrics                            │              │          │                                                                  │
│   area                               │ A            │ D        │ Calculate the 2D area of a geometry                              │
│   area_3d                            │ A            │ D        │ Calculate the 3D surface area of a geometry                      │
│   distance                           │ A, B         │ D        │ Calculate the 2D minimum distance between two geometries         │
│   distance_3d                        │ A, B         │ D        │ Calculate the 3D minimum distance between two geometries         │
│   length                             │ A            │ D        │ Calculate the 2D length of linear geometries                     │
│   length_3d                          │ A            │ D        │ Calculate the 3D length of linear geometries                     │
│   volume                             │ A            │ D        │ Calculate the 3D volume of a solid geometry                      │
│                                      │              │          │                                                                  │
│ ▶ Predicates                         │              │          │                                                                  │
│   covers                             │ A, B         │ B        │ Test if geometry A completely covers geometry B                  │
│   intersects                         │ A, B         │ B        │ Test if two geometries intersect in 2D                           │
│   intersects_3d                      │ A, B         │ B        │ Test if two geometries intersect in 3D                           │
│   is_3d                              │ A            │ B        │ Test if geometry has Z coordinates                               │
│   is_closed                          │ A            │ B        │ Test if linear geometry forms a closed ring                      │
│   is_empty                           │ A            │ B        │ Test if geometry contains no points                              │
│   is_measured                        │ A            │ B        │ Test if geometry has measure (M) coordinates                     │
│   is_simple                          │ A            │ B        │ Test if geometry has no self-intersections                       │
│   is_valid                           │ A            │ B        │ Test if geometry is topologically valid                          │
│                                      │              │          │                                                                  │
│ ▶ Processing                         │              │          │                                                                  │
│   chamfer                            │ B            │          │ Apply a chamfer or fillet operation to a solid along an edge     │
│                                      │              │          │                                                                  │
│ ▶ Boolean Operations                 │              │          │                                                                  │
│   difference                         │ A, B         │ G        │ Compute geometry A minus geometry B                              │
│   difference_3d                      │ A, B         │ G        │ Compute 3D geometry A minus geometry B                           │
│   intersection                       │ A, B         │ G        │ Compute the geometric intersection of two geometries             │
│   intersection_3d                    │ A, B         │ G        │ Compute the 3D geometric intersection of two geometries          │
│   union                              │ A, B         │ G        │ Compute the geometric union of two geometries                    │
│   union_3d                           │ A, B         │ G        │ Compute the 3D geometric union of two geometries                 │
│                                      │              │          │                                                                  │
│ ▶ Transformations                    │              │          │                                                                  │
│   force_2d                           │ A            │ G        │ Remove Z coordinates to create 2D geometry                       │
│   force_3d                           │ A, params    │ G        │ Add Z coordinates to create 3D geometry                          │
│   force_measured                     │ A, params    │ G        │ Add measure coordinates to geometry                              │
│   force_lhr                          │ A            │ G        │ Force a Left Handed Rule on the given Geometry                   │
│   force_rhr                          │ A            │ G        │ Force a Right Handed Rule on the given Geometry                  │
│   rotate                             │ A, params    │ G        │ Rotate geometry around specified axis                            │
│   scale                              │ A, params    │ G        │ Scale geometry by specified factors                              │
│   simplify                           │ A, params    │ G        │ Simplify geometry by removing vertices within tolerance          │
│   split_3d                           │ A, params    │ G        │ Split geometry with a plane                                      │
│   surface_simplification             │ A, params    │ G        │ Simplify a 3D surface mesh using edge collapse                   │
│   translate                          │ A, params    │ G        │ Translate geometry by specified offset                           │
╰──────────────────────────────────────┴──────────────┴──────────┴──────────────────────────────────────────────────────────────────╯

ℹ Use --help-op=<operation> for detailed help on a specific operation
```

### explenation

`sfcgalop --help-op=extrude`

```
Operation: extrude
Category: Construction
Description: Extrude a 2D geometry to create a 3D solid

Parameters:
  dx=X: X-axis extrusion distance
  dy=Y: Y-axis extrusion distance
  dz=Z: Z-axis extrusion distance (default: 1.0)

Example:
  sfcgalop -a "POLYGON((0 0,1 0,1 1,0 1,0 0))" extrude "dx=0,dy=0,dz=2"
```

Input:
  - A: Geometry A
  - B: Geometry B
  - params: parameters required by the algorithm

Output:
  - G: Geometry
  - B: boolean (true/false)
  - D: double (number)

## Error Handling

The tool provides detailed error messages with context:

```bash
#Invalid geometry
sfcgalop --validate -a "POLYGON((0 0, 0 10, 10 10, 10 0, 1 0))" area
#Output:
# ⚠ Geometry A validation issues:
#Validation Error : ring 0 is not closed
#Details:
#Geometry type : Polygon
#WKT : POLYGON((0 0, 0 10, 10 10, 10 0, 1 0))
```

## Naming Conventions and Compatibility

SFCGAL supports multiple naming conventions for better compatibility and flexibility through a normalization approach:

### Normalization Approach
SFCGAL uses a normalization system that allows multiple naming conventions for the same operation by converting names to lowercase and removing underscores for matching:

- `alpha_wrapping_3d`, `alphawrapping3d`, and `alphaWrapping3D` all refer to the same operation
- `line_substring`, `linesubstring`, and `LineSubstring` all refer to the same operation
- `convex_hull_3d`, `convexhull3d`, and `ConvexHull3D` all refer to the same operation

### Canonical Names
The preferred canonical names follow the underscore convention (e.g., `minkowski_sum`, `alpha_shapes`, `area_3d`).
The `sfcgalop --list` command displays operations using these canonical names.

### Case Insensitivity
All operation names are case-insensitive, so `Area`, `AREA`, and `area` all work the same.

Example usage:
```bash
# Using different naming conventions for the same operation
sfcgalop -a "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))" area_3d
sfcgalop -a "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))" area3d
sfcgalop -a "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))" AREA3D
# All of these work the same way
```

## License

This project is part of SFCGAL and follows the same license terms.

## Authors

SFCGAL Contributors
