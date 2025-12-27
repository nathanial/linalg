/-
  Tests for Polygon2D operations.
-/

import Linalg
import Crucible

namespace LinalgTests.Polygon2DTests

open Crucible
open Linalg

testSuite "Polygon2D Basic"

test "triangle has 3 vertices" := do
  let p := Polygon2D.triangle Vec2.zero (Vec2.mk 1.0 0.0) (Vec2.mk 0.5 1.0)
  ensure (p.vertexCount == 3) "should have 3 vertices"

test "rectangle has 4 vertices" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 2.0 3.0)
  ensure (p.vertexCount == 4) "should have 4 vertices"

test "regular hexagon has 6 vertices" := do
  let p := Polygon2D.regular 6 1.0
  ensure (p.vertexCount == 6) "should have 6 vertices"

test "vertex access with wraparound" := do
  let p := Polygon2D.triangle Vec2.zero (Vec2.mk 1.0 0.0) (Vec2.mk 0.5 1.0)
  let v0 := p.vertex 0
  let v3 := p.vertex 3  -- Should wrap to vertex 0
  ensure (floatNear v0.x v3.x 0.0001) "wrapped vertex should match"

test "isValid requires 3+ vertices" := do
  let validPoly := Polygon2D.triangle Vec2.zero (Vec2.mk 1.0 0.0) (Vec2.mk 0.5 1.0)
  let invalidPoly := Polygon2D.mk #[Vec2.zero, Vec2.mk 1.0 0.0]
  ensure validPoly.isValid "triangle should be valid"
  ensure (!invalidPoly.isValid) "2-vertex polygon should be invalid"

testSuite "Polygon2D Area"

test "unit square has area 1" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 1.0 1.0)
  ensure (floatNear p.area 1.0 0.0001) "area should be 1"

test "2x3 rectangle has area 6" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 2.0 3.0)
  ensure (floatNear p.area 6.0 0.0001) "area should be 6"

test "triangle area formula" := do
  -- Right triangle with base 4 and height 3
  let p := Polygon2D.triangle Vec2.zero (Vec2.mk 4.0 0.0) (Vec2.mk 0.0 3.0)
  ensure (floatNear p.area 6.0 0.0001) "area should be 6 (0.5 * 4 * 3)"

testSuite "Polygon2D Winding Order"

test "counter-clockwise rectangle" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 1.0 1.0)
  ensure p.isCounterClockwise "rectangle should be CCW"

test "reverse makes clockwise" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 1.0 1.0)
  let rev := p.reverse
  ensure rev.isClockwise "reversed should be CW"

test "makeCounterClockwise normalizes" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 1.0 1.0)
  let rev := p.reverse
  let normalized := rev.makeCounterClockwise
  ensure normalized.isCounterClockwise "should be CCW after normalization"

test "windingOrder enum values" := do
  let ccw := Polygon2D.rectangle Vec2.zero (Vec2.mk 1.0 1.0)
  let cw := ccw.reverse
  ensure (ccw.windingOrder == WindingOrder.counterClockwise) "should be CCW"
  ensure (cw.windingOrder == WindingOrder.clockwise) "should be CW"

testSuite "Polygon2D Centroid"

test "unit square centroid is at center" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 2.0 2.0)
  let c := p.centroid
  ensure (floatNear c.x 1.0 0.0001) "centroid.x should be 1"
  ensure (floatNear c.y 1.0 0.0001) "centroid.y should be 1"

test "symmetric triangle centroid" := do
  let p := Polygon2D.triangle (Vec2.mk 0.0 0.0) (Vec2.mk 3.0 0.0) (Vec2.mk 1.5 3.0)
  let c := p.centroid
  ensure (floatNear c.x 1.5 0.0001) "centroid.x should be 1.5"
  ensure (floatNear c.y 1.0 0.0001) "centroid.y should be 1"

testSuite "Polygon2D Point Tests"

test "point inside rectangle" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 10.0 10.0)
  ensure (p.containsPoint (Vec2.mk 5.0 5.0)) "center should be inside"

test "point outside rectangle" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 10.0 10.0)
  ensure (!p.containsPoint (Vec2.mk 15.0 5.0)) "outside point should not be inside"

test "point on boundary" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 10.0 10.0)
  ensure (p.isOnBoundary (Vec2.mk 5.0 0.0) 0.0001) "point on edge should be on boundary"

test "containsPointInclusive includes boundary" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 10.0 10.0)
  ensure (p.containsPointInclusive (Vec2.mk 0.0 5.0) 0.0001) "boundary point should be included"

test "point inside triangle" := do
  let p := Polygon2D.triangle Vec2.zero (Vec2.mk 10.0 0.0) (Vec2.mk 5.0 10.0)
  ensure (p.containsPoint (Vec2.mk 5.0 3.0)) "point should be inside triangle"

test "point outside triangle" := do
  let p := Polygon2D.triangle Vec2.zero (Vec2.mk 10.0 0.0) (Vec2.mk 5.0 10.0)
  ensure (!p.containsPoint (Vec2.mk 0.0 5.0)) "point should be outside triangle"

testSuite "Polygon2D Bounding Box"

test "rectangle bounding box" := do
  let p := Polygon2D.rectangle (Vec2.mk 1.0 2.0) (Vec2.mk 5.0 7.0)
  let (minV, maxV) := p.boundingBox
  ensure (floatNear minV.x 1.0 0.0001) "min.x should be 1"
  ensure (floatNear minV.y 2.0 0.0001) "min.y should be 2"
  ensure (floatNear maxV.x 5.0 0.0001) "max.x should be 5"
  ensure (floatNear maxV.y 7.0 0.0001) "max.y should be 7"

test "perimeter of unit square" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 1.0 1.0)
  ensure (floatNear p.perimeter 4.0 0.0001) "perimeter should be 4"

testSuite "Polygon2D Convexity"

test "rectangle is convex" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 5.0 5.0)
  ensure p.isConvex "rectangle should be convex"

test "regular hexagon is convex" := do
  let p := Polygon2D.regular 6 1.0
  ensure p.isConvex "regular hexagon should be convex"

test "concave L-shape is not convex" := do
  -- L-shaped polygon (concave)
  let p := Polygon2D.mk #[
    Vec2.mk 0.0 0.0,
    Vec2.mk 2.0 0.0,
    Vec2.mk 2.0 1.0,
    Vec2.mk 1.0 1.0,
    Vec2.mk 1.0 2.0,
    Vec2.mk 0.0 2.0
  ]
  ensure (!p.isConvex) "L-shape should not be convex"

testSuite "Polygon2D Transformations"

test "translate moves polygon" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 1.0 1.0)
  let translated := p.translate (Vec2.mk 5.0 5.0)
  let (minV, _) := translated.boundingBox
  ensure (floatNear minV.x 5.0 0.0001) "min.x should be 5"
  ensure (floatNear minV.y 5.0 0.0001) "min.y should be 5"

test "scale from origin" := do
  let p := Polygon2D.rectangle (Vec2.mk 1.0 1.0) (Vec2.mk 2.0 2.0)
  let scaled := p.scale 2.0
  let (minV, maxV) := scaled.boundingBox
  ensure (floatNear minV.x 2.0 0.0001) "scaled min.x should be 2"
  ensure (floatNear maxV.x 4.0 0.0001) "scaled max.x should be 4"

test "rotate 90 degrees" := do
  let p := Polygon2D.mk #[Vec2.mk 1.0 0.0, Vec2.mk 2.0 0.0, Vec2.mk 1.5 1.0]
  let rotated := p.rotate (Float.pi / 2.0)
  let v0 := rotated.vertex 0
  ensure (floatNear v0.x 0.0 0.0001) "rotated x should be ~0"
  ensure (floatNear v0.y 1.0 0.0001) "rotated y should be 1"

test "scaleFromCenter preserves centroid" := do
  let p := Polygon2D.rectangle (Vec2.mk 2.0 2.0) (Vec2.mk 4.0 4.0)
  let c1 := p.centroid
  let scaled := p.scaleFromCenter 2.0
  let c2 := scaled.centroid
  ensure (floatNear c1.x c2.x 0.0001) "centroid.x should be preserved"
  ensure (floatNear c1.y c2.y 0.0001) "centroid.y should be preserved"

testSuite "Polygon2D Edge Queries"

test "edges returns correct count" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 2.0 2.0)
  let edges := p.edges
  ensure (edges.size == 4) "rectangle should have 4 edges"

test "closestPointOnBoundary" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 10.0 10.0)
  let closest := p.closestPointOnBoundary (Vec2.mk 5.0 15.0)
  ensure (floatNear closest.x 5.0 0.0001) "closest.x should be 5"
  ensure (floatNear closest.y 10.0 0.0001) "closest.y should be 10"

test "distanceToBoundary from outside" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 10.0 10.0)
  let d := p.distanceToBoundary (Vec2.mk 5.0 15.0)
  ensure (floatNear d 5.0 0.0001) "distance should be 5"

testSuite "Polygon2D Convex Hull"

test "convex hull of triangle is the triangle" := do
  let points := #[Vec2.mk 0.0 0.0, Vec2.mk 4.0 0.0, Vec2.mk 2.0 3.0]
  let hull := Polygon2D.convexHull points
  ensure (hull.vertexCount == 3) "hull should have 3 vertices"

test "convex hull excludes interior points" := do
  -- Square with center point
  let points := #[
    Vec2.mk 0.0 0.0,
    Vec2.mk 4.0 0.0,
    Vec2.mk 4.0 4.0,
    Vec2.mk 0.0 4.0,
    Vec2.mk 2.0 2.0  -- interior point
  ]
  let hull := Polygon2D.convexHull points
  ensure (hull.vertexCount == 4) "hull should have 4 vertices (interior point excluded)"

test "convex hull of many points on circle" := do
  -- 8 points on a circle plus some interior points
  let points := #[
    Vec2.mk 1.0 0.0, Vec2.mk 0.707 0.707, Vec2.mk 0.0 1.0, Vec2.mk (-0.707) 0.707,
    Vec2.mk (-1.0) 0.0, Vec2.mk (-0.707) (-0.707), Vec2.mk 0.0 (-1.0), Vec2.mk 0.707 (-0.707),
    Vec2.mk 0.0 0.0, Vec2.mk 0.1 0.1, Vec2.mk (-0.1) (-0.1)  -- interior
  ]
  let hull := Polygon2D.convexHull points
  ensure (hull.vertexCount == 8) "hull should have 8 vertices from circle"
  ensure hull.isConvex "hull should be convex"

test "convex hull is counter-clockwise" := do
  let points := #[Vec2.mk 0.0 0.0, Vec2.mk 4.0 0.0, Vec2.mk 4.0 4.0, Vec2.mk 0.0 4.0]
  let hull := Polygon2D.convexHull points
  ensure hull.isCounterClockwise "hull should be CCW"

test "toConvexHull on concave polygon" := do
  -- L-shaped polygon (concave)
  let p := Polygon2D.mk #[
    Vec2.mk 0.0 0.0,
    Vec2.mk 2.0 0.0,
    Vec2.mk 2.0 1.0,
    Vec2.mk 1.0 1.0,
    Vec2.mk 1.0 2.0,
    Vec2.mk 0.0 2.0
  ]
  let hull := p.toConvexHull
  ensure hull.isConvex "convex hull of L-shape should be convex"

testSuite "Polygon2D Triangulation"

test "triangle triangulates to itself" := do
  let p := Polygon2D.triangle (Vec2.mk 0.0 0.0) (Vec2.mk 4.0 0.0) (Vec2.mk 2.0 3.0)
  let triangles := p.triangulate
  ensure (triangles.size == 1) "triangle should produce 1 triangle"
  let tri := triangles[0]!
  ensure (tri.i0 == 0 && tri.i1 == 1 && tri.i2 == 2) "indices should be 0,1,2"

test "rectangle triangulates to 2 triangles" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 3.0)
  let triangles := p.triangulate
  ensure (triangles.size == 2) "rectangle should produce 2 triangles"

test "pentagon triangulates to 3 triangles" := do
  let p := Polygon2D.regular 5 1.0
  let triangles := p.triangulate
  ensure (triangles.size == 3) "pentagon should produce 3 triangles"

test "hexagon triangulates to 4 triangles" := do
  let p := Polygon2D.regular 6 1.0
  let triangles := p.triangulate
  ensure (triangles.size == 4) "hexagon should produce 4 triangles"

test "n-gon triangulates to n-2 triangles" := do
  for n in [3, 4, 5, 6, 7, 8] do
    let p := Polygon2D.regular n 1.0
    let triangles := p.triangulate
    ensure (triangles.size == n - 2) s!"regular {n}-gon should produce {n-2} triangles"

test "triangulated area equals original area" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 3.0)
  let triangles := p.triangulateToPolygons
  let triangleAreaSum := triangles.foldl (fun acc tri => acc + tri.area) 0.0
  ensure (floatNear triangleAreaSum p.area 0.001) "sum of triangle areas should equal polygon area"

test "concave L-shape triangulates correctly" := do
  -- L-shaped polygon (concave)
  let p := Polygon2D.mk #[
    Vec2.mk 0.0 0.0,
    Vec2.mk 2.0 0.0,
    Vec2.mk 2.0 1.0,
    Vec2.mk 1.0 1.0,
    Vec2.mk 1.0 2.0,
    Vec2.mk 0.0 2.0
  ]
  let triangles := p.triangulate
  ensure (triangles.size == 4) "6-vertex L-shape should produce 4 triangles"
  -- Check area is preserved
  let trianglePolys := p.triangulateToPolygons
  let triangleAreaSum := trianglePolys.foldl (fun acc tri => acc + tri.area) 0.0
  ensure (floatNear triangleAreaSum p.area 0.001) "triangle areas should sum to L-shape area"

test "triangulateToVertices returns vertex data" := do
  let p := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 3.0)
  let triangles := p.triangulateToVertices
  ensure (triangles.size == 2) "should have 2 triangles"
  ensure (triangles[0]!.size == 3) "each triangle should have 3 vertices"

#generate_tests

end LinalgTests.Polygon2DTests
