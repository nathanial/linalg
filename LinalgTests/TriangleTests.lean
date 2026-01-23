/-
  Tests for Triangle primitive.
-/

import Linalg
import Crucible

namespace LinalgTests.TriangleTests

open Crucible
open Linalg

testSuite "Triangle Basic"

test "triangle area calculation" := do
  -- Right triangle with legs of length 1
  let tri := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  ensure (floatNear tri.area 0.5 0.0001) "area should be 0.5"

test "triangle centroid" := do
  let tri := Triangle.mk' Vec3.zero (Vec3.mk 3.0 0.0 0.0) (Vec3.mk 0.0 3.0 0.0)
  let c := tri.centroid
  ensure (floatNear c.x 1.0 0.0001) "centroid x should be 1"
  ensure (floatNear c.y 1.0 0.0001) "centroid y should be 1"
  ensure (floatNear c.z 0.0 0.0001) "centroid z should be 0"

test "triangle normal pointing up for XY plane triangle" := do
  let tri := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  let n := tri.unitNormal
  ensure (floatNear n.x 0.0 0.0001) "normal x should be 0"
  ensure (floatNear n.y 0.0 0.0001) "normal y should be 0"
  ensure (floatNear n.z 1.0 0.0001) "normal z should be 1"

testSuite "Triangle Barycentric"

test "barycentric of vertex v0" := do
  let tri := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  let bc := tri.barycentric Vec3.zero
  ensure (floatNear bc.u 1.0 0.0001) "u should be 1"
  ensure (floatNear bc.v 0.0 0.0001) "v should be 0"
  ensure (floatNear bc.w 0.0 0.0001) "w should be 0"

test "barycentric of centroid" := do
  let tri := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  let c := tri.centroid
  let bc := tri.barycentric c
  ensure (floatNear bc.u 0.3333 0.01) "u should be ~0.333"
  ensure (floatNear bc.v 0.3333 0.01) "v should be ~0.333"
  ensure (floatNear bc.w 0.3333 0.01) "w should be ~0.333"

test "point inside triangle" := do
  let tri := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  ensure (tri.containsPoint (Vec3.mk 0.25 0.25 0.0)) "should contain point"

test "point outside triangle" := do
  let tri := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  ensure (!(tri.containsPoint (Vec3.mk 1.0 1.0 0.0))) "should not contain point"

testSuite "Ray-Triangle Intersection"

test "ray hits triangle" := do
  let tri := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  let ray := Ray.mk' (Vec3.mk 0.25 0.25 (-1.0)) Vec3.unitZ
  match Intersection.rayTriangle ray tri with
  | some hit => ensure (floatNear hit.t 1.0 0.0001) "t should be 1"
  | none => ensure false "expected ray to hit triangle"

test "ray misses triangle" := do
  let tri := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  let ray := Ray.mk' (Vec3.mk 2.0 2.0 (-1.0)) Vec3.unitZ
  match Intersection.rayTriangle ray tri with
  | some _ => ensure false "expected ray to miss triangle"
  | none => pure ()

test "ray parallel to triangle misses" := do
  let tri := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  let ray := Ray.mk' (Vec3.mk 0.0 0.0 1.0) Vec3.unitX
  match Intersection.rayTriangle ray tri with
  | some _ => ensure false "expected parallel ray to miss"
  | none => pure ()

testSuite "Triangle-Triangle Intersection"

test "overlapping triangles intersect" := do
  -- Two triangles in XY plane that share the origin
  let tri1 := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  let tri2 := Triangle.mk' Vec3.zero (Vec3.mk 0.5 0.5 0.0) (Vec3.mk (-0.5) 0.5 0.0)
  ensure (Intersection.triangleTriangle tri1 tri2) "overlapping triangles should intersect"

test "same triangle intersects itself" := do
  let tri := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  ensure (Intersection.triangleTriangle tri tri) "triangle should intersect itself"

test "separated triangles do not intersect" := do
  let tri1 := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  let tri2 := Triangle.mk' (Vec3.mk 10.0 0.0 0.0) (Vec3.mk 11.0 0.0 0.0) (Vec3.mk 10.0 1.0 0.0)
  ensure (!Intersection.triangleTriangle tri1 tri2) "separated triangles should not intersect"

test "triangles in parallel planes do not intersect" := do
  let tri1 := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  let tri2 := Triangle.mk' (Vec3.mk 0.0 0.0 5.0) (Vec3.mk 1.0 0.0 5.0) (Vec3.mk 0.0 1.0 5.0)
  ensure (!Intersection.triangleTriangle tri1 tri2) "parallel plane triangles should not intersect"

test "crossing triangles intersect" := do
  -- Triangle in XY plane
  let tri1 := Triangle.mk' (Vec3.mk (-1.0) (-1.0) 0.0) (Vec3.mk 1.0 (-1.0) 0.0) (Vec3.mk 0.0 1.0 0.0)
  -- Triangle in XZ plane, crossing through center
  let tri2 := Triangle.mk' (Vec3.mk (-1.0) 0.0 (-1.0)) (Vec3.mk 1.0 0.0 (-1.0)) (Vec3.mk 0.0 0.0 1.0)
  ensure (Intersection.triangleTriangle tri1 tri2) "crossing triangles should intersect"

test "edge-touching triangles intersect" := do
  -- Two triangles sharing an edge
  let tri1 := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  let tri2 := Triangle.mk' Vec3.zero Vec3.unitX (Vec3.mk 0.5 0.5 1.0)
  ensure (Intersection.triangleTriangle tri1 tri2) "edge-touching triangles should intersect"

test "vertex-touching triangles intersect" := do
  -- Two triangles sharing only a vertex
  let tri1 := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  let tri2 := Triangle.mk' Vec3.zero (Vec3.mk (-1.0) 0.0 0.0) (Vec3.mk 0.0 (-1.0) 0.0)
  ensure (Intersection.triangleTriangle tri1 tri2) "vertex-touching triangles should intersect"

test "coplanar non-overlapping triangles do not intersect" := do
  let tri1 := Triangle.mk' Vec3.zero Vec3.unitX Vec3.unitY
  let tri2 := Triangle.mk' (Vec3.mk 5.0 0.0 0.0) (Vec3.mk 6.0 0.0 0.0) (Vec3.mk 5.0 1.0 0.0)
  ensure (!Intersection.triangleTriangle tri1 tri2) "coplanar non-overlapping should not intersect"

test "triangleTriangleContact returns intersection segment" := do
  -- Triangle in XY plane
  let tri1 := Triangle.mk' (Vec3.mk (-1.0) (-1.0) 0.0) (Vec3.mk 1.0 (-1.0) 0.0) (Vec3.mk 0.0 1.0 0.0)
  -- Triangle in XZ plane
  let tri2 := Triangle.mk' (Vec3.mk (-1.0) 0.0 (-1.0)) (Vec3.mk 1.0 0.0 (-1.0)) (Vec3.mk 0.0 0.0 1.0)
  match Intersection.triangleTriangleContact tri1 tri2 with
  | some (p1, p2) =>
    -- Both points should be on the intersection line (y=0, z=0 in this case)
    ensure (floatNear p1.y 0.0 0.01) "p1.y should be ~0"
    ensure (floatNear p1.z 0.0 0.01) "p1.z should be ~0"
  | none => ensure false "expected contact info"



end LinalgTests.TriangleTests
