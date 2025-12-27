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

#generate_tests

end LinalgTests.TriangleTests
