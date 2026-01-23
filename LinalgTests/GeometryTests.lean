/-
  Tests for geometry primitives and intersection tests.
-/

import Linalg
import Crucible

namespace LinalgTests.GeometryTests

open Crucible
open Linalg

testSuite "Ray"

test "ray pointAt works correctly" := do
  let ray := Ray.mk' Vec3.zero Vec3.unitZ
  let p := ray.pointAt 5.0
  ensure (floatNear p.z 5.0 0.0001) "z should be 5"

testSuite "AABB"

test "AABB center is correct" := do
  let aabb := AABB.fromMinMax (Vec3.mk 0.0 0.0 0.0) (Vec3.mk 10.0 10.0 10.0)
  let c := aabb.center
  ensure (floatNear c.x 5.0 0.0001) "center x should be 5"
  ensure (floatNear c.y 5.0 0.0001) "center y should be 5"
  ensure (floatNear c.z 5.0 0.0001) "center z should be 5"

test "AABB contains center point" := do
  let aabb := AABB.fromMinMax Vec3.zero (Vec3.mk 10.0 10.0 10.0)
  ensure (aabb.containsPoint (Vec3.mk 5.0 5.0 5.0)) "should contain center"

test "AABB does not contain point outside" := do
  let aabb := AABB.fromMinMax Vec3.zero (Vec3.mk 10.0 10.0 10.0)
  ensure (!(aabb.containsPoint (Vec3.mk 15.0 5.0 5.0))) "should not contain outside point"

testSuite "Sphere"

test "sphere contains its center" := do
  let s := Sphere.mk' (Vec3.mk 5.0 5.0 5.0) 3.0
  ensure (s.containsPoint s.center) "should contain center"

test "sphere does not contain distant point" := do
  let s := Sphere.mk' Vec3.zero 1.0
  ensure (!(s.containsPoint (Vec3.mk 10.0 0.0 0.0))) "should not contain distant point"

testSuite "Plane"

test "point on XY plane has zero distance" := do
  let p := Plane.xy
  let point := Vec3.mk 5.0 5.0 0.0
  ensure (floatNear (p.distanceToPoint point) 0.0 0.0001) "distance should be 0"

test "point above XY plane has positive signed distance" := do
  let p := Plane.xy
  let point := Vec3.mk 0.0 0.0 5.0
  ensure (floatNear (p.signedDistance point) 5.0 0.0001) "signed distance should be 5"

testSuite "Intersection"

test "ray hits sphere from outside" := do
  let ray := Ray.mk' (Vec3.mk 0.0 0.0 (-10.0)) Vec3.unitZ
  let sphere := Sphere.mk' Vec3.zero 1.0
  match Intersection.raySphere ray sphere with
  | some hit => ensure (floatNear hit.t 9.0 0.0001) "t should be 9"
  | none => ensure false "expected ray to hit sphere"

test "ray misses sphere" := do
  let ray := Ray.mk' (Vec3.mk 10.0 0.0 (-10.0)) Vec3.unitZ
  let sphere := Sphere.mk' Vec3.zero 1.0
  match Intersection.raySphere ray sphere with
  | some _ => ensure false "expected ray to miss sphere"
  | none => pure ()

test "ray hits AABB" := do
  let ray := Ray.mk' (Vec3.mk 0.0 0.0 (-10.0)) Vec3.unitZ
  let aabb := AABB.fromMinMax (Vec3.mk (-1.0) (-1.0) (-1.0)) (Vec3.mk 1.0 1.0 1.0)
  match Intersection.rayAABB ray aabb with
  | some (tMin, _) => ensure (floatNear tMin 9.0 0.0001) "tMin should be 9"
  | none => ensure false "expected ray to hit AABB"

test "ray hits plane" := do
  let ray := Ray.mk' (Vec3.mk 0.0 0.0 (-5.0)) Vec3.unitZ
  let plane := Plane.xy
  match Intersection.rayPlane ray plane with
  | some hit => ensure (floatNear hit.t 5.0 0.0001) "t should be 5"
  | none => ensure false "expected ray to hit plane"

test "sphere-sphere intersection" := do
  let a := Sphere.mk' Vec3.zero 2.0
  let b := Sphere.mk' (Vec3.mk 3.0 0.0 0.0) 2.0
  ensure (Intersection.sphereSphere a b) "spheres should intersect"

test "sphere-sphere no intersection" := do
  let a := Sphere.mk' Vec3.zero 1.0
  let b := Sphere.mk' (Vec3.mk 10.0 0.0 0.0) 1.0
  ensure (!(Intersection.sphereSphere a b)) "spheres should not intersect"

test "AABB-AABB intersection" := do
  let a := AABB.fromMinMax Vec3.zero (Vec3.mk 5.0 5.0 5.0)
  let b := AABB.fromMinMax (Vec3.mk 3.0 3.0 3.0) (Vec3.mk 8.0 8.0 8.0)
  ensure (Intersection.aabbAABB a b) "AABBs should intersect"

test "AABB-AABB no intersection" := do
  let a := AABB.fromMinMax Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let b := AABB.fromMinMax (Vec3.mk 5.0 5.0 5.0) (Vec3.mk 6.0 6.0 6.0)
  ensure (!(Intersection.aabbAABB a b)) "AABBs should not intersect"



end LinalgTests.GeometryTests
