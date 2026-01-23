/-
  Tests for Capsule primitive.
-/

import Linalg
import Crucible

namespace LinalgTests.CapsuleTests

open Crucible
open Linalg

testSuite "Capsule Basic"

test "capsule from endpoints" := do
  let c := Capsule.fromEndpoints (Vec3.mk 0.0 0.0 0.0) (Vec3.mk 0.0 2.0 0.0) 0.5
  ensure (floatNear c.a.y 0.0 0.0001) "endpoint a"
  ensure (floatNear c.b.y 2.0 0.0001) "endpoint b"
  ensure (floatNear c.radius 0.5 0.0001) "radius"

test "vertical capsule" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  ensure (floatNear c.a.y (-1.0) 0.0001) "bottom endpoint"
  ensure (floatNear c.b.y 1.0 0.0001) "top endpoint"

test "capsule center" := do
  let c := Capsule.fromEndpoints (Vec3.mk 0.0 0.0 0.0) (Vec3.mk 0.0 4.0 0.0) 1.0
  let center := c.center
  ensure (floatNear center.y 2.0 0.0001) "center y"

test "capsule segment length" := do
  let c := Capsule.fromEndpoints Vec3.zero (Vec3.mk 3.0 4.0 0.0) 1.0
  let len := c.segmentLength
  ensure (floatNear len 5.0 0.0001) "segment length = 5"

test "capsule total height" := do
  let c := Capsule.fromEndpoints Vec3.zero (Vec3.mk 0.0 2.0 0.0) 0.5
  let h := c.totalHeight
  ensure (floatNear h 3.0 0.0001) "total height = segment + 2*radius"

test "capsule contains center point" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  ensure (c.containsPoint Vec3.zero) "center should be inside"

test "capsule contains point on axis" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  ensure (c.containsPoint (Vec3.mk 0.0 0.5 0.0)) "point on axis should be inside"

test "capsule contains point off axis within radius" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  ensure (c.containsPoint (Vec3.mk 0.3 0.0 0.0)) "point within radius should be inside"

test "capsule does not contain point outside" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  ensure (!c.containsPoint (Vec3.mk 2.0 0.0 0.0)) "point far away should be outside"

test "capsule contains point in hemisphere" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  ensure (c.containsPoint (Vec3.mk 0.0 1.3 0.0)) "point in top hemisphere"

test "closest point on segment" := do
  let c := Capsule.fromEndpoints Vec3.zero (Vec3.mk 4.0 0.0 0.0) 1.0
  let closest := c.closestPointOnSegment (Vec3.mk 2.0 3.0 0.0)
  ensure (floatNear closest.x 2.0 0.0001) "closest x"
  ensure (floatNear closest.y 0.0 0.0001) "closest y"

test "closest point clamped to endpoint" := do
  let c := Capsule.fromEndpoints Vec3.zero (Vec3.mk 4.0 0.0 0.0) 1.0
  let closest := c.closestPointOnSegment (Vec3.mk (-2.0) 0.0 0.0)
  ensure (floatNear closest.x 0.0 0.0001) "should clamp to start"

test "capsule distance to point" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  let dist := c.distance (Vec3.mk 2.0 0.0 0.0)
  ensure (floatNear dist 1.5 0.0001) "distance = 2 - 0.5"

test "capsule signed distance inside" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  let sd := c.signedDistance Vec3.zero
  ensure (sd < 0.0) "signed distance inside should be negative"

test "capsule to AABB" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  let aabb := c.toAABB
  ensure (floatNear aabb.min.y (-1.5) 0.0001) "min y includes radius"
  ensure (floatNear aabb.max.y 1.5 0.0001) "max y includes radius"

test "capsule volume" := do
  let c := Capsule.fromEndpoints Vec3.zero (Vec3.mk 0.0 2.0 0.0) 1.0
  let vol := c.volume
  -- Volume = π*r²*h + (4/3)*π*r³ = π*1*2 + (4/3)*π*1 = π*(2 + 4/3) = π*10/3 ≈ 10.47
  ensure (floatNear vol (Float.pi * 10.0 / 3.0) 0.01) "volume formula"

testSuite "Capsule Intersection"

test "capsule-sphere intersects" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  let s := Sphere.mk Vec3.zero 0.5
  ensure (c.intersectsSphere s) "sphere at center"

test "capsule-sphere touching" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  let s := Sphere.mk (Vec3.mk 1.0 0.0 0.0) 0.5
  ensure (c.intersectsSphere s) "sphere touching"

test "capsule-sphere not intersecting" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  let s := Sphere.mk (Vec3.mk 5.0 0.0 0.0) 0.5
  ensure (!c.intersectsSphere s) "sphere far away"

test "capsule-capsule parallel intersects" := do
  let c1 := Capsule.vertical Vec3.zero 1.0 0.5
  let c2 := Capsule.vertical (Vec3.mk 0.8 0.0 0.0) 1.0 0.5
  ensure (c1.intersectsCapsule c2) "parallel capsules overlapping"

test "capsule-capsule parallel not intersecting" := do
  let c1 := Capsule.vertical Vec3.zero 1.0 0.5
  let c2 := Capsule.vertical (Vec3.mk 5.0 0.0 0.0) 1.0 0.5
  ensure (!c1.intersectsCapsule c2) "parallel capsules far apart"

test "capsule-capsule perpendicular intersects" := do
  let c1 := Capsule.vertical Vec3.zero 1.0 0.5
  let c2 := Capsule.fromEndpoints (Vec3.mk (-1.0) 0.0 0.0) (Vec3.mk 1.0 0.0 0.0) 0.5
  ensure (c1.intersectsCapsule c2) "perpendicular capsules"

testSuite "Ray-Capsule"

test "ray hits capsule cylinder body" := do
  let c := Capsule.vertical Vec3.zero 2.0 1.0
  let ray := Ray.mk (Vec3.mk (-5.0) 0.0 0.0) Vec3.unitX
  match Intersection.rayCapsule ray c with
  | some hit =>
    ensure (floatNear hit.t 4.0 0.0001) "hit distance should be 4"
    ensure (floatNear hit.normal.x (-1.0) 0.0001) "normal should point -X"
  | none => ensure false "should hit cylinder body"

test "ray hits capsule hemisphere" := do
  let c := Capsule.vertical Vec3.zero 1.0 1.0
  let ray := Ray.mk (Vec3.mk 0.0 5.0 0.0) Vec3.down
  match Intersection.rayCapsule ray c with
  | some hit =>
    ensure (hit.t > 0.0) "should have positive t"
    ensure (hit.point.y > 1.0) "hit point should be on top hemisphere"
  | none => ensure false "should hit top hemisphere"

test "ray misses capsule" := do
  let c := Capsule.vertical Vec3.zero 1.0 0.5
  let ray := Ray.mk (Vec3.mk (-5.0) 5.0 0.0) Vec3.unitX
  match Intersection.rayCapsule ray c with
  | some _ => ensure false "should miss"
  | none => ensure true "correctly missed"

test "ray parallel to capsule axis hits" := do
  let c := Capsule.vertical Vec3.zero 2.0 1.0
  let ray := Ray.mk (Vec3.mk 0.5 (-5.0) 0.0) Vec3.unitY
  match Intersection.rayCapsule ray c with
  | some hit =>
    ensure (hit.t > 0.0) "should have positive t"
  | none => ensure false "parallel ray should hit"



end LinalgTests.CapsuleTests
