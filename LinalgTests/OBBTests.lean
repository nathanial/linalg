/-
  Tests for OBB (Oriented Bounding Box).
-/

import Linalg
import Crucible

namespace LinalgTests.OBBTests

open Crucible
open Linalg

testSuite "OBB Basic"

test "OBB from center and extents" := do
  let obb := OBB.fromCenterExtents (Vec3.mk 1.0 2.0 3.0) (Vec3.mk 0.5 1.0 1.5)
  ensure (floatNear obb.center.x 1.0 0.0001) "center x"
  ensure (floatNear obb.halfExtents.x 0.5 0.0001) "halfExtents x"

test "axis-aligned OBB contains center" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  ensure (obb.containsPoint Vec3.zero) "center should be inside"

test "axis-aligned OBB contains point inside" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  ensure (obb.containsPoint (Vec3.mk 0.5 0.5 0.5)) "point should be inside"

test "axis-aligned OBB does not contain point outside" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  ensure (!obb.containsPoint (Vec3.mk 2.0 0.0 0.0)) "point should be outside"

test "rotated OBB contains point" := do
  -- OBB rotated 45 degrees around Y axis
  let angle := Float.pi / 4.0
  let obb : OBB := { center := Vec3.zero, halfExtents := Vec3.mk 1.0 1.0 1.0, orientation := Quat.fromAxisAngle Vec3.unitY angle }
  -- Point at (0.7, 0, 0.7) should be inside the rotated box
  ensure (obb.containsPoint (Vec3.mk 0.5 0.0 0.5)) "point should be inside rotated OBB"

test "OBB closest point on surface" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let closest := obb.closestPoint (Vec3.mk 3.0 0.0 0.0)
  ensure (floatNear closest.x 1.0 0.0001) "closest x should be 1"
  ensure (floatNear closest.y 0.0 0.0001) "closest y should be 0"

test "OBB distance to point outside" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let dist := obb.distance (Vec3.mk 3.0 0.0 0.0)
  ensure (floatNear dist 2.0 0.0001) "distance should be 2"

test "OBB distance to point inside is 0" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let dist := obb.distance (Vec3.mk 0.5 0.0 0.0)
  ensure (floatNear dist 0.0 0.0001) "distance should be 0"

test "OBB corners count" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let cs := obb.corners
  ensure (cs.size == 8) "should have 8 corners"

test "OBB to AABB" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let aabb := obb.toAABB
  ensure (floatNear aabb.min.x (-1.0) 0.0001) "min x"
  ensure (floatNear aabb.max.x 1.0 0.0001) "max x"

test "rotated OBB to AABB is larger" := do
  let angle := Float.pi / 4.0
  let obb : OBB := { center := Vec3.zero, halfExtents := Vec3.mk 1.0 1.0 1.0, orientation := Quat.fromAxisAngle Vec3.unitY angle }
  let aabb := obb.toAABB
  -- Rotated box should have larger AABB
  ensure (aabb.max.x > 1.0) "rotated AABB should be larger"

test "OBB volume" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 2.0 3.0)
  let vol := obb.volume
  ensure (floatNear vol 48.0 0.0001) "volume = 8 * 1 * 2 * 3 = 48"

testSuite "OBB Intersection"

test "OBB-OBB same box intersects" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  ensure (obb.intersectsOBB obb) "same OBB should intersect"

test "OBB-OBB overlapping intersects" := do
  let a := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let b := OBB.fromCenterExtents (Vec3.mk 1.0 0.0 0.0) (Vec3.mk 1.0 1.0 1.0)
  ensure (a.intersectsOBB b) "overlapping OBBs should intersect"

test "OBB-OBB separated does not intersect" := do
  let a := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let b := OBB.fromCenterExtents (Vec3.mk 5.0 0.0 0.0) (Vec3.mk 1.0 1.0 1.0)
  ensure (!a.intersectsOBB b) "separated OBBs should not intersect"

test "OBB-Sphere intersects" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  ensure (obb.intersectsSphere Vec3.zero 0.5) "sphere at center"
  ensure (obb.intersectsSphere (Vec3.mk 1.5 0.0 0.0) 1.0) "sphere touching"
  ensure (!obb.intersectsSphere (Vec3.mk 5.0 0.0 0.0) 1.0) "sphere far away"

testSuite "Ray-OBB"

test "ray hits axis-aligned OBB" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let ray := Ray.mk (Vec3.mk (-5.0) 0.0 0.0) Vec3.unitX
  match Intersection.rayOBB ray obb with
  | some hit =>
    ensure (floatNear hit.t 4.0 0.0001) "hit distance should be 4"
    ensure (floatNear hit.normal.x (-1.0) 0.0001) "normal should point -X"
  | none => ensure false "should hit"

test "ray misses OBB" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let ray := Ray.mk (Vec3.mk (-5.0) 5.0 0.0) Vec3.unitX
  match Intersection.rayOBB ray obb with
  | some _ => ensure false "should miss"
  | none => ensure true "correctly missed"

test "ray hits rotated OBB" := do
  let angle := Float.pi / 4.0
  let obb : OBB := { center := Vec3.zero, halfExtents := Vec3.mk 1.0 1.0 1.0, orientation := Quat.fromAxisAngle Vec3.unitY angle }
  let ray := Ray.mk (Vec3.mk (-5.0) 0.0 0.0) Vec3.unitX
  match Intersection.rayOBB ray obb with
  | some hit =>
    ensure (hit.t > 0.0) "should have positive t"
    ensure (hit.point.x < 0.0) "hit point should be on -X side"
  | none => ensure false "should hit rotated OBB"

test "ray inside OBB" := do
  let obb := OBB.fromCenterExtents Vec3.zero (Vec3.mk 2.0 2.0 2.0)
  let ray := Ray.mk Vec3.zero Vec3.unitX
  match Intersection.rayOBB ray obb with
  | some hit =>
    ensure (hit.t >= 0.0) "t should be non-negative"
  | none => ensure false "ray from inside should hit"



end LinalgTests.OBBTests
