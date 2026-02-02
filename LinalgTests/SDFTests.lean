/-
  Tests for signed distance field helpers.
-/

import Linalg
import Crucible

namespace LinalgTests.SDFTests

open Crucible
open Linalg
open Linalg.SDF

testSuite "SDF Primitives"

test "sphere sdf basic distances" := do
  let center := Vec3.zero
  let r := 1.0
  ensure (floatNear (sphere (Vec3.mk 0.0 0.0 0.0) center r) (-1.0) 0.0001) "inside negative"
  ensure (floatNear (sphere (Vec3.mk 1.0 0.0 0.0) center r) 0.0 0.0001) "on surface"
  ensure (floatNear (sphere (Vec3.mk 2.0 0.0 0.0) center r) 1.0 0.0001) "outside positive"

test "box sdf basic distances" := do
  let center := Vec3.zero
  let half := Vec3.mk 1.0 1.0 1.0
  let dInside := box (Vec3.mk 0.0 0.0 0.0) center half
  ensure (floatNear dInside (-1.0) 0.0001) "inside distance"
  let dOutside := box (Vec3.mk 2.0 0.0 0.0) center half
  ensure (floatNear dOutside 1.0 0.0001) "outside distance"

test "capsule sdf basic distances" := do
  let a := Vec3.mk 0.0 0.0 0.0
  let b := Vec3.mk 0.0 0.0 2.0
  let r := 0.5
  let dInside := capsule (Vec3.mk 0.0 0.0 1.0) a b r
  ensure (floatNear dInside (-0.5) 0.0001) "capsule inside"
  let dOutside := capsule (Vec3.mk 1.5 0.0 1.0) a b r
  ensure (floatNear dOutside 1.0 0.0001) "capsule outside"

test "circle/box2d sdf basic distances" := do
  let c := Vec2.zero
  let r := 1.0
  ensure (floatNear (circle (Vec2.mk 0.0 0.0) c r) (-1.0) 0.0001) "circle inside"
  ensure (floatNear (circle (Vec2.mk 1.0 0.0) c r) 0.0 0.0001) "circle surface"
  let half := Vec2.mk 1.0 1.0
  ensure (floatNear (box2D (Vec2.mk 0.0 0.0) c half) (-1.0) 0.0001) "box2d inside"

testSuite "SDF Operations"

test "union/intersection/subtract" := do
  let c1 := Vec3.mk (-1.0) 0.0 0.0
  let c2 := Vec3.mk (1.0) 0.0 0.0
  let r := 1.0
  let p := Vec3.zero
  let d1 := sphere p c1 r
  let d2 := sphere p c2 r
  let du := opUnion d1 d2
  let di := opIntersection d1 d2
  let ds := opSubtract d1 d2
  ensure (floatNear du (-0.0) 0.0001) "union"
  ensure (di >= du) "intersection should be >= union"
  ensure (ds >= d1) "subtract should be >= d1"

end LinalgTests.SDFTests
