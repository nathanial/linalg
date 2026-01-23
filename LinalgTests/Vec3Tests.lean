/-
  Tests for Vec3 operations.
-/

import Linalg
import Crucible

namespace LinalgTests.Vec3Tests

open Crucible
open Linalg

testSuite "Vec3 Basic"

test "zero vector has all zero components" := do
  Vec3.zero.x ≡ 0.0
  Vec3.zero.y ≡ 0.0
  Vec3.zero.z ≡ 0.0

test "addition works correctly" := do
  let a := Vec3.mk 1.0 2.0 3.0
  let b := Vec3.mk 4.0 5.0 6.0
  let c := a + b
  c.x ≡ 5.0
  c.y ≡ 7.0
  c.z ≡ 9.0

test "scalar multiplication works" := do
  let a := Vec3.mk 1.0 2.0 3.0
  let b := a * 2.0
  b.x ≡ 2.0
  b.y ≡ 4.0
  b.z ≡ 6.0

testSuite "Vec3 Math"

test "dot product of perpendicular vectors is zero" := do
  let a := Vec3.unitX
  let b := Vec3.unitY
  ensure (floatNear (a.dot b) 0.0 0.0001) "perpendicular dot should be 0"

test "dot product of parallel vectors is product of lengths" := do
  let a := Vec3.mk 2.0 0.0 0.0
  let b := Vec3.mk 3.0 0.0 0.0
  ensure (floatNear (a.dot b) 6.0 0.0001) "parallel dot should be 6"

test "cross product of X and Y is Z" := do
  let c := Vec3.unitX.cross Vec3.unitY
  ensure (floatNear c.x 0.0 0.0001) "cross x should be 0"
  ensure (floatNear c.y 0.0 0.0001) "cross y should be 0"
  ensure (floatNear c.z 1.0 0.0001) "cross z should be 1"

test "cross product of Y and X is -Z" := do
  let c := Vec3.unitY.cross Vec3.unitX
  ensure (floatNear c.x 0.0 0.0001) "cross x should be 0"
  ensure (floatNear c.y 0.0 0.0001) "cross y should be 0"
  ensure (floatNear c.z (-1.0) 0.0001) "cross z should be -1"

test "normalize produces unit length" := do
  let v := Vec3.mk 3.0 4.0 0.0
  let n := v.normalize
  ensure (floatNear n.length 1.0 0.0001) "normalized length should be 1"

test "distance between points" := do
  let a := Vec3.mk 0.0 0.0 0.0
  let b := Vec3.mk 3.0 4.0 0.0
  ensure (floatNear (a.distance b) 5.0 0.0001) "distance should be 5"

test "reflect off horizontal surface" := do
  let v := Vec3.mk 1.0 (-1.0) 0.0
  let n := Vec3.unitY
  let r := v.reflect n
  ensure (floatNear r.x 1.0 0.0001) "reflected x should be 1"
  ensure (floatNear r.y 1.0 0.0001) "reflected y should be 1"
  ensure (floatNear r.z 0.0 0.0001) "reflected z should be 0"



end LinalgTests.Vec3Tests
