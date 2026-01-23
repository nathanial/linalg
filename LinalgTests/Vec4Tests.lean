/-
  Tests for Vec4 operations.
-/

import Linalg
import Crucible

namespace LinalgTests.Vec4Tests

open Crucible
open Linalg

testSuite "Vec4 Basic"

test "zero vector has all zero components" := do
  Vec4.zero.x ≡ 0.0
  Vec4.zero.y ≡ 0.0
  Vec4.zero.z ≡ 0.0
  Vec4.zero.w ≡ 0.0

test "addition works correctly" := do
  let a := Vec4.mk 1.0 2.0 3.0 4.0
  let b := Vec4.mk 5.0 6.0 7.0 8.0
  let c := a + b
  c.x ≡ 6.0
  c.y ≡ 8.0
  c.z ≡ 10.0
  c.w ≡ 12.0

testSuite "Vec4 Homogeneous"

test "fromPoint sets w to 1" := do
  let v3 := Vec3.mk 1.0 2.0 3.0
  let v4 := Vec4.fromPoint v3
  v4.x ≡ 1.0
  v4.y ≡ 2.0
  v4.z ≡ 3.0
  v4.w ≡ 1.0

test "fromDirection sets w to 0" := do
  let v3 := Vec3.mk 1.0 2.0 3.0
  let v4 := Vec4.fromDirection v3
  v4.x ≡ 1.0
  v4.y ≡ 2.0
  v4.z ≡ 3.0
  v4.w ≡ 0.0

test "toVec3Normalized performs perspective divide" := do
  let v := Vec4.mk 2.0 4.0 6.0 2.0
  let v3 := v.toVec3Normalized
  ensure (floatNear v3.x 1.0 0.0001) "x should be 1"
  ensure (floatNear v3.y 2.0 0.0001) "y should be 2"
  ensure (floatNear v3.z 3.0 0.0001) "z should be 3"



end LinalgTests.Vec4Tests
