/-
  Tests for Mat2 operations.
-/

import Linalg
import Crucible

namespace LinalgTests.Mat2Tests

open Crucible
open Linalg

testSuite "Mat2 Basic"

test "identity matrix diagonal is ones" := do
  let m := Mat2.identity
  m.get 0 0 ≡ 1.0
  m.get 1 1 ≡ 1.0

test "identity matrix off-diagonal is zeros" := do
  let m := Mat2.identity
  m.get 0 1 ≡ 0.0
  m.get 1 0 ≡ 0.0

test "identity times vector returns vector" := do
  let v := Vec2.mk 3.0 4.0
  let result := Mat2.identity * v
  result.x ≡ 3.0
  result.y ≡ 4.0

testSuite "Mat2 Operations"

test "multiply with identity returns same matrix" := do
  let m := Mat2.fromColumns (Vec2.mk 1.0 2.0) (Vec2.mk 3.0 4.0)
  let result := m * Mat2.identity
  result.get 0 0 ≡ m.get 0 0
  result.get 1 1 ≡ m.get 1 1

test "determinant of identity is 1" := do
  Mat2.identity.determinant ≡ 1.0

test "determinant of scaling matrix" := do
  let m := Mat2.scaling 2.0 3.0
  ensure (floatNear m.determinant 6.0 0.0001) "determinant should be 6"

test "inverse of identity is identity" := do
  match Mat2.identity.inverse with
  | some inv => ensure (inv.approxEq Mat2.identity) "inverse should equal identity"
  | none => ensure false "identity should be invertible"

test "rotation by 0 is identity" := do
  let m := Mat2.rotation 0.0
  ensure (m.approxEq Mat2.identity) "rotation by 0 should be identity"

#generate_tests

end LinalgTests.Mat2Tests
