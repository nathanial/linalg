/-
  Tests for Mat4 operations.
-/

import Linalg
import Crucible

namespace LinalgTests.Mat4Tests

open Crucible
open Linalg

testSuite "Mat4 Basic"

test "identity matrix works correctly" := do
  let m := Mat4.identity
  m.get 0 0 ≡ 1.0
  m.get 1 1 ≡ 1.0
  m.get 2 2 ≡ 1.0
  m.get 3 3 ≡ 1.0

test "identity times vec4 returns vec4" := do
  let v := Vec4.mk 1.0 2.0 3.0 4.0
  let result := Mat4.identity * v
  result.x ≡ 1.0
  result.y ≡ 2.0
  result.z ≡ 3.0
  result.w ≡ 4.0

testSuite "Mat4 Transforms"

test "translation matrix moves point" := do
  let m := Mat4.translation 10.0 20.0 30.0
  let p := Vec3.zero
  let result := m.transformPoint p
  ensure (floatNear result.x 10.0 0.0001) "x should be 10"
  ensure (floatNear result.y 20.0 0.0001) "y should be 20"
  ensure (floatNear result.z 30.0 0.0001) "z should be 30"

test "translation does not affect direction" := do
  let m := Mat4.translation 10.0 20.0 30.0
  let d := Vec3.unitX
  let result := m.transformDirection d
  ensure (floatNear result.x 1.0 0.0001) "x should be 1"
  ensure (floatNear result.y 0.0 0.0001) "y should be 0"
  ensure (floatNear result.z 0.0 0.0001) "z should be 0"

test "scaling matrix scales point" := do
  let m := Mat4.scaling 2.0 3.0 4.0
  let p := Vec3.mk 1.0 1.0 1.0
  let result := m.transformPoint p
  ensure (floatNear result.x 2.0 0.0001) "x should be 2"
  ensure (floatNear result.y 3.0 0.0001) "y should be 3"
  ensure (floatNear result.z 4.0 0.0001) "z should be 4"

test "rotation X by 90 degrees rotates Y to Z" := do
  let m := Mat4.rotationX Float.halfPi
  let v := Vec3.unitY
  let result := m.transformDirection v
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y 0.0 0.0001) "y should be 0"
  ensure (floatNear result.z 1.0 0.0001) "z should be 1"

test "getTranslation extracts translation" := do
  let m := Mat4.translation 5.0 10.0 15.0
  let t := m.getTranslation
  ensure (floatNear t.x 5.0 0.0001) "x should be 5"
  ensure (floatNear t.y 10.0 0.0001) "y should be 10"
  ensure (floatNear t.z 15.0 0.0001) "z should be 15"

testSuite "Mat4 Inverse"

test "identity inverse is identity" := do
  match Mat4.identity.inverse with
  | some inv => ensure (inv.approxEq Mat4.identity) "inverse should equal identity"
  | none => ensure false "identity should be invertible"

test "translation inverse negates translation" := do
  let m := Mat4.translation 5.0 10.0 15.0
  match m.inverse with
  | some inv =>
    let t := inv.getTranslation
    ensure (floatNear t.x (-5.0) 0.0001) "x should be -5"
    ensure (floatNear t.y (-10.0) 0.0001) "y should be -10"
    ensure (floatNear t.z (-15.0) 0.0001) "z should be -15"
  | none => ensure false "translation should be invertible"

#generate_tests

end LinalgTests.Mat4Tests
