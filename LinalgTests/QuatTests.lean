/-
  Tests for Quat operations.
-/

import Linalg
import Crucible

namespace LinalgTests.QuatTests

open Crucible
open Linalg

testSuite "Quat Basic"

test "identity quaternion has w=1 and xyz=0" := do
  let q := Quat.identity
  q.x ≡ 0.0
  q.y ≡ 0.0
  q.z ≡ 0.0
  q.w ≡ 1.0

test "identity quaternion is normalized" := do
  ensure (floatNear Quat.identity.length 1.0 0.0001) "length should be 1"

test "identity rotation leaves vector unchanged" := do
  let v := Vec3.mk 1.0 2.0 3.0
  let result := Quat.identity * v
  ensure (floatNear result.x 1.0 0.0001) "x should be 1"
  ensure (floatNear result.y 2.0 0.0001) "y should be 2"
  ensure (floatNear result.z 3.0 0.0001) "z should be 3"

testSuite "Quat Axis-Angle"

test "rotation around Y by 90 degrees rotates X to -Z" := do
  let q := Quat.fromAxisAngle Vec3.unitY Float.halfPi
  let v := Vec3.unitX
  let result := q * v
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y 0.0 0.0001) "y should be 0"
  ensure (floatNear result.z (-1.0) 0.0001) "z should be -1"

test "rotation around X by 90 degrees rotates Y to Z" := do
  let q := Quat.fromAxisAngle Vec3.unitX Float.halfPi
  let v := Vec3.unitY
  let result := q * v
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y 0.0 0.0001) "y should be 0"
  ensure (floatNear result.z 1.0 0.0001) "z should be 1"

testSuite "Quat Operations"

test "quaternion times its conjugate gives identity-like result" := do
  let q := Quat.fromAxisAngle Vec3.unitY (Float.pi / 3.0)
  let result := q * q.conjugate
  ensure (floatNear result.w 1.0 0.0001) "w should be 1"
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y 0.0 0.0001) "y should be 0"
  ensure (floatNear result.z 0.0 0.0001) "z should be 0"

test "slerp at 0 returns first quaternion" := do
  let a := Quat.identity
  let b := Quat.fromAxisAngle Vec3.unitY Float.halfPi
  let result := Quat.slerp a b 0.0
  ensure (floatNear result.w a.w 0.0001) "w should match"

test "slerp at 1 returns second quaternion" := do
  let a := Quat.identity
  let b := Quat.fromAxisAngle Vec3.unitY Float.halfPi
  let result := Quat.slerp a b 1.0
  ensure (result.sameRotation b) "should be same rotation"

test "toMat4 produces valid rotation matrix" := do
  let q := Quat.fromAxisAngle Vec3.unitY Float.halfPi
  let m := q.toMat4
  -- Determinant should be 1 for rotation matrix
  ensure (floatNear m.determinant 1.0 0.0001) "determinant should be 1"

#generate_tests

end LinalgTests.QuatTests
