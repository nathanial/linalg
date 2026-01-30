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

testSuite "Quat SQUAD"

test "log and exp are inverses" := do
  let q := (Quat.fromAxisAngle (Vec3.mk 1 1 1).normalize (Float.pi / 3.0)).normalize
  let roundTrip := q.log.exp
  ensure (roundTrip.sameRotation q) "exp(log(q)) should equal q"

test "log of identity is zero" := do
  let result := Quat.identity.log
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y 0.0 0.0001) "y should be 0"
  ensure (floatNear result.z 0.0 0.0001) "z should be 0"
  ensure (floatNear result.w 0.0 0.0001) "w should be 0"

test "exp of zero is identity" := do
  let result := (Quat.mk 0 0 0 0).exp
  ensure (result.sameRotation Quat.identity) "exp(0) should be identity"

test "squad at t=0 returns first quaternion" := do
  let q1 := Quat.identity
  let q2 := Quat.fromAxisAngle Vec3.unitY Float.halfPi
  let s1 := q1
  let s2 := q2
  let result := Quat.squad q1 q2 s1 s2 0.0
  ensure (result.sameRotation q1) "t=0 should match q1"

test "squad at t=1 returns second quaternion" := do
  let q1 := Quat.identity
  let q2 := Quat.fromAxisAngle Vec3.unitY Float.halfPi
  let s1 := q1
  let s2 := q2
  let result := Quat.squad q1 q2 s1 s2 1.0
  ensure (result.sameRotation q2) "t=1 should match q2"

test "squad midpoint is normalized" := do
  let q1 := Quat.identity
  let q2 := Quat.fromAxisAngle Vec3.unitY Float.pi
  let s1 := Quat.squadIntermediate q1 q1 q2
  let s2 := Quat.squadIntermediate q1 q2 q2
  let mid := Quat.squad q1 q2 s1 s2 0.5
  ensure (floatNear mid.length 1.0 0.0001) "midpoint should be normalized"

test "squadPath endpoints match keyframes" := do
  let k0 := Quat.identity
  let k1 := Quat.fromAxisAngle Vec3.unitY Float.halfPi
  let k2 := Quat.fromAxisAngle Vec3.unitY Float.pi
  let keyframes := #[k0, k1, k2]
  ensure ((Quat.squadPath keyframes 0.0).sameRotation k0) "t=0 should match k0"
  ensure ((Quat.squadPath keyframes 1.0).sameRotation k1) "t=1 should match k1"

test "squadPath with 2 keyframes at midpoint" := do
  let q1 := Quat.identity
  let q2 := Quat.fromAxisAngle Vec3.unitY Float.pi
  let path := #[q1, q2]
  let squadResult := Quat.squadPath path 0.5
  let slerpResult := Quat.slerp q1 q2 0.5
  ensure (squadResult.sameRotation slerpResult) "should match slerp"

end LinalgTests.QuatTests
