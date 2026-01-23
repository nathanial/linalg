/-
  Tests for Euler angle utilities.
-/

import Linalg
import Crucible

namespace LinalgTests.EulerTests

open Crucible
open Linalg

testSuite "Euler Constructors"

test "xyz creates correct angles" := do
  let e := Euler.xyz 0.1 0.2 0.3
  ensure (floatNear e.a1 0.1 0.0001) "a1 should be 0.1"
  ensure (floatNear e.a2 0.2 0.0001) "a2 should be 0.2"
  ensure (floatNear e.a3 0.3 0.0001) "a3 should be 0.3"
  ensure (e.order == EulerOrder.XYZ) "order should be XYZ"

test "yxz creates correct angles" := do
  let e := Euler.yxz 0.1 0.2 0.3
  ensure (e.order == EulerOrder.YXZ) "order should be YXZ"

test "zyx creates correct angles" := do
  let e := Euler.zyx 0.1 0.2 0.3
  ensure (e.order == EulerOrder.ZYX) "order should be ZYX"

test "fromDegrees converts correctly" := do
  let e := Euler.fromDegrees 90.0 0.0 0.0
  ensure (floatNear e.a1 Float.halfPi 0.0001) "90 degrees should be π/2"

test "toDegrees converts correctly" := do
  let e := Euler.xyz Float.halfPi 0.0 0.0
  let deg := e.toDegrees
  ensure (floatNear deg.x 90.0 0.0001) "π/2 should be 90 degrees"

testSuite "Euler to Quaternion"

test "identity euler gives identity quaternion" := do
  let e := Euler.zero
  let q := e.toQuat
  ensure (q.sameRotation Quat.identity) "zero euler should give identity quat"

test "90 degree X rotation" := do
  let e := Euler.xyz Float.halfPi 0.0 0.0
  let q := e.toQuat
  -- Rotate unit Y around X by 90° should give unit Z
  let v := q * Vec3.unitY
  ensure (floatNear v.x 0.0 0.01) "x should be 0"
  ensure (floatNear v.y 0.0 0.01) "y should be 0"
  ensure (floatNear v.z 1.0 0.01) "z should be 1"

test "90 degree Y rotation" := do
  let e := Euler.xyz 0.0 Float.halfPi 0.0
  let q := e.toQuat
  -- Rotate unit X around Y by 90° should give -unit Z
  let v := q * Vec3.unitX
  ensure (floatNear v.x 0.0 0.01) "x should be 0"
  ensure (floatNear v.z (-1.0) 0.01) "z should be -1"

test "90 degree Z rotation" := do
  let e := Euler.xyz 0.0 0.0 Float.halfPi
  let q := e.toQuat
  -- Rotate unit X around Z by 90° should give unit Y
  let v := q * Vec3.unitX
  ensure (floatNear v.x 0.0 0.01) "x should be 0"
  ensure (floatNear v.y 1.0 0.01) "y should be 1"

testSuite "Euler Round Trip"

test "XYZ round trip through quaternion produces same rotation" := do
  let e := Euler.xyz 0.3 0.2 0.1
  let q := e.toQuat
  let e2 := Euler.fromQuat q .XYZ
  let q2 := e2.toQuat
  -- Compare resulting quaternions (rotations), not raw angles
  ensure (q.sameRotation q2 0.01) "rotations should match"

test "YXZ round trip through quaternion produces same rotation" := do
  let e := Euler.yxz 0.3 0.2 0.1
  let q := e.toQuat
  let e2 := Euler.fromQuat q .YXZ
  let q2 := e2.toQuat
  ensure (q.sameRotation q2 0.01) "rotations should match"

test "ZYX round trip through quaternion produces same rotation" := do
  let e := Euler.zyx 0.3 0.2 0.1
  let q := e.toQuat
  let e2 := Euler.fromQuat q .ZYX
  let q2 := e2.toQuat
  ensure (q.sameRotation q2 0.01) "rotations should match"

test "round trip through matrix produces same rotation" := do
  let e := Euler.xyz 0.3 0.2 0.1
  let m := e.toMat3
  let e2 := Euler.fromMat3 m .XYZ
  let q1 := e.toQuat
  let q2 := e2.toQuat
  ensure (q1.sameRotation q2 0.01) "rotations should match"

testSuite "Euler Utilities"

test "clamp limits angles" := do
  let e := Euler.xyz 1.0 2.0 3.0
  let clamped := e.clamp (Vec3.mk (-0.5) (-0.5) (-0.5)) (Vec3.mk 0.5 0.5 0.5)
  ensure (floatNear clamped.a1 0.5 0.0001) "a1 should be clamped to 0.5"
  ensure (floatNear clamped.a2 0.5 0.0001) "a2 should be clamped to 0.5"
  ensure (floatNear clamped.a3 0.5 0.0001) "a3 should be clamped to 0.5"

test "lerp at 0 returns first" := do
  let a := Euler.xyz 0.0 0.0 0.0
  let b := Euler.xyz 1.0 1.0 1.0
  let result := Euler.lerp a b 0.0
  ensure (floatNear result.a1 0.0 0.0001) "a1 should be 0"

test "lerp at 1 returns second" := do
  let a := Euler.xyz 0.0 0.0 0.0
  let b := Euler.xyz 1.0 1.0 1.0
  let result := Euler.lerp a b 1.0
  ensure (floatNear result.a1 1.0 0.0001) "a1 should be 1"

test "lerp at 0.5 returns midpoint" := do
  let a := Euler.xyz 0.0 0.0 0.0
  let b := Euler.xyz 1.0 1.0 1.0
  let result := Euler.lerp a b 0.5
  ensure (floatNear result.a1 0.5 0.0001) "a1 should be 0.5"

testSuite "Quat Euler Extensions"

test "fromEulerXYZ works" := do
  let q := Quat.fromEulerXYZ 0.1 0.2 0.3
  let e := Euler.xyz 0.1 0.2 0.3
  ensure (q.sameRotation e.toQuat) "should match Euler.toQuat"

test "fromEulerYXZ works" := do
  let q := Quat.fromEulerYXZ 0.1 0.2 0.3
  let e := Euler.yxz 0.1 0.2 0.3
  ensure (q.sameRotation e.toQuat) "should match Euler.toQuat"

test "toEulerXYZ extracts angles" := do
  let q := Quat.fromAxisAngle Vec3.unitY 0.5
  let angles := q.toEulerXYZ
  -- For a pure Y rotation, the Y component should be the angle
  ensure (angles.y.abs > 0.4) "y angle should be significant"



end LinalgTests.EulerTests
