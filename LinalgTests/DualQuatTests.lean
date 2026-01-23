/-
  Tests for Dual Quaternion type.
-/

import Linalg
import Crucible

namespace LinalgTests.DualQuatTests

open Crucible
open Linalg

testSuite "DualQuat Construction"

test "identity dual quaternion" := do
  let dq := DualQuat.identity
  ensure (dq.real.sameRotation Quat.identity) "real should be identity"
  let t := dq.getTranslation
  ensure (floatNear t.x 0.0 0.0001) "translation x should be 0"
  ensure (floatNear t.y 0.0 0.0001) "translation y should be 0"
  ensure (floatNear t.z 0.0 0.0001) "translation z should be 0"

test "fromRotation creates pure rotation" := do
  let rot := Quat.fromAxisAngle Vec3.unitY Float.halfPi
  let dq := DualQuat.fromRotation rot
  ensure (dq.real.sameRotation rot) "rotation should match"
  let t := dq.getTranslation
  ensure (floatNear t.length 0.0 0.0001) "translation should be zero"

test "fromTranslation creates pure translation" := do
  let trans := Vec3.mk 5.0 10.0 15.0
  let dq := DualQuat.fromTranslation trans
  ensure (dq.real.sameRotation Quat.identity) "rotation should be identity"
  let t := dq.getTranslation
  ensure (floatNear t.x 5.0 0.01) "translation x"
  ensure (floatNear t.y 10.0 0.01) "translation y"
  ensure (floatNear t.z 15.0 0.01) "translation z"

test "fromRotationTranslation combines both" := do
  let rot := Quat.fromAxisAngle Vec3.unitZ Float.halfPi
  let trans := Vec3.mk 1.0 2.0 3.0
  let dq := DualQuat.fromRotationTranslation rot trans
  ensure (dq.real.sameRotation rot) "rotation should match"
  let t := dq.getTranslation
  ensure (floatNear t.x 1.0 0.01) "translation x"
  ensure (floatNear t.y 2.0 0.01) "translation y"
  ensure (floatNear t.z 3.0 0.01) "translation z"

testSuite "DualQuat Transform"

test "identity transforms point unchanged" := do
  let dq := DualQuat.identity
  let p := Vec3.mk 1.0 2.0 3.0
  let result := dq.transformPoint p
  ensure (floatNear result.x 1.0 0.0001) "x"
  ensure (floatNear result.y 2.0 0.0001) "y"
  ensure (floatNear result.z 3.0 0.0001) "z"

test "pure translation moves point" := do
  let dq := DualQuat.fromTranslation (Vec3.mk 10.0 0.0 0.0)
  let p := Vec3.mk 1.0 0.0 0.0
  let result := dq.transformPoint p
  ensure (floatNear result.x 11.0 0.01) "x should be 11"

test "pure rotation rotates point" := do
  let dq := DualQuat.fromRotation (Quat.fromAxisAngle Vec3.unitZ Float.halfPi)
  let p := Vec3.unitX
  let result := dq.transformPoint p
  ensure (floatNear result.x 0.0 0.01) "x should be 0"
  ensure (floatNear result.y 1.0 0.01) "y should be 1"

test "transformVector ignores translation" := do
  let dq := DualQuat.fromRotationTranslation Quat.identity (Vec3.mk 100.0 100.0 100.0)
  let v := Vec3.mk 1.0 0.0 0.0
  let result := dq.transformVector v
  ensure (floatNear result.x 1.0 0.0001) "x should be 1"
  ensure (floatNear result.y 0.0 0.0001) "y should be 0"

testSuite "DualQuat Composition"

test "multiply identity gives same" := do
  let dq := DualQuat.fromRotationTranslation
    (Quat.fromAxisAngle Vec3.unitY 0.5)
    (Vec3.mk 1.0 2.0 3.0)
  let result := dq * DualQuat.identity
  ensure (dq.sameTransform result) "should be same transform"

test "composition applies transforms in order" := do
  let translate := DualQuat.fromTranslation (Vec3.mk 5.0 0.0 0.0)
  let rotate := DualQuat.fromRotation (Quat.fromAxisAngle Vec3.unitZ Float.halfPi)
  -- rotate * translate means: first translate, then rotate
  let combined := rotate * translate
  let p := Vec3.zero
  let result := combined.transformPoint p
  -- Translate (0,0,0) by (5,0,0) -> (5,0,0)
  -- Rotate (5,0,0) by 90° around Z -> (0,5,0)
  ensure (floatNear result.x 0.0 0.01) "x should be 0"
  ensure (floatNear result.y 5.0 0.01) "y should be 5"

testSuite "DualQuat Inverse"

test "inverse of identity is identity" := do
  let dq := DualQuat.identity
  let inv := dq.inverse
  ensure (inv.sameTransform DualQuat.identity) "inverse of identity is identity"

test "inverse undoes transform" := do
  let dq := DualQuat.fromRotationTranslation
    (Quat.fromAxisAngle Vec3.unitY 0.7)
    (Vec3.mk 3.0 4.0 5.0)
  let inv := dq.inverse
  let p := Vec3.mk 1.0 2.0 3.0
  let transformed := dq.transformPoint p
  let result := inv.transformPoint transformed
  ensure (floatNear result.x p.x 0.01) "x should return"
  ensure (floatNear result.y p.y 0.01) "y should return"
  ensure (floatNear result.z p.z 0.01) "z should return"

test "dq * inverse = identity" := do
  let dq := DualQuat.fromRotationTranslation
    (Quat.fromAxisAngle Vec3.unitX 0.5)
    (Vec3.mk 1.0 2.0 3.0)
  let inv := dq.inverse
  let result := dq * inv
  ensure (result.sameTransform DualQuat.identity 0.01) "should be identity"

testSuite "DualQuat Interpolation"

test "lerp at 0 returns first" := do
  let a := DualQuat.fromTranslation Vec3.zero
  let b := DualQuat.fromTranslation (Vec3.mk 10.0 0.0 0.0)
  let result := DualQuat.lerp a b 0.0
  let t := result.getTranslation
  ensure (floatNear t.x 0.0 0.1) "should be at first position"

test "lerp at 1 returns second" := do
  let a := DualQuat.fromTranslation Vec3.zero
  let b := DualQuat.fromTranslation (Vec3.mk 10.0 0.0 0.0)
  let result := DualQuat.lerp a b 1.0
  let t := result.getTranslation
  ensure (floatNear t.x 10.0 0.1) "should be at second position"

test "lerp at 0.5 returns midpoint" := do
  let a := DualQuat.fromTranslation Vec3.zero
  let b := DualQuat.fromTranslation (Vec3.mk 10.0 0.0 0.0)
  let result := DualQuat.lerp a b 0.5
  let t := result.getTranslation
  ensure (floatNear t.x 5.0 0.5) "should be at midpoint"

testSuite "DualQuat Blending"

test "blend with single element returns same" := do
  let dq := DualQuat.fromRotationTranslation
    (Quat.fromAxisAngle Vec3.unitY 0.5)
    (Vec3.mk 1.0 2.0 3.0)
  let result := DualQuat.blend #[dq] #[1.0]
  ensure (result.sameTransform dq 0.1) "should be same transform"

test "blend with equal weights averages" := do
  let a := DualQuat.fromTranslation Vec3.zero
  let b := DualQuat.fromTranslation (Vec3.mk 10.0 0.0 0.0)
  let result := DualQuat.blend #[a, b] #[0.5, 0.5]
  let t := result.getTranslation
  ensure (floatNear t.x 5.0 0.5) "should be at midpoint"

testSuite "DualQuat Matrix Conversion"

test "toMat4 round trip" := do
  let dq := DualQuat.fromRotationTranslation
    (Quat.fromAxisAngle Vec3.unitY 0.5)
    (Vec3.mk 1.0 2.0 3.0)
  let m := dq.toMat4
  let dq2 := DualQuat.fromMat4 m
  ensure (dq.sameTransform dq2 0.01) "should round trip"

test "toMat4 transforms point correctly" := do
  let dq := DualQuat.fromRotationTranslation
    (Quat.fromAxisAngle Vec3.unitZ Float.halfPi)
    (Vec3.mk 5.0 0.0 0.0)
  let p := Vec3.unitX
  let dqResult := dq.transformPoint p
  let m := dq.toMat4
  let mResult := m.transformPoint p
  ensure (floatNear dqResult.x mResult.x 0.01) "x should match"
  ensure (floatNear dqResult.y mResult.y 0.01) "y should match"
  ensure (floatNear dqResult.z mResult.z 0.01) "z should match"



end LinalgTests.DualQuatTests
