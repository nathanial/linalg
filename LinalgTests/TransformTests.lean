/-
  Tests for Transform type.
-/

import Linalg
import Crucible

namespace LinalgTests.TransformTests

open Crucible
open Linalg

testSuite "Transform Basic"

test "identity transform" := do
  let t := Transform.identity
  ensure (floatNear t.position.x 0.0 0.0001) "position x should be 0"
  ensure (floatNear t.position.y 0.0 0.0001) "position y should be 0"
  ensure (floatNear t.position.z 0.0 0.0001) "position z should be 0"
  ensure (floatNear t.scale.x 1.0 0.0001) "scale x should be 1"
  ensure (t.rotation.sameRotation Quat.identity) "rotation should be identity"

test "fromPosition creates correct transform" := do
  let t := Transform.fromPosition (Vec3.mk 5.0 10.0 15.0)
  ensure (floatNear t.position.x 5.0 0.0001) "position x"
  ensure (floatNear t.position.y 10.0 0.0001) "position y"
  ensure (floatNear t.position.z 15.0 0.0001) "position z"

test "fromUniformScale creates correct transform" := do
  let t := Transform.fromUniformScale 2.0
  ensure (floatNear t.scale.x 2.0 0.0001) "scale x"
  ensure (floatNear t.scale.y 2.0 0.0001) "scale y"
  ensure (floatNear t.scale.z 2.0 0.0001) "scale z"

testSuite "Transform Point"

test "identity transform leaves point unchanged" := do
  let t := Transform.identity
  let p := Vec3.mk 1.0 2.0 3.0
  let result := t.transformPoint p
  ensure (floatNear result.x 1.0 0.0001) "x"
  ensure (floatNear result.y 2.0 0.0001) "y"
  ensure (floatNear result.z 3.0 0.0001) "z"

test "translation moves point" := do
  let t := Transform.fromPosition (Vec3.mk 10.0 0.0 0.0)
  let p := Vec3.zero
  let result := t.transformPoint p
  ensure (floatNear result.x 10.0 0.0001) "x should be 10"

test "scale scales point" := do
  let t := Transform.fromUniformScale 2.0
  let p := Vec3.mk 1.0 1.0 1.0
  let result := t.transformPoint p
  ensure (floatNear result.x 2.0 0.0001) "x should be 2"
  ensure (floatNear result.y 2.0 0.0001) "y should be 2"
  ensure (floatNear result.z 2.0 0.0001) "z should be 2"

test "rotation rotates point" := do
  let t := Transform.fromRotation (Quat.fromAxisAngle Vec3.unitY Float.halfPi)
  let p := Vec3.unitX
  let result := t.transformPoint p
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.z (-1.0) 0.0001) "z should be -1"

testSuite "Transform Inverse"

test "inverse of identity is identity" := do
  let t := Transform.identity
  let inv := t.inverse
  ensure (inv.approxEq Transform.identity) "inverse of identity should be identity"

test "inverse undoes transform" := do
  let t := Transform.mk' (Vec3.mk 5.0 10.0 15.0) Quat.identity (Vec3.mk 2.0 2.0 2.0)
  let inv := t.inverse
  let p := Vec3.mk 1.0 2.0 3.0
  let transformed := t.transformPoint p
  let result := inv.transformPoint transformed
  ensure (floatNear result.x p.x 0.001) "x should return to original"
  ensure (floatNear result.y p.y 0.001) "y should return to original"
  ensure (floatNear result.z p.z 0.001) "z should return to original"

testSuite "Transform Composition"

test "composing with identity gives same transform" := do
  let t := Transform.fromPosition (Vec3.mk 5.0 0.0 0.0)
  let result := t * Transform.identity
  ensure (result.approxEq t) "composing with identity should give same"

test "translation composition" := do
  let t1 := Transform.fromPosition (Vec3.mk 5.0 0.0 0.0)
  let t2 := Transform.fromPosition (Vec3.mk 0.0 10.0 0.0)
  let result := t1 * t2
  ensure (floatNear result.position.x 5.0 0.0001) "x should be 5"
  ensure (floatNear result.position.y 10.0 0.0001) "y should be 10"

testSuite "Transform Lerp"

test "lerp at 0 returns first transform" := do
  let a := Transform.fromPosition Vec3.zero
  let b := Transform.fromPosition (Vec3.mk 10.0 0.0 0.0)
  let result := Transform.lerp a b 0.0
  ensure (floatNear result.position.x 0.0 0.0001) "should be at first position"

test "lerp at 1 returns second transform" := do
  let a := Transform.fromPosition Vec3.zero
  let b := Transform.fromPosition (Vec3.mk 10.0 0.0 0.0)
  let result := Transform.lerp a b 1.0
  ensure (floatNear result.position.x 10.0 0.0001) "should be at second position"

test "lerp at 0.5 returns midpoint" := do
  let a := Transform.fromPosition Vec3.zero
  let b := Transform.fromPosition (Vec3.mk 10.0 0.0 0.0)
  let result := Transform.lerp a b 0.5
  ensure (floatNear result.position.x 5.0 0.0001) "should be at midpoint"

#generate_tests

end LinalgTests.TransformTests
