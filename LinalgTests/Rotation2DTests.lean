/-
  Tests for Rotation2D operations.
-/

import Linalg
import Crucible

namespace LinalgTests.Rotation2DTests

open Crucible
open Linalg

testSuite "Rotation2D Identity"

test "identity rotation is zero" := do
  let r := Rotation2D.identity
  ensure (floatNear r.angle 0.0 0.0001) "identity angle should be 0"

test "identity leaves vector unchanged" := do
  let v := Vec2.mk 3.0 4.0
  let result := Rotation2D.identity * v
  ensure (floatNear result.x 3.0 0.0001) "x should be unchanged"
  ensure (floatNear result.y 4.0 0.0001) "y should be unchanged"

testSuite "Rotation2D Construction"

test "fromDegrees and degrees round-trip" := do
  let r := Rotation2D.fromDegrees 45.0
  ensure (floatNear r.degrees 45.0 0.0001) "should be 45 degrees"

test "fromRadians creates correct rotation" := do
  let r := Rotation2D.fromRadians (Float.pi / 2.0)
  ensure (floatNear r.angle (Float.pi / 2.0) 0.0001) "angle should be pi/2"

test "fromDirection creates rotation pointing in direction" := do
  let r := Rotation2D.fromDirection (Vec2.mk 0.0 1.0)
  ensure (floatNear r.angle (Float.pi / 2.0) 0.0001) "should be 90 degrees"

test "fromTo creates rotation between directions" := do
  let fromDir := Vec2.mk 1.0 0.0
  let toDir := Vec2.mk 0.0 1.0
  let r := Rotation2D.fromTo fromDir toDir
  ensure (floatNear r.angle (Float.pi / 2.0) 0.0001) "should be 90 degrees"

test "direction gives unit vector" := do
  let r := Rotation2D.fromDegrees 90.0
  let d := r.direction
  ensure (floatNear d.x 0.0 0.0001) "x should be 0"
  ensure (floatNear d.y 1.0 0.0001) "y should be 1"

testSuite "Rotation2D Operations"

test "rotate 90 degrees" := do
  let r := Rotation2D.fromDegrees 90.0
  let v := Vec2.mk 1.0 0.0
  let result := r * v
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y 1.0 0.0001) "y should be 1"

test "rotate 180 degrees" := do
  let r := Rotation2D.fromDegrees 180.0
  let v := Vec2.mk 1.0 0.0
  let result := r * v
  ensure (floatNear result.x (-1.0) 0.0001) "x should be -1"
  ensure (floatNear result.y 0.0 0.001) "y should be 0"

test "rotate 270 degrees" := do
  let r := Rotation2D.fromDegrees 270.0
  let v := Vec2.mk 1.0 0.0
  let result := r * v
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y (-1.0) 0.0001) "y should be -1"

test "rotate preserves length" := do
  let r := Rotation2D.fromDegrees 37.0
  let v := Vec2.mk 3.0 4.0
  let result := r * v
  ensure (floatNear result.length v.length 0.0001) "length should be preserved"

testSuite "Rotation2D Composition"

test "inverse rotation" := do
  let r := Rotation2D.fromDegrees 45.0
  let inv := r.inverse
  let composed := r * inv
  ensure (floatNear composed.angle 0.0 0.0001) "r * r^-1 should be identity"

test "compose adds angles" := do
  let r1 := Rotation2D.fromDegrees 30.0
  let r2 := Rotation2D.fromDegrees 60.0
  let composed := r1 * r2
  ensure (floatNear composed.degrees 90.0 0.0001) "30 + 60 = 90"

test "rotation then inverse returns original" := do
  let r := Rotation2D.fromDegrees 45.0
  let v := Vec2.mk 1.0 0.0
  let rotated := r * v
  let restored := r.inverse * rotated
  ensure (floatNear restored.x 1.0 0.0001) "x should be restored"
  ensure (floatNear restored.y 0.0 0.0001) "y should be restored"

testSuite "Rotation2D Interpolation"

test "lerp at 0 returns first" := do
  let r1 := Rotation2D.fromDegrees 0.0
  let r2 := Rotation2D.fromDegrees 90.0
  let result := Rotation2D.lerp r1 r2 0.0
  ensure (floatNear result.degrees 0.0 0.0001) "should be 0 degrees"

test "lerp at 1 returns second" := do
  let r1 := Rotation2D.fromDegrees 0.0
  let r2 := Rotation2D.fromDegrees 90.0
  let result := Rotation2D.lerp r1 r2 1.0
  ensure (floatNear result.degrees 90.0 0.0001) "should be 90 degrees"

test "lerp at 0.5 returns midpoint" := do
  let r1 := Rotation2D.fromDegrees 0.0
  let r2 := Rotation2D.fromDegrees 90.0
  let result := Rotation2D.lerp r1 r2 0.5
  ensure (floatNear result.degrees 45.0 0.0001) "should be 45 degrees"

testSuite "Rotation2D Conversion"

test "toMat2 produces correct rotation matrix" := do
  let r := Rotation2D.fromDegrees 90.0
  let m := r.toMat2
  let v := Vec2.mk 1.0 0.0
  let result := m * v
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y 1.0 0.0001) "y should be 1"

test "toAffine2D produces correct transform" := do
  let r := Rotation2D.fromDegrees 90.0
  let a := r.toAffine2D
  let v := Vec2.mk 1.0 0.0
  let result := a * v
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y 1.0 0.0001) "y should be 1"

testSuite "Rotation2D Normalization"

test "normalize brings angle to [-pi, pi]" := do
  let r := Rotation2D.fromRadians (3.0 * Float.pi)
  let normalized := r.normalize
  ensure (normalized.angle >= -Float.pi) "should be >= -pi"
  ensure (normalized.angle <= Float.pi) "should be <= pi"

test "cos and sin accessors" := do
  let r := Rotation2D.fromDegrees 60.0
  ensure (floatNear r.cos 0.5 0.0001) "cos(60) should be 0.5"
  ensure (floatNear r.sin (Float.sqrt 3.0 / 2.0) 0.0001) "sin(60) should be sqrt(3)/2"



end LinalgTests.Rotation2DTests
