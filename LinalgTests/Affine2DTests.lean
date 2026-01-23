/-
  Tests for Affine2D operations.
-/

import Linalg
import Crucible

namespace LinalgTests.Affine2DTests

open Crucible
open Linalg

testSuite "Affine2D Identity"

test "identity transforms point to itself" := do
  let p := Vec2.mk 3.0 4.0
  let result := Affine2D.identity * p
  ensure (floatNear result.x 3.0 0.0001) "x should be unchanged"
  ensure (floatNear result.y 4.0 0.0001) "y should be unchanged"

test "identity isIdentity returns true" := do
  ensure Affine2D.identity.isIdentity "identity should be identity"

test "identity determinant is 1" := do
  ensure (floatNear Affine2D.identity.determinant 1.0 0.0001) "det should be 1"

testSuite "Affine2D Translation"

test "translation moves point" := do
  let t := Affine2D.translationMat 10.0 20.0
  let p := Vec2.mk 1.0 2.0
  let result := t * p
  ensure (floatNear result.x 11.0 0.0001) "x should be 11"
  ensure (floatNear result.y 22.0 0.0001) "y should be 22"

test "translationV works with vector" := do
  let t := Affine2D.translationV (Vec2.mk 5.0 (-3.0))
  let p := Vec2.mk 0.0 0.0
  let result := t * p
  ensure (floatNear result.x 5.0 0.0001) "x should be 5"
  ensure (floatNear result.y (-3.0) 0.0001) "y should be -3"

test "translation component extraction" := do
  let t := Affine2D.translationMat 7.0 8.0
  let trans := t.translation
  ensure (floatNear trans.x 7.0 0.0001) "tx should be 7"
  ensure (floatNear trans.y 8.0 0.0001) "ty should be 8"

test "transformVector ignores translation" := do
  let t := Affine2D.translationMat 100.0 200.0
  let v := Vec2.mk 1.0 0.0
  let result := t.transformVector v
  ensure (floatNear result.x 1.0 0.0001) "x should be 1 (no translation)"
  ensure (floatNear result.y 0.0 0.0001) "y should be 0 (no translation)"

testSuite "Affine2D Rotation"

test "rotation by 90 degrees" := do
  let r := Affine2D.rotation (Float.pi / 2.0)
  let p := Vec2.mk 1.0 0.0
  let result := r * p
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y 1.0 0.0001) "y should be 1"

test "rotation by 180 degrees" := do
  let r := Affine2D.rotation Float.pi
  let p := Vec2.mk 1.0 0.0
  let result := r * p
  ensure (floatNear result.x (-1.0) 0.0001) "x should be -1"
  ensure (floatNear result.y 0.0 0.001) "y should be 0"

test "rotation by 45 degrees" := do
  let r := Affine2D.rotation (Float.pi / 4.0)
  let p := Vec2.mk 1.0 0.0
  let result := r * p
  let expected := Float.sqrt 2.0 / 2.0
  ensure (floatNear result.x expected 0.0001) "x should be sqrt(2)/2"
  ensure (floatNear result.y expected 0.0001) "y should be sqrt(2)/2"

test "extract rotation angle" := do
  let angle := Float.pi / 3.0
  let r := Affine2D.rotation angle
  let extracted := r.extractRotation
  ensure (floatNear extracted angle 0.0001) "extracted angle should match"

testSuite "Affine2D Scaling"

test "uniform scaling" := do
  let s := Affine2D.uniformScaling 2.0
  let p := Vec2.mk 3.0 4.0
  let result := s * p
  ensure (floatNear result.x 6.0 0.0001) "x should be 6"
  ensure (floatNear result.y 8.0 0.0001) "y should be 8"

test "non-uniform scaling" := do
  let s := Affine2D.scaling 2.0 3.0
  let p := Vec2.mk 1.0 1.0
  let result := s * p
  ensure (floatNear result.x 2.0 0.0001) "x should be 2"
  ensure (floatNear result.y 3.0 0.0001) "y should be 3"

test "extract scale factors" := do
  let s := Affine2D.scaling 2.0 3.0
  let extracted := s.extractScale
  ensure (floatNear extracted.x 2.0 0.0001) "sx should be 2"
  ensure (floatNear extracted.y 3.0 0.0001) "sy should be 3"

testSuite "Affine2D Composition"

test "compose translation with rotation" := do
  -- First translate, then rotate
  let t := Affine2D.translationMat 1.0 0.0
  let r := Affine2D.rotation (Float.pi / 2.0)
  let combined := r * t  -- rotate after translate
  let p := Vec2.mk 0.0 0.0
  let result := combined * p
  -- Point starts at origin, translated to (1,0), then rotated 90 to (0,1)
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y 1.0 0.0001) "y should be 1"

test "compose scaling with translation" := do
  let s := Affine2D.scaling 2.0 2.0
  let t := Affine2D.translationMat 10.0 10.0
  let combined := t * s  -- translate after scale
  let p := Vec2.mk 1.0 1.0
  let result := combined * p
  -- Scale (1,1) to (2,2), then translate to (12,12)
  ensure (floatNear result.x 12.0 0.0001) "x should be 12"
  ensure (floatNear result.y 12.0 0.0001) "y should be 12"

test "identity composed with transform equals transform" := do
  let t := Affine2D.translationMat 5.0 10.0
  let combined := Affine2D.identity * t
  let p := Vec2.mk 1.0 2.0
  let r1 := t * p
  let r2 := combined * p
  ensure (floatNear r1.x r2.x 0.0001) "x should match"
  ensure (floatNear r1.y r2.y 0.0001) "y should match"

testSuite "Affine2D Inverse"

test "inverse of identity is identity" := do
  match Affine2D.identity.inverse with
  | some inv => ensure inv.isIdentity "inverse of identity should be identity"
  | none => ensure false "identity should be invertible"

test "transform then inverse returns original" := do
  let t := Affine2D.translationMat 10.0 20.0
  match t.inverse with
  | some inv =>
    let p := Vec2.mk 5.0 7.0
    let transformed := t * p
    let restored := inv * transformed
    ensure (floatNear restored.x 5.0 0.0001) "x should be restored"
    ensure (floatNear restored.y 7.0 0.0001) "y should be restored"
  | none => ensure false "translation should be invertible"

test "rotation inverse works" := do
  let r := Affine2D.rotation (Float.pi / 4.0)
  match r.inverse with
  | some inv =>
    let p := Vec2.mk 1.0 0.0
    let transformed := r * p
    let restored := inv * transformed
    ensure (floatNear restored.x 1.0 0.0001) "x should be restored"
    ensure (floatNear restored.y 0.0 0.0001) "y should be restored"
  | none => ensure false "rotation should be invertible"

test "singular matrix returns None" := do
  -- Zero scale is singular
  let singular := Affine2D.scaling 0.0 0.0
  match singular.inverse with
  | some _ => ensure false "singular matrix should not have inverse"
  | none => ensure true "correctly returned None"

testSuite "Affine2D Shear"

test "shear along x" := do
  let sh := Affine2D.shear 1.0 0.0
  let p := Vec2.mk 0.0 1.0
  let result := sh * p
  ensure (floatNear result.x 1.0 0.0001) "x should be 1 (sheared by y)"
  ensure (floatNear result.y 1.0 0.0001) "y should stay 1"

test "shear along y" := do
  let sh := Affine2D.shear 0.0 1.0
  let p := Vec2.mk 1.0 0.0
  let result := sh * p
  ensure (floatNear result.x 1.0 0.0001) "x should stay 1"
  ensure (floatNear result.y 1.0 0.0001) "y should be 1 (sheared by x)"

testSuite "Affine2D Pivot Operations"

test "rotation around pivot" := do
  let pivot := Vec2.mk 1.0 0.0
  let r := Affine2D.rotationAround pivot (Float.pi / 2.0)
  -- Rotate point (2,0) around pivot (1,0) by 90 degrees
  -- Relative to pivot: (1,0), after rotation: (0,1), absolute: (1,1)
  let p := Vec2.mk 2.0 0.0
  let result := r * p
  ensure (floatNear result.x 1.0 0.0001) "x should be 1"
  ensure (floatNear result.y 1.0 0.0001) "y should be 1"

test "scaling around pivot" := do
  let pivot := Vec2.mk 1.0 1.0
  let s := Affine2D.scalingAround pivot 2.0 2.0
  -- Scale point (2,2) around pivot (1,1) by 2x
  -- Relative to pivot: (1,1), after scale: (2,2), absolute: (3,3)
  let p := Vec2.mk 2.0 2.0
  let result := s * p
  ensure (floatNear result.x 3.0 0.0001) "x should be 3"
  ensure (floatNear result.y 3.0 0.0001) "y should be 3"

testSuite "Affine2D Conversion"

test "toMat2 extracts linear part" := do
  let a := Affine2D.rotation (Float.pi / 4.0)
  let m := a.toMat2
  let c := Float.sqrt 2.0 / 2.0
  ensure (floatNear (m.get 0 0) c 0.0001) "m00 should be cos"
  ensure (floatNear (m.get 1 0) c 0.0001) "m10 should be sin"

test "basisX and basisY extraction" := do
  let a := Affine2D.scaling 2.0 3.0
  let bx := a.basisX
  let basisY := a.basisY
  ensure (floatNear bx.x 2.0 0.0001) "basisX.x should be 2"
  ensure (floatNear bx.y 0.0 0.0001) "basisX.y should be 0"
  ensure (floatNear basisY.x 0.0 0.0001) "basisY.x should be 0"
  ensure (floatNear basisY.y 3.0 0.0001) "basisY.y should be 3"



end LinalgTests.Affine2DTests
