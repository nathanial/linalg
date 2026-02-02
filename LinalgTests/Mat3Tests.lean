/-
  Tests for Mat3 operations.
-/

import Linalg
import Crucible

namespace LinalgTests.Mat3Tests

open Crucible
open Linalg

testSuite "Mat3 Basic"

test "identity matrix works correctly" := do
  let m := Mat3.identity
  m.get 0 0 ≡ 1.0
  m.get 1 1 ≡ 1.0
  m.get 2 2 ≡ 1.0
  m.get 0 1 ≡ 0.0

test "identity times vector returns vector" := do
  let v := Vec3.mk 1.0 2.0 3.0
  let result := Mat3.identity * v
  result.x ≡ 1.0
  result.y ≡ 2.0
  result.z ≡ 3.0

testSuite "Mat3 Operations"

test "determinant of identity is 1" := do
  Mat3.identity.determinant ≡ 1.0

test "scaling matrix determinant" := do
  let m := Mat3.scaling 2.0 3.0 4.0
  ensure (floatNear m.determinant 24.0 0.0001) "determinant should be 24"

test "rotation X preserves X axis" := do
  let m := Mat3.rotationX (Float.pi / 4.0)
  let v := Vec3.unitX
  let result := m * v
  ensure (floatNear result.x 1.0 0.0001) "x should be 1"
  ensure (floatNear result.y 0.0 0.0001) "y should be 0"
  ensure (floatNear result.z 0.0 0.0001) "z should be 0"

test "rotation Y preserves Y axis" := do
  let m := Mat3.rotationY (Float.pi / 4.0)
  let v := Vec3.unitY
  let result := m * v
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y 1.0 0.0001) "y should be 1"
  ensure (floatNear result.z 0.0 0.0001) "z should be 0"

test "rotation Z preserves Z axis" := do
  let m := Mat3.rotationZ (Float.pi / 4.0)
  let v := Vec3.unitZ
  let result := m * v
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.y 0.0 0.0001) "y should be 0"
  ensure (floatNear result.z 1.0 0.0001) "z should be 1"

testSuite "Mat3 Solvers"

test "LU/QR/Cholesky solve" := do
  let a := Mat3.fromColumns
    (Vec3.mk 4.0 1.0 2.0)
    (Vec3.mk 1.0 3.0 0.0)
    (Vec3.mk 2.0 0.0 5.0)
  let xTrue := Vec3.mk 1.0 2.0 3.0
  let b := a * xTrue
  match a.solve b with
  | some x =>
      ensure (floatNear x.x xTrue.x 0.0001) "lu x"
      ensure (floatNear x.y xTrue.y 0.0001) "lu y"
      ensure (floatNear x.z xTrue.z 0.0001) "lu z"
  | none => ensure false "lu solve failed"
  match a.qrDecompose with
  | some qr =>
      match Mat3.solveQR qr b with
      | some x =>
          ensure (floatNear x.x xTrue.x 0.0001) "qr x"
          ensure (floatNear x.y xTrue.y 0.0001) "qr y"
          ensure (floatNear x.z xTrue.z 0.0001) "qr z"
      | none => ensure false "qr solve failed"
  | none => ensure false "qr decompose failed"
  match a.choleskyDecompose with
  | some chol =>
      match Mat3.solveCholesky chol b with
      | some x =>
          ensure (floatNear x.x xTrue.x 0.0001) "chol x"
          ensure (floatNear x.y xTrue.y 0.0001) "chol y"
          ensure (floatNear x.z xTrue.z 0.0001) "chol z"
      | none => ensure false "chol solve failed"
  | none => ensure false "chol decompose failed"



end LinalgTests.Mat3Tests
