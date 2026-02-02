/-
  Tests for analysis helpers (eigen, SVD, PCA).
-/

import Linalg
import Crucible

namespace LinalgTests.AnalysisTests

open Crucible
open Linalg

testSuite "Analysis Decompositions"

test "Mat2 symmetric eigen reconstruction" := do
  let m := Mat2.fromRows (Vec2.mk 3.0 1.0) (Vec2.mk 1.0 2.0)
  let e := Mat2.eigenDecomposeSymmetric m
  let diag := Mat2.scaling e.values.x e.values.y
  let recon := e.vectors * diag * e.vectors.transpose
  ensure (recon.approxEq m 0.001) "eigen reconstruction"

test "Mat2 SVD reconstruction" := do
  let m := Mat2.fromRows (Vec2.mk 3.0 1.0) (Vec2.mk 0.0 2.0)
  let s := Mat2.svd m
  let diag := Mat2.scaling s.S.x s.S.y
  let recon := s.U * diag * s.V.transpose
  ensure (recon.approxEq m 0.001) "svd reconstruction"

test "Mat3 symmetric eigen reconstruction" := do
  let r := Mat3.rotationZ (Float.pi / 4.0)
  let d := Mat3.scaling 5.0 3.0 2.0
  let m := r * d * r.transpose
  let e := Mat3.eigenDecomposeSymmetric m
  let diag := Mat3.scaling e.values.x e.values.y e.values.z
  let recon := e.vectors * diag * e.vectors.transpose
  ensure (recon.approxEq m 0.001) "eigen reconstruction"

test "Mat3 SVD reconstruction" := do
  let m := Mat3.scaling 4.0 2.0 1.0
  let s := Mat3.svd m
  let diag := Mat3.scaling s.S.x s.S.y s.S.z
  let recon := s.U * diag * s.V.transpose
  ensure (recon.approxEq m 0.001) "svd reconstruction"

testSuite "PCA"

test "PCA2 principal axis aligns with line" := do
  let points := #[
    Vec2.mk 0.0 0.0,
    Vec2.mk 1.0 2.0,
    Vec2.mk 2.0 4.0,
    Vec2.mk 3.0 6.0
  ]
  match pca2 points with
  | none => ensure false "pca2 returned none"
  | some pca =>
      ensure (floatNear pca.mean.x 1.5 0.0001) "mean x"
      ensure (floatNear pca.mean.y 3.0 0.0001) "mean y"
      let dir := pca.vectors.column 0
      let expected := (Vec2.mk 1.0 2.0).normalize
      ensure (Float.abs' (dir.dot expected) > 0.999) "principal axis"

test "PCA3 normal aligns with plane" := do
  let points := #[
    Vec3.mk 0.0 0.0 0.0,
    Vec3.mk 1.0 0.0 0.0,
    Vec3.mk 0.0 1.0 0.0,
    Vec3.mk 1.0 1.0 0.0,
    Vec3.mk 2.0 0.0 0.0
  ]
  match pca3 points with
  | none => ensure false "pca3 returned none"
  | some pca =>
      let normal := pca.vectors.column 2
      ensure (Float.abs' (normal.dot Vec3.unitZ) > 0.999) "plane normal"

end LinalgTests.AnalysisTests
