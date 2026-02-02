/-
  Analysis helpers: PCA, covariance, and fitting utilities.
-/

import Linalg.Vec2
import Linalg.Vec3
import Linalg.Mat2
import Linalg.Mat3

namespace Linalg

/-- PCA result for 2D points. -/
structure PCA2 where
  mean : Vec2
  values : Vec2
  vectors : Mat2
deriving Repr, Inhabited

/-- PCA result for 3D points. -/
structure PCA3 where
  mean : Vec3
  values : Vec3
  vectors : Mat3
deriving Repr, Inhabited

/-- Mean of a set of 2D points. -/
def mean2 (points : Array Vec2) : Vec2 :=
  if points.isEmpty then
    Vec2.zero
  else
    let sum := points.foldl (fun acc p => acc + p) Vec2.zero
    sum / Float.ofNat points.size

/-- Mean of a set of 3D points. -/
def mean3 (points : Array Vec3) : Vec3 :=
  if points.isEmpty then
    Vec3.zero
  else
    let sum := points.foldl (fun acc p => acc + p) Vec3.zero
    sum / Float.ofNat points.size

private def covariance2WithMean (points : Array Vec2) : Vec2 × Mat2 := Id.run do
  let mean := mean2 points
  let mut s00 := 0.0
  let mut s01 := 0.0
  let mut s11 := 0.0
  for p in points do
    let dx := p.x - mean.x
    let dy := p.y - mean.y
    s00 := s00 + dx * dx
    s01 := s01 + dx * dy
    s11 := s11 + dy * dy
  let denom :=
    if points.size > 1 then Float.ofNat (points.size - 1) else 1.0
  let inv := 1.0 / denom
  let cov := Mat2.fromRows
    (Vec2.mk (s00 * inv) (s01 * inv))
    (Vec2.mk (s01 * inv) (s11 * inv))
  return (mean, cov)

private def covariance3WithMean (points : Array Vec3) : Vec3 × Mat3 := Id.run do
  let mean := mean3 points
  let mut s00 := 0.0
  let mut s01 := 0.0
  let mut s02 := 0.0
  let mut s11 := 0.0
  let mut s12 := 0.0
  let mut s22 := 0.0
  for p in points do
    let dx := p.x - mean.x
    let dy := p.y - mean.y
    let dz := p.z - mean.z
    s00 := s00 + dx * dx
    s01 := s01 + dx * dy
    s02 := s02 + dx * dz
    s11 := s11 + dy * dy
    s12 := s12 + dy * dz
    s22 := s22 + dz * dz
  let denom :=
    if points.size > 1 then Float.ofNat (points.size - 1) else 1.0
  let inv := 1.0 / denom
  let cov :=
    Mat3.fromColumns
      (Vec3.mk (s00 * inv) (s01 * inv) (s02 * inv))
      (Vec3.mk (s01 * inv) (s11 * inv) (s12 * inv))
      (Vec3.mk (s02 * inv) (s12 * inv) (s22 * inv))
  return (mean, cov)

/-- Covariance matrix of 2D points (sample covariance when size > 1). -/
def covariance2 (points : Array Vec2) : Mat2 :=
  (covariance2WithMean points).2

/-- Covariance matrix of 3D points (sample covariance when size > 1). -/
def covariance3 (points : Array Vec3) : Mat3 :=
  (covariance3WithMean points).2

/-- Principal component analysis for 2D points. -/
def pca2 (points : Array Vec2) : Option PCA2 :=
  if points.isEmpty then
    none
  else
    let (mean, cov) := covariance2WithMean points
    let eigen := Mat2.eigenDecomposeSymmetric cov
    some { mean := mean, values := eigen.values, vectors := eigen.vectors }

/-- Principal component analysis for 3D points. -/
def pca3 (points : Array Vec3) : Option PCA3 :=
  if points.isEmpty then
    none
  else
    let (mean, cov) := covariance3WithMean points
    let eigen := Mat3.eigenDecomposeSymmetric cov
    some { mean := mean, values := eigen.values, vectors := eigen.vectors }

/-- Project a 2D point into PCA coordinates. -/
def pca2Project (pca : PCA2) (p : Vec2) : Vec2 :=
  let d := p - pca.mean
  let v0 := pca.vectors.column 0
  let v1 := pca.vectors.column 1
  Vec2.mk (d.dot v0) (d.dot v1)

/-- Project a 3D point into PCA coordinates. -/
def pca3Project (pca : PCA3) (p : Vec3) : Vec3 :=
  let d := p - pca.mean
  let v0 := pca.vectors.column 0
  let v1 := pca.vectors.column 1
  let v2 := pca.vectors.column 2
  Vec3.mk (d.dot v0) (d.dot v1) (d.dot v2)

/-- Best-fit line for 2D points (point on line, direction). -/
def bestFitLine2 (points : Array Vec2) : Option (Vec2 × Vec2) :=
  match pca2 points with
  | none => none
  | some pca => some (pca.mean, pca.vectors.column 0)

/-- Best-fit line for 3D points (point on line, direction). -/
def bestFitLine3 (points : Array Vec3) : Option (Vec3 × Vec3) :=
  match pca3 points with
  | none => none
  | some pca => some (pca.mean, pca.vectors.column 0)

/-- Best-fit plane for 3D points (point on plane, normal). -/
def bestFitPlane3 (points : Array Vec3) : Option (Vec3 × Vec3) :=
  match pca3 points with
  | none => none
  | some pca => some (pca.mean, pca.vectors.column 2)

end Linalg
