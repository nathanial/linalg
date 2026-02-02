/-
  2x2 Matrix type and operations.
  Column-major storage for consistency with GPU conventions.
-/

import Linalg.Core
import Linalg.Vec2

namespace Linalg

/-- 2x2 matrix in column-major order: [m00, m10, m01, m11] -/
structure Mat2 where
  data : Array Float
deriving Inhabited, Repr

namespace Mat2

/-- Identity matrix. -/
def identity : Mat2 := { data := #[1, 0, 0, 1] }

/-- Zero matrix. -/
def zero : Mat2 := { data := #[0, 0, 0, 0] }

/-- Create matrix from column vectors. -/
def fromColumns (c0 c1 : Vec2) : Mat2 :=
  { data := #[c0.x, c0.y, c1.x, c1.y] }

/-- Create matrix from row vectors. -/
def fromRows (r0 r1 : Vec2) : Mat2 :=
  { data := #[r0.x, r1.x, r0.y, r1.y] }

/-- Get element at (row, col). -/
@[inline]
def get (m : Mat2) (row col : Nat) : Float :=
  m.data.getD (col * 2 + row) 0.0

/-- Set element at (row, col). -/
def set (m : Mat2) (row col : Nat) (v : Float) : Mat2 :=
  { data := m.data.set! (col * 2 + row) v }

/-- Get column as Vec2. -/
def column (m : Mat2) (c : Nat) : Vec2 :=
  ⟨m.get 0 c, m.get 1 c⟩

/-- Get row as Vec2. -/
def row (m : Mat2) (r : Nat) : Vec2 :=
  ⟨m.get r 0, m.get r 1⟩

/-- Component-wise addition. -/
def add (a b : Mat2) : Mat2 :=
  let result := Id.run do
    let mut arr : Array Float := #[]
    for i in [:a.data.size] do
      arr := arr.push (a.data.getD i 0.0 + b.data.getD i 0.0)
    return arr
  { data := result }

/-- Component-wise subtraction. -/
def sub (a b : Mat2) : Mat2 :=
  let result := Id.run do
    let mut arr : Array Float := #[]
    for i in [:a.data.size] do
      arr := arr.push (a.data.getD i 0.0 - b.data.getD i 0.0)
    return arr
  { data := result }

/-- Scale all elements by a scalar. -/
def scale (m : Mat2) (s : Float) : Mat2 :=
  { data := m.data.map (· * s) }

/-- Matrix multiplication. -/
def multiply (a b : Mat2) : Mat2 :=
  let m00 := a.get 0 0 * b.get 0 0 + a.get 0 1 * b.get 1 0
  let m10 := a.get 1 0 * b.get 0 0 + a.get 1 1 * b.get 1 0
  let m01 := a.get 0 0 * b.get 0 1 + a.get 0 1 * b.get 1 1
  let m11 := a.get 1 0 * b.get 0 1 + a.get 1 1 * b.get 1 1
  { data := #[m00, m10, m01, m11] }

/-- Transpose the matrix. -/
def transpose (m : Mat2) : Mat2 :=
  { data := #[m.get 0 0, m.get 0 1, m.get 1 0, m.get 1 1] }

/-- Compute the determinant. -/
def determinant (m : Mat2) : Float :=
  m.get 0 0 * m.get 1 1 - m.get 0 1 * m.get 1 0

/-- Compute the inverse. Returns none if matrix is singular. -/
def inverse (m : Mat2) : Option Mat2 :=
  let det := m.determinant
  if Float.abs' det < Float.epsilon then none
  else
    let invDet := 1.0 / det
    some { data := #[
      m.get 1 1 * invDet, -m.get 1 0 * invDet,
      -m.get 0 1 * invDet, m.get 0 0 * invDet
    ]}

/-- LU decomposition result. -/
structure LU where
  L : Mat2
  U : Mat2
  perm : Array Nat
  permSign : Int
deriving Repr, Inhabited

/-- QR decomposition result. -/
structure QR where
  Q : Mat2
  R : Mat2
deriving Repr, Inhabited

/-- Cholesky decomposition result (lower-triangular). -/
structure Cholesky where
  L : Mat2
deriving Repr, Inhabited

/-- Eigen decomposition result for symmetric matrices. -/
structure Eigen where
  values : Vec2
  vectors : Mat2
deriving Repr, Inhabited

/-- Singular value decomposition result. -/
structure SVD where
  U : Mat2
  S : Vec2
  V : Mat2
deriving Repr, Inhabited

private def idx (row col : Nat) : Nat := col * 2 + row

private def getData (data : Array Float) (row col : Nat) : Float :=
  data.getD (idx row col) 0.0

private def setData (data : Array Float) (row col : Nat) (v : Float) : Array Float :=
  data.set! (idx row col) v

private def swapRows (data : Array Float) (i j : Nat) : Array Float := Id.run do
  let mut d := data
  for col in [:2] do
    let idx1 := idx i col
    let idx2 := idx j col
    let tmp := d[idx1]!
    d := d.set! idx1 d[idx2]!
    d := d.set! idx2 tmp
  return d

private def swapLowerRows (data : Array Float) (i j cols : Nat) : Array Float := Id.run do
  let mut d := data
  for col in [:cols] do
    let idx1 := idx i col
    let idx2 := idx j col
    let tmp := d[idx1]!
    d := d.set! idx1 d[idx2]!
    d := d.set! idx2 tmp
  return d

/-- LU decomposition with partial pivoting. -/
def luDecompose (m : Mat2) : Option LU := Id.run do
  let mut u := m.data
  let mut l := Mat2.identity.data
  let mut perm := #[0, 1]
  let mut sign : Int := 1
  for k in [:2] do
    let mut pivot := k
    let mut maxVal := Float.abs' (getData u k k)
    for i in [k+1:2] do
      let v := Float.abs' (getData u i k)
      if v > maxVal then
        maxVal := v
        pivot := i
    if maxVal < Float.epsilon then
      return none
    if pivot != k then
      u := swapRows u k pivot
      l := swapLowerRows l k pivot k
      let pk := perm[k]!
      let pp := perm[pivot]!
      perm := perm.set! k pp
      perm := perm.set! pivot pk
      sign := -sign
    for i in [k+1:2] do
      let ukk := getData u k k
      let factor := getData u i k / ukk
      l := setData l i k factor
      for j in [k:2] do
        let value := getData u i j - factor * getData u k j
        u := setData u i j value
  return some { L := { data := l }, U := { data := u }, perm := perm, permSign := sign }

private def permuteVec (perm : Array Nat) (b : Vec2) : Array Float :=
  let arr := #[b.x, b.y]
  #[arr[perm[0]!]!, arr[perm[1]!]!]

private def forwardSub (l : Array Float) (b : Array Float) : Array Float :=
  let y0 := b[0]!
  let y1 := b[1]! - getData l 1 0 * y0
  #[y0, y1]

private def backSub (u : Array Float) (y : Array Float) : Option (Array Float) :=
  let u11 := getData u 1 1
  if Float.abs' u11 < Float.epsilon then none else
    let x1 := y[1]! / u11
    let u00 := getData u 0 0
    if Float.abs' u00 < Float.epsilon then none else
      let x0 := (y[0]! - getData u 0 1 * x1) / u00
      some #[x0, x1]

/-- Solve Ax=b using LU decomposition. -/
def solveLU (lu : LU) (b : Vec2) : Option Vec2 :=
  let bp := permuteVec lu.perm b
  let y := forwardSub lu.L.data bp
  match backSub lu.U.data y with
  | none => none
  | some x => some (Vec2.mk x[0]! x[1]!)

/-- Solve Ax=b (uses LU). -/
def solve (m : Mat2) (b : Vec2) : Option Vec2 :=
  match m.luDecompose with
  | none => none
  | some lu => solveLU lu b

/-- QR decomposition (modified Gram-Schmidt). -/
def qrDecompose (m : Mat2) : Option QR :=
  let a0 := m.column 0
  let a1 := m.column 1
  let r00 := a0.length
  if r00 < Float.epsilon then none
  else
    let q0 := a0.scale (1.0 / r00)
    let r01 := q0.dot a1
    let v1 := a1 - q0.scale r01
    let r11 := v1.length
    if r11 < Float.epsilon then none
    else
      let q1 := v1.scale (1.0 / r11)
      let q := Mat2.fromColumns q0 q1
      let r := Mat2.fromRows (Vec2.mk r00 r01) (Vec2.mk 0.0 r11)
      some { Q := q, R := r }

/-- Solve Ax=b using QR decomposition. -/
def solveQR (qr : QR) (b : Vec2) : Option Vec2 :=
  let q0 := qr.Q.column 0
  let q1 := qr.Q.column 1
  let y0 := q0.dot b
  let y1 := q1.dot b
  let u := qr.R.data
  match backSub u #[y0, y1] with
  | none => none
  | some x => some (Vec2.mk x[0]! x[1]!)

/-- Solve least squares (uses QR). -/
def leastSquares (m : Mat2) (b : Vec2) : Option Vec2 :=
  match m.qrDecompose with
  | none => none
  | some qr => solveQR qr b

/-- Cholesky decomposition (A = L * Lᵀ). -/
def choleskyDecompose (m : Mat2) : Option Cholesky :=
  let a00 := m.get 0 0
  if a00 <= 0.0 then none
  else
    let l00 := Float.sqrt a00
    let l10 := m.get 1 0 / l00
    let a11 := m.get 1 1 - l10 * l10
    if a11 <= 0.0 then none
    else
      let l11 := Float.sqrt a11
      let l := Mat2.fromRows (Vec2.mk l00 0.0) (Vec2.mk l10 l11)
      some { L := l }

/-- Solve Ax=b using Cholesky decomposition. -/
def solveCholesky (chol : Cholesky) (b : Vec2) : Option Vec2 :=
  let l := chol.L
  let y0 := b.x / l.get 0 0
  let y1 := (b.y - l.get 1 0 * y0) / l.get 1 1
  let x1 := y1 / l.get 1 1
  let x0 := (y0 - l.get 1 0 * x1) / l.get 0 0
  some (Vec2.mk x0 x1)

/-- Eigen decomposition of a symmetric 2x2 matrix. -/
def eigenDecomposeSymmetric (m : Mat2) : Eigen :=
  let a := m.get 0 0
  let b := m.get 0 1
  let d := m.get 1 1
  let t := (a + d) * 0.5
  let diff := (a - d) * 0.5
  let s := Float.sqrt (diff * diff + b * b)
  let l0 := t + s
  let l1 := t - s
  if Float.abs' b < Float.epsilon then
    let v0 := if a >= d then Vec2.unitX else Vec2.unitY
    let v1 := if a >= d then Vec2.unitY else Vec2.unitX
    { values := Vec2.mk l0 l1, vectors := Mat2.fromColumns v0 v1 }
  else
    let x := b
    let y := l0 - a
    let len := Float.sqrt (x * x + y * y)
    let v0 :=
      if len > Float.epsilon then Vec2.mk (x / len) (y / len) else Vec2.unitX
    let v1 := Vec2.perpendicular v0
    { values := Vec2.mk l0 l1, vectors := Mat2.fromColumns v0 v1 }

/-- Singular value decomposition (A = U * diag(S) * Vᵀ). -/
def svd (m : Mat2) : SVD :=
  let ata := Mat2.multiply m.transpose m
  let eigen := ata.eigenDecomposeSymmetric
  let s0 := Float.sqrt (Float.max 0.0 eigen.values.x)
  let s1 := Float.sqrt (Float.max 0.0 eigen.values.y)
  let v0 := eigen.vectors.column 0
  let v1 := eigen.vectors.column 1
  let mulVec := fun v : Vec2 =>
    Vec2.mk
      (m.get 0 0 * v.x + m.get 0 1 * v.y)
      (m.get 1 0 * v.x + m.get 1 1 * v.y)
  let u0 :=
    if s0 > Float.epsilon then (mulVec v0).scale (1.0 / s0) else Vec2.unitX
  let u0n := u0.normalize
  let u1 :=
    if s1 > Float.epsilon then (mulVec v1).scale (1.0 / s1) else Vec2.perpendicular u0n
  let u1n := u1.normalize
  { U := Mat2.fromColumns u0n u1n, S := Vec2.mk s0 s1, V := eigen.vectors }

/-- Transform a vector by this matrix. -/
def transformVec2 (m : Mat2) (v : Vec2) : Vec2 :=
  ⟨m.get 0 0 * v.x + m.get 0 1 * v.y,
   m.get 1 0 * v.x + m.get 1 1 * v.y⟩

/-- Create a rotation matrix (angle in radians). -/
def rotation (angle : Float) : Mat2 :=
  let c := Float.cos angle
  let s := Float.sin angle
  { data := #[c, s, -s, c] }

/-- Create a scale matrix. -/
def scaling (sx sy : Float) : Mat2 :=
  { data := #[sx, 0, 0, sy] }

/-- Create a uniform scale matrix. -/
def scalingUniform (s : Float) : Mat2 := scaling s s

/-- Check if two matrices are approximately equal. -/
def approxEq (a b : Mat2) (eps : Float := Float.epsilon) : Bool :=
  a.data.size == b.data.size &&
  (Array.zip a.data b.data).all fun (x, y) => Float.approxEq x y eps

/-- Convert to flat array (column-major). -/
def toArray (m : Mat2) : Array Float := m.data

instance : Add Mat2 := ⟨add⟩
instance : Sub Mat2 := ⟨sub⟩
instance : HMul Mat2 Mat2 Mat2 := ⟨multiply⟩
instance : HMul Mat2 Vec2 Vec2 := ⟨transformVec2⟩
instance : HMul Mat2 Float Mat2 := ⟨scale⟩

end Mat2

end Linalg
