/-
  3x3 Matrix type and operations.
  Column-major storage for consistency with GPU conventions.
-/

import Linalg.Core
import Linalg.Vec3
import Linalg.Mat2

namespace Linalg

/-- 3x3 matrix in column-major order. -/
structure Mat3 where
  data : Array Float
deriving Inhabited, Repr

namespace Mat3

/-- Identity matrix. -/
def identity : Mat3 := { data := #[1, 0, 0, 0, 1, 0, 0, 0, 1] }

/-- Zero matrix. -/
def zero : Mat3 := { data := #[0, 0, 0, 0, 0, 0, 0, 0, 0] }

/-- Create matrix from column vectors. -/
def fromColumns (c0 c1 c2 : Vec3) : Mat3 :=
  { data := #[c0.x, c0.y, c0.z, c1.x, c1.y, c1.z, c2.x, c2.y, c2.z] }

/-- Get element at (row, col). -/
@[inline]
def get (m : Mat3) (row col : Nat) : Float :=
  m.data.getD (col * 3 + row) 0.0

/-- Set element at (row, col). -/
def set (m : Mat3) (row col : Nat) (v : Float) : Mat3 :=
  { data := m.data.set! (col * 3 + row) v }

/-- Get column as Vec3. -/
def column (m : Mat3) (c : Nat) : Vec3 :=
  ⟨m.get 0 c, m.get 1 c, m.get 2 c⟩

/-- Get row as Vec3. -/
def row (m : Mat3) (r : Nat) : Vec3 :=
  ⟨m.get r 0, m.get r 1, m.get r 2⟩

/-- Component-wise addition. -/
def add (a b : Mat3) : Mat3 :=
  let result := Id.run do
    let mut arr : Array Float := #[]
    for i in [:a.data.size] do
      arr := arr.push (a.data.getD i 0.0 + b.data.getD i 0.0)
    return arr
  { data := result }

/-- Component-wise subtraction. -/
def sub (a b : Mat3) : Mat3 :=
  let result := Id.run do
    let mut arr : Array Float := #[]
    for i in [:a.data.size] do
      arr := arr.push (a.data.getD i 0.0 - b.data.getD i 0.0)
    return arr
  { data := result }

/-- Scale all elements by a scalar. -/
def scale (m : Mat3) (s : Float) : Mat3 :=
  { data := m.data.map (· * s) }

/-- Matrix multiplication. -/
def multiply (a b : Mat3) : Mat3 :=
  let result := Id.run do
    let mut arr : Array Float := #[]
    for col in [:3] do
      for row in [:3] do
        let mut sum := 0.0
        for k in [:3] do
          sum := sum + a.get row k * b.get k col
        arr := arr.push sum
    return arr
  { data := result }

/-- Transpose the matrix. -/
def transpose (m : Mat3) : Mat3 :=
  { data := #[
    m.get 0 0, m.get 0 1, m.get 0 2,
    m.get 1 0, m.get 1 1, m.get 1 2,
    m.get 2 0, m.get 2 1, m.get 2 2
  ]}

/-- Compute the determinant. -/
def determinant (m : Mat3) : Float :=
  m.get 0 0 * (m.get 1 1 * m.get 2 2 - m.get 1 2 * m.get 2 1) -
  m.get 0 1 * (m.get 1 0 * m.get 2 2 - m.get 1 2 * m.get 2 0) +
  m.get 0 2 * (m.get 1 0 * m.get 2 1 - m.get 1 1 * m.get 2 0)

/-- Compute the inverse using adjugate/determinant method. -/
def inverse (m : Mat3) : Option Mat3 :=
  let det := m.determinant
  if Float.abs' det < Float.epsilon then none
  else
    let invDet := 1.0 / det
    -- Compute cofactor matrix (transposed = adjugate)
    let c00 := m.get 1 1 * m.get 2 2 - m.get 1 2 * m.get 2 1
    let c01 := -(m.get 1 0 * m.get 2 2 - m.get 1 2 * m.get 2 0)
    let c02 := m.get 1 0 * m.get 2 1 - m.get 1 1 * m.get 2 0
    let c10 := -(m.get 0 1 * m.get 2 2 - m.get 0 2 * m.get 2 1)
    let c11 := m.get 0 0 * m.get 2 2 - m.get 0 2 * m.get 2 0
    let c12 := -(m.get 0 0 * m.get 2 1 - m.get 0 1 * m.get 2 0)
    let c20 := m.get 0 1 * m.get 1 2 - m.get 0 2 * m.get 1 1
    let c21 := -(m.get 0 0 * m.get 1 2 - m.get 0 2 * m.get 1 0)
    let c22 := m.get 0 0 * m.get 1 1 - m.get 0 1 * m.get 1 0
    -- Adjugate is transpose of cofactor, then scale by 1/det
    some { data := #[
      c00 * invDet, c01 * invDet, c02 * invDet,
      c10 * invDet, c11 * invDet, c12 * invDet,
      c20 * invDet, c21 * invDet, c22 * invDet
    ]}

/-- LU decomposition result. -/
structure LU where
  L : Mat3
  U : Mat3
  perm : Array Nat
  permSign : Int
deriving Repr, Inhabited

/-- QR decomposition result. -/
structure QR where
  Q : Mat3
  R : Mat3
deriving Repr, Inhabited

/-- Cholesky decomposition result (lower-triangular). -/
structure Cholesky where
  L : Mat3
deriving Repr, Inhabited

/-- Eigen decomposition result for symmetric matrices. -/
structure Eigen where
  values : Vec3
  vectors : Mat3
deriving Repr, Inhabited

/-- Singular value decomposition result. -/
structure SVD where
  U : Mat3
  S : Vec3
  V : Mat3
deriving Repr, Inhabited

private def idx (row col : Nat) : Nat := col * 3 + row

private def getData (data : Array Float) (row col : Nat) : Float :=
  data.getD (idx row col) 0.0

private def setData (data : Array Float) (row col : Nat) (v : Float) : Array Float :=
  data.set! (idx row col) v

private def swapRows (data : Array Float) (i j : Nat) : Array Float := Id.run do
  let mut d := data
  for col in [:3] do
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
def luDecompose (m : Mat3) : Option LU := Id.run do
  let mut u := m.data
  let mut l := Mat3.identity.data
  let mut perm := #[0, 1, 2]
  let mut sign : Int := 1
  for k in [:3] do
    let mut pivot := k
    let mut maxVal := Float.abs' (getData u k k)
    for i in [k+1:3] do
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
    for i in [k+1:3] do
      let ukk := getData u k k
      let factor := getData u i k / ukk
      l := setData l i k factor
      for j in [k:3] do
        let value := getData u i j - factor * getData u k j
        u := setData u i j value
  return some { L := { data := l }, U := { data := u }, perm := perm, permSign := sign }

private def permuteVec (perm : Array Nat) (b : Vec3) : Array Float :=
  let arr := #[b.x, b.y, b.z]
  #[arr[perm[0]!]!, arr[perm[1]!]!, arr[perm[2]!]!]

private def forwardSub (l : Array Float) (b : Array Float) : Array Float := Id.run do
  let mut y := Array.replicate 3 0.0
  for i in [:3] do
    let mut sum := b[i]!
    for j in [:i] do
      sum := sum - getData l i j * y[j]!
    y := y.set! i sum
  return y

private def backSub (u : Array Float) (y : Array Float) : Option (Array Float) := Id.run do
  let mut x := Array.replicate 3 0.0
  for offset in [:3] do
    let i := 2 - offset
    let mut sum := y[i]!
    for j in [i+1:3] do
      sum := sum - getData u i j * x[j]!
    let diag := getData u i i
    if Float.abs' diag < Float.epsilon then
      return none
    x := x.set! i (sum / diag)
  return some x

/-- Solve Ax=b using LU decomposition. -/
def solveLU (lu : LU) (b : Vec3) : Option Vec3 :=
  let bp := permuteVec lu.perm b
  let y := forwardSub lu.L.data bp
  match backSub lu.U.data y with
  | none => none
  | some x => some (Vec3.mk x[0]! x[1]! x[2]!)

/-- Solve Ax=b (uses LU). -/
def solve (m : Mat3) (b : Vec3) : Option Vec3 :=
  match m.luDecompose with
  | none => none
  | some lu => solveLU lu b

/-- QR decomposition (modified Gram-Schmidt). -/
def qrDecompose (m : Mat3) : Option QR :=
  let a0 := m.column 0
  let a1 := m.column 1
  let a2 := m.column 2
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
      let r02 := q0.dot a2
      let r12 := q1.dot a2
      let v2 := a2 - q0.scale r02 - q1.scale r12
      let r22 := v2.length
      if r22 < Float.epsilon then none
      else
      let q2 := v2.scale (1.0 / r22)
      let q := Mat3.fromColumns q0 q1 q2
      let r := Id.run do
        let mut r := Mat3.zero
        r := r.set 0 0 r00
        r := r.set 0 1 r01
        r := r.set 0 2 r02
        r := r.set 1 1 r11
        r := r.set 1 2 r12
        r := r.set 2 2 r22
        return r
      some { Q := q, R := r }

/-- Solve Ax=b using QR decomposition. -/
def solveQR (qr : QR) (b : Vec3) : Option Vec3 :=
  let q0 := qr.Q.column 0
  let q1 := qr.Q.column 1
  let q2 := qr.Q.column 2
  let y0 := q0.dot b
  let y1 := q1.dot b
  let y2 := q2.dot b
  match backSub qr.R.data #[y0, y1, y2] with
  | none => none
  | some x => some (Vec3.mk x[0]! x[1]! x[2]!)

/-- Solve least squares (uses QR). -/
def leastSquares (m : Mat3) (b : Vec3) : Option Vec3 :=
  match m.qrDecompose with
  | none => none
  | some qr => solveQR qr b

/-- Cholesky decomposition (A = L * Lᵀ). -/
def choleskyDecompose (m : Mat3) : Option Cholesky := Id.run do
  let mut l := Mat3.zero.data
  for i in [:3] do
    for j in [:i+1] do
      let mut sum := m.get i j
      for k in [:j] do
        sum := sum - getData l i k * getData l j k
      if i == j then
        if sum <= 0.0 then return none
        l := setData l i j (Float.sqrt sum)
      else
        let ljj := getData l j j
        if Float.abs' ljj < Float.epsilon then return none
        l := setData l i j (sum / ljj)
  return some { L := { data := l } }

/-- Solve Ax=b using Cholesky decomposition. -/
def solveCholesky (chol : Cholesky) (b : Vec3) : Option Vec3 := Id.run do
  let l := chol.L.data
  let bArr := #[b.x, b.y, b.z]
  let mut y := Array.replicate 3 0.0
  for i in [:3] do
    let mut sum := bArr[i]!
    for j in [:i] do
      sum := sum - getData l i j * y[j]!
    let diag := getData l i i
    if Float.abs' diag < Float.epsilon then return none
    y := y.set! i (sum / diag)
  let mut x := Array.replicate 3 0.0
  for offset in [:3] do
    let i := 2 - offset
    let mut sum := y[i]!
    for j in [i+1:3] do
      sum := sum - getData l j i * x[j]!
    let diag := getData l i i
    if Float.abs' diag < Float.epsilon then return none
    x := x.set! i (sum / diag)
  return some (Vec3.mk x[0]! x[1]! x[2]!)

/-- Eigen decomposition of a symmetric 3x3 matrix (Jacobi iterations). -/
def eigenDecomposeSymmetric (m : Mat3) (iterations : Nat := 20) : Eigen := Id.run do
  let mut a := m.data
  let mut v := Mat3.identity.data
  for _ in [:iterations] do
    let mut p := 0
    let mut q := 1
    let mut maxVal := Float.abs' (getData a 0 1)
    for i in [:3] do
      for j in [i+1:3] do
        let value := Float.abs' (getData a i j)
        if value > maxVal then
          maxVal := value
          p := i
          q := j
    if maxVal < Float.epsilon then
      break
    let app := getData a p p
    let aqq := getData a q q
    let apq := getData a p q
    let tau := (aqq - app) / (2.0 * apq)
    let t :=
      if tau >= 0.0 then
        1.0 / (tau + Float.sqrt (1.0 + tau * tau))
      else
        -1.0 / (-tau + Float.sqrt (1.0 + tau * tau))
    let c := 1.0 / Float.sqrt (1.0 + t * t)
    let s := t * c
    for k in [:3] do
      if k != p && k != q then
        let akp := getData a k p
        let akq := getData a k q
        let newAkp := c * akp - s * akq
        let newAkq := s * akp + c * akq
        a := setData a k p newAkp
        a := setData a p k newAkp
        a := setData a k q newAkq
        a := setData a q k newAkq
    let appNew := c * c * app - 2.0 * s * c * apq + s * s * aqq
    let aqqNew := s * s * app + 2.0 * s * c * apq + c * c * aqq
    a := setData a p p appNew
    a := setData a q q aqqNew
    a := setData a p q 0.0
    a := setData a q p 0.0
    for k in [:3] do
      let vkp := getData v k p
      let vkq := getData v k q
      let newVkp := c * vkp - s * vkq
      let newVkq := s * vkp + c * vkq
      v := setData v k p newVkp
      v := setData v k q newVkq
  let vMat : Mat3 := { data := v }
  let cols := #[
    (getData a 0 0, vMat.column 0),
    (getData a 1 1, vMat.column 1),
    (getData a 2 2, vMat.column 2)
  ]
  let sorted := cols.qsort (fun a b => a.1 > b.1)
  let v0 := sorted[0]!.2
  let v1 := sorted[1]!.2
  let v2 := sorted[2]!.2
  let values := Vec3.mk sorted[0]!.1 sorted[1]!.1 sorted[2]!.1
  return { values := values, vectors := Mat3.fromColumns v0 v1 v2 }

private def orthonormalVector (v : Vec3) : Vec3 :=
  let axis := if Float.abs' v.x < 0.9 then Vec3.unitX else Vec3.unitY
  (v.cross axis).normalize

/-- Singular value decomposition (A = U * diag(S) * Vᵀ). -/
def svd (m : Mat3) : SVD :=
  let ata := Mat3.multiply m.transpose m
  let eigen := ata.eigenDecomposeSymmetric
  let s0 := Float.sqrt (Float.max 0.0 eigen.values.x)
  let s1 := Float.sqrt (Float.max 0.0 eigen.values.y)
  let s2 := Float.sqrt (Float.max 0.0 eigen.values.z)
  let v0 := eigen.vectors.column 0
  let v1 := eigen.vectors.column 1
  let v2 := eigen.vectors.column 2
  let mulVec := fun v : Vec3 =>
    Vec3.mk
      (m.get 0 0 * v.x + m.get 0 1 * v.y + m.get 0 2 * v.z)
      (m.get 1 0 * v.x + m.get 1 1 * v.y + m.get 1 2 * v.z)
      (m.get 2 0 * v.x + m.get 2 1 * v.y + m.get 2 2 * v.z)
  let u0 :=
    if s0 > Float.epsilon then (mulVec v0).scale (1.0 / s0) else Vec3.unitX
  let u0n := u0.normalize
  let u1raw :=
    if s1 > Float.epsilon then (mulVec v1).scale (1.0 / s1) else orthonormalVector u0n
  let u1orth := u1raw - u0n.scale (u1raw.dot u0n)
  let u1n := u1orth.normalize
  let u2raw :=
    if s2 > Float.epsilon then (mulVec v2).scale (1.0 / s2) else u0n.cross u1n
  let u2orth := u2raw - u0n.scale (u2raw.dot u0n) - u1n.scale (u2raw.dot u1n)
  let u2n := u2orth.normalize
  { U := Mat3.fromColumns u0n u1n u2n, S := Vec3.mk s0 s1 s2, V := eigen.vectors }

/-- Transform a vector by this matrix. -/
def transformVec3 (m : Mat3) (v : Vec3) : Vec3 :=
  ⟨m.get 0 0 * v.x + m.get 0 1 * v.y + m.get 0 2 * v.z,
   m.get 1 0 * v.x + m.get 1 1 * v.y + m.get 1 2 * v.z,
   m.get 2 0 * v.x + m.get 2 1 * v.y + m.get 2 2 * v.z⟩

/-- Create rotation matrix around X axis (angle in radians). -/
def rotationX (angle : Float) : Mat3 :=
  let c := Float.cos angle
  let s := Float.sin angle
  { data := #[
    1, 0, 0,
    0, c, s,
    0, -s, c
  ]}

/-- Create rotation matrix around Y axis (angle in radians). -/
def rotationY (angle : Float) : Mat3 :=
  let c := Float.cos angle
  let s := Float.sin angle
  { data := #[
    c, 0, -s,
    0, 1, 0,
    s, 0, c
  ]}

/-- Create rotation matrix around Z axis (angle in radians). -/
def rotationZ (angle : Float) : Mat3 :=
  let c := Float.cos angle
  let s := Float.sin angle
  { data := #[
    c, s, 0,
    -s, c, 0,
    0, 0, 1
  ]}

/-- Create rotation matrix around arbitrary axis (angle in radians). -/
def rotationAxisAngle (axis : Vec3) (angle : Float) : Mat3 :=
  let n := axis.normalize
  let c := Float.cos angle
  let s := Float.sin angle
  let t := 1.0 - c
  { data := #[
    t * n.x * n.x + c,       t * n.x * n.y + s * n.z, t * n.x * n.z - s * n.y,
    t * n.x * n.y - s * n.z, t * n.y * n.y + c,       t * n.y * n.z + s * n.x,
    t * n.x * n.z + s * n.y, t * n.y * n.z - s * n.x, t * n.z * n.z + c
  ]}

/-- Create a scale matrix. -/
def scaling (sx sy sz : Float) : Mat3 :=
  { data := #[sx, 0, 0, 0, sy, 0, 0, 0, sz] }

/-- Create a uniform scale matrix. -/
def scalingUniform (s : Float) : Mat3 := scaling s s s

/-- Check if two matrices are approximately equal. -/
def approxEq (a b : Mat3) (eps : Float := Float.epsilon) : Bool :=
  a.data.size == b.data.size &&
  (Array.zip a.data b.data).all fun (x, y) => Float.approxEq x y eps

/-- Convert to flat array (column-major). -/
def toArray (m : Mat3) : Array Float := m.data

/-- Extract upper-left 2x2 submatrix. -/
def toMat2 (m : Mat3) : Mat2 :=
  { data := #[m.get 0 0, m.get 1 0, m.get 0 1, m.get 1 1] }

instance : Add Mat3 := ⟨add⟩
instance : Sub Mat3 := ⟨sub⟩
instance : HMul Mat3 Mat3 Mat3 := ⟨multiply⟩
instance : HMul Mat3 Vec3 Vec3 := ⟨transformVec3⟩
instance : HMul Mat3 Float Mat3 := ⟨scale⟩

end Mat3

end Linalg
