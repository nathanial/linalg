/-
  4x4 Matrix type and operations.
  Column-major storage for GPU compatibility (OpenGL/Metal convention).
-/

import Linalg.Core
import Linalg.Vec3
import Linalg.Vec4
import Linalg.Mat3

namespace Linalg

/-- 4x4 matrix in column-major order (like OpenGL/Metal). -/
structure Mat4 where
  data : Array Float
deriving Inhabited, Repr

namespace Mat4

/-- Identity matrix. -/
def identity : Mat4 := { data := #[
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
]}

/-- Zero matrix. -/
def zero : Mat4 := { data := #[
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0
]}

/-- Create matrix from column vectors. -/
def fromColumns (c0 c1 c2 c3 : Vec4) : Mat4 :=
  { data := #[
    c0.x, c0.y, c0.z, c0.w,
    c1.x, c1.y, c1.z, c1.w,
    c2.x, c2.y, c2.z, c2.w,
    c3.x, c3.y, c3.z, c3.w
  ]}

/-- Get element at (row, col). Column-major indexing. -/
@[inline]
def get (m : Mat4) (row col : Nat) : Float :=
  m.data.getD (col * 4 + row) 0.0

/-- Set element at (row, col). -/
def set (m : Mat4) (row col : Nat) (v : Float) : Mat4 :=
  { data := m.data.set! (col * 4 + row) v }

/-- Get column as Vec4. -/
def column (m : Mat4) (c : Nat) : Vec4 :=
  ⟨m.get 0 c, m.get 1 c, m.get 2 c, m.get 3 c⟩

/-- Get row as Vec4. -/
def row (m : Mat4) (r : Nat) : Vec4 :=
  ⟨m.get r 0, m.get r 1, m.get r 2, m.get r 3⟩

/-- Component-wise addition. -/
def add (a b : Mat4) : Mat4 :=
  let result := Id.run do
    let mut arr : Array Float := #[]
    for i in [:a.data.size] do
      arr := arr.push (a.data.getD i 0.0 + b.data.getD i 0.0)
    return arr
  { data := result }

/-- Component-wise subtraction. -/
def sub (a b : Mat4) : Mat4 :=
  let result := Id.run do
    let mut arr : Array Float := #[]
    for i in [:a.data.size] do
      arr := arr.push (a.data.getD i 0.0 - b.data.getD i 0.0)
    return arr
  { data := result }

/-- Scale all elements by a scalar. -/
def scale (m : Mat4) (s : Float) : Mat4 :=
  { data := m.data.map (· * s) }

/-- Matrix multiplication. -/
def multiply (a b : Mat4) : Mat4 :=
  let result := Id.run do
    let mut arr : Array Float := #[]
    for col in [:4] do
      for row in [:4] do
        let mut sum := 0.0
        for k in [:4] do
          sum := sum + a.get row k * b.get k col
        arr := arr.push sum
    return arr
  { data := result }

/-- Transpose the matrix. -/
def transpose (m : Mat4) : Mat4 :=
  let result := Id.run do
    let mut arr : Array Float := #[]
    for col in [:4] do
      for row in [:4] do
        arr := arr.push (m.get col row)
    return arr
  { data := result }

/-- Compute the determinant. -/
def determinant (m : Mat4) : Float :=
  let a := m.get 0 0; let b := m.get 0 1; let c := m.get 0 2; let d := m.get 0 3
  let e := m.get 1 0; let f := m.get 1 1; let g := m.get 1 2; let h := m.get 1 3
  let i := m.get 2 0; let j := m.get 2 1; let k := m.get 2 2; let l := m.get 2 3
  let mm := m.get 3 0; let n := m.get 3 1; let o := m.get 3 2; let p := m.get 3 3

  let kp_lo := k * p - l * o
  let jp_ln := j * p - l * n
  let jo_kn := j * o - k * n
  let ip_lm := i * p - l * mm
  let io_km := i * o - k * mm
  let in_jm := i * n - j * mm

  a * (f * kp_lo - g * jp_ln + h * jo_kn) -
  b * (e * kp_lo - g * ip_lm + h * io_km) +
  c * (e * jp_ln - f * ip_lm + h * in_jm) -
  d * (e * jo_kn - f * io_km + g * in_jm)

/-- Compute the inverse. Returns none if matrix is singular. -/
def inverse (m : Mat4) : Option Mat4 :=
  let det := m.determinant
  if Float.abs' det < Float.epsilon then none
  else
    let a := m.get 0 0; let b := m.get 0 1; let c := m.get 0 2; let d := m.get 0 3
    let e := m.get 1 0; let f := m.get 1 1; let g := m.get 1 2; let h := m.get 1 3
    let i := m.get 2 0; let j := m.get 2 1; let k := m.get 2 2; let l := m.get 2 3
    let mm := m.get 3 0; let n := m.get 3 1; let o := m.get 3 2; let p := m.get 3 3

    let kp_lo := k * p - l * o
    let jp_ln := j * p - l * n
    let jo_kn := j * o - k * n
    let ip_lm := i * p - l * mm
    let io_km := i * o - k * mm
    let in_jm := i * n - j * mm
    let gp_ho := g * p - h * o
    let fp_hn := f * p - h * n
    let fo_gn := f * o - g * n
    let ep_hm := e * p - h * mm
    let eo_gm := e * o - g * mm
    let en_fm := e * n - f * mm
    let gl_hk := g * l - h * k
    let fl_hj := f * l - h * j
    let fk_gj := f * k - g * j
    let el_hi := e * l - h * i
    let ek_gi := e * k - g * i
    let ej_fi := e * j - f * i

    let invDet := 1.0 / det

    some { data := #[
      (f * kp_lo - g * jp_ln + h * jo_kn) * invDet,
      -(e * kp_lo - g * ip_lm + h * io_km) * invDet,
      (e * jp_ln - f * ip_lm + h * in_jm) * invDet,
      -(e * jo_kn - f * io_km + g * in_jm) * invDet,

      -(b * kp_lo - c * jp_ln + d * jo_kn) * invDet,
      (a * kp_lo - c * ip_lm + d * io_km) * invDet,
      -(a * jp_ln - b * ip_lm + d * in_jm) * invDet,
      (a * jo_kn - b * io_km + c * in_jm) * invDet,

      (b * gp_ho - c * fp_hn + d * fo_gn) * invDet,
      -(a * gp_ho - c * ep_hm + d * eo_gm) * invDet,
      (a * fp_hn - b * ep_hm + d * en_fm) * invDet,
      -(a * fo_gn - b * eo_gm + c * en_fm) * invDet,

      -(b * gl_hk - c * fl_hj + d * fk_gj) * invDet,
      (a * gl_hk - c * el_hi + d * ek_gi) * invDet,
      -(a * fl_hj - b * el_hi + d * ej_fi) * invDet,
      (a * fk_gj - b * ek_gi + c * ej_fi) * invDet
    ]}

/-- LU decomposition result. -/
structure LU where
  L : Mat4
  U : Mat4
  perm : Array Nat
  permSign : Int
deriving Repr, Inhabited

/-- QR decomposition result. -/
structure QR where
  Q : Mat4
  R : Mat4
deriving Repr, Inhabited

/-- Cholesky decomposition result (lower-triangular). -/
structure Cholesky where
  L : Mat4
deriving Repr, Inhabited

private def idx (row col : Nat) : Nat := col * 4 + row

private def getData (data : Array Float) (row col : Nat) : Float :=
  data.getD (idx row col) 0.0

private def setData (data : Array Float) (row col : Nat) (v : Float) : Array Float :=
  data.set! (idx row col) v

private def swapRows (data : Array Float) (i j : Nat) : Array Float := Id.run do
  let mut d := data
  for col in [:4] do
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
def luDecompose (m : Mat4) : Option LU := Id.run do
  let mut u := m.data
  let mut l := Mat4.identity.data
  let mut perm := #[0, 1, 2, 3]
  let mut sign : Int := 1
  for k in [:4] do
    let mut pivot := k
    let mut maxVal := Float.abs' (getData u k k)
    for i in [k+1:4] do
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
    for i in [k+1:4] do
      let ukk := getData u k k
      let factor := getData u i k / ukk
      l := setData l i k factor
      for j in [k:4] do
        let value := getData u i j - factor * getData u k j
        u := setData u i j value
  return some { L := { data := l }, U := { data := u }, perm := perm, permSign := sign }

private def permuteVec (perm : Array Nat) (b : Vec4) : Array Float :=
  let arr := #[b.x, b.y, b.z, b.w]
  #[arr[perm[0]!]!, arr[perm[1]!]!, arr[perm[2]!]!, arr[perm[3]!]!]

private def forwardSub (l : Array Float) (b : Array Float) : Array Float := Id.run do
  let mut y := Array.replicate 4 0.0
  for i in [:4] do
    let mut sum := b[i]!
    for j in [:i] do
      sum := sum - getData l i j * y[j]!
    y := y.set! i sum
  return y

private def backSub (u : Array Float) (y : Array Float) : Option (Array Float) := Id.run do
  let mut x := Array.replicate 4 0.0
  for offset in [:4] do
    let i := 3 - offset
    let mut sum := y[i]!
    for j in [i+1:4] do
      sum := sum - getData u i j * x[j]!
    let diag := getData u i i
    if Float.abs' diag < Float.epsilon then
      return none
    x := x.set! i (sum / diag)
  return some x

/-- Solve Ax=b using LU decomposition. -/
def solveLU (lu : LU) (b : Vec4) : Option Vec4 :=
  let bp := permuteVec lu.perm b
  let y := forwardSub lu.L.data bp
  match backSub lu.U.data y with
  | none => none
  | some x => some (Vec4.mk x[0]! x[1]! x[2]! x[3]!)

/-- Solve Ax=b (uses LU). -/
def solve (m : Mat4) (b : Vec4) : Option Vec4 :=
  match m.luDecompose with
  | none => none
  | some lu => solveLU lu b

/-- QR decomposition (modified Gram-Schmidt). -/
def qrDecompose (m : Mat4) : Option QR :=
  let a0 := m.column 0
  let a1 := m.column 1
  let a2 := m.column 2
  let a3 := m.column 3
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
        let r03 := q0.dot a3
        let r13 := q1.dot a3
        let r23 := q2.dot a3
        let v3 := a3 - q0.scale r03 - q1.scale r13 - q2.scale r23
        let r33 := v3.length
        if r33 < Float.epsilon then none
        else
          let q3 := v3.scale (1.0 / r33)
          let q := Mat4.fromColumns q0 q1 q2 q3
          let r := Id.run do
            let mut r := Mat4.zero
            r := r.set 0 0 r00
            r := r.set 0 1 r01
            r := r.set 0 2 r02
            r := r.set 0 3 r03
            r := r.set 1 1 r11
            r := r.set 1 2 r12
            r := r.set 1 3 r13
            r := r.set 2 2 r22
            r := r.set 2 3 r23
            r := r.set 3 3 r33
            return r
          some { Q := q, R := r }

/-- Solve Ax=b using QR decomposition. -/
def solveQR (qr : QR) (b : Vec4) : Option Vec4 :=
  let q0 := qr.Q.column 0
  let q1 := qr.Q.column 1
  let q2 := qr.Q.column 2
  let q3 := qr.Q.column 3
  let y0 := q0.dot b
  let y1 := q1.dot b
  let y2 := q2.dot b
  let y3 := q3.dot b
  match backSub qr.R.data #[y0, y1, y2, y3] with
  | none => none
  | some x => some (Vec4.mk x[0]! x[1]! x[2]! x[3]!)

/-- Solve least squares (uses QR). -/
def leastSquares (m : Mat4) (b : Vec4) : Option Vec4 :=
  match m.qrDecompose with
  | none => none
  | some qr => solveQR qr b

/-- Cholesky decomposition (A = L * Lᵀ). -/
def choleskyDecompose (m : Mat4) : Option Cholesky := Id.run do
  let mut l := Mat4.zero.data
  for i in [:4] do
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
def solveCholesky (chol : Cholesky) (b : Vec4) : Option Vec4 := Id.run do
  let l := chol.L.data
  let bArr := #[b.x, b.y, b.z, b.w]
  let mut y := Array.replicate 4 0.0
  for i in [:4] do
    let mut sum := bArr[i]!
    for j in [:i] do
      sum := sum - getData l i j * y[j]!
    let diag := getData l i i
    if Float.abs' diag < Float.epsilon then return none
    y := y.set! i (sum / diag)
  let mut x := Array.replicate 4 0.0
  for offset in [:4] do
    let i := 3 - offset
    let mut sum := y[i]!
    for j in [i+1:4] do
      sum := sum - getData l j i * x[j]!
    let diag := getData l i i
    if Float.abs' diag < Float.epsilon then return none
    x := x.set! i (sum / diag)
  return some (Vec4.mk x[0]! x[1]! x[2]! x[3]!)

/-- Transform a Vec4 by this matrix. -/
def transformVec4 (m : Mat4) (v : Vec4) : Vec4 :=
  ⟨m.get 0 0 * v.x + m.get 0 1 * v.y + m.get 0 2 * v.z + m.get 0 3 * v.w,
   m.get 1 0 * v.x + m.get 1 1 * v.y + m.get 1 2 * v.z + m.get 1 3 * v.w,
   m.get 2 0 * v.x + m.get 2 1 * v.y + m.get 2 2 * v.z + m.get 2 3 * v.w,
   m.get 3 0 * v.x + m.get 3 1 * v.y + m.get 3 2 * v.z + m.get 3 3 * v.w⟩

/-- Transform a point (w=1) with perspective divide. -/
def transformPoint (m : Mat4) (v : Vec3) : Vec3 :=
  let result := m.transformVec4 (Vec4.fromPoint v)
  result.toVec3Normalized

/-- Transform a direction vector (w=0, no translation). -/
def transformDirection (m : Mat4) (v : Vec3) : Vec3 :=
  let result := m.transformVec4 (Vec4.fromDirection v)
  result.toVec3

/-- Create a translation matrix. -/
def translation (x y z : Float) : Mat4 := { data := #[
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  x, y, z, 1
]}

/-- Create a translation matrix from a vector. -/
def translationV (v : Vec3) : Mat4 := translation v.x v.y v.z

/-- Create a rotation matrix around X axis (angle in radians). -/
def rotationX (angle : Float) : Mat4 :=
  let c := Float.cos angle
  let s := Float.sin angle
  { data := #[
    1, 0, 0, 0,
    0, c, s, 0,
    0, -s, c, 0,
    0, 0, 0, 1
  ]}

/-- Create a rotation matrix around Y axis (angle in radians). -/
def rotationY (angle : Float) : Mat4 :=
  let c := Float.cos angle
  let s := Float.sin angle
  { data := #[
    c, 0, -s, 0,
    0, 1, 0, 0,
    s, 0, c, 0,
    0, 0, 0, 1
  ]}

/-- Create a rotation matrix around Z axis (angle in radians). -/
def rotationZ (angle : Float) : Mat4 :=
  let c := Float.cos angle
  let s := Float.sin angle
  { data := #[
    c, s, 0, 0,
    -s, c, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  ]}

/-- Create a rotation matrix around arbitrary axis (angle in radians). -/
def rotationAxisAngle (axis : Vec3) (angle : Float) : Mat4 :=
  let n := axis.normalize
  let c := Float.cos angle
  let s := Float.sin angle
  let t := 1.0 - c
  { data := #[
    t * n.x * n.x + c,       t * n.x * n.y + s * n.z, t * n.x * n.z - s * n.y, 0,
    t * n.x * n.y - s * n.z, t * n.y * n.y + c,       t * n.y * n.z + s * n.x, 0,
    t * n.x * n.z + s * n.y, t * n.y * n.z - s * n.x, t * n.z * n.z + c,       0,
    0,                       0,                       0,                       1
  ]}

/-- Create a scale matrix. -/
def scaling (x y z : Float) : Mat4 := { data := #[
  x, 0, 0, 0,
  0, y, 0, 0,
  0, 0, z, 0,
  0, 0, 0, 1
]}

/-- Create a scale matrix from a vector. -/
def scalingV (v : Vec3) : Mat4 := scaling v.x v.y v.z

/-- Create a uniform scale matrix. -/
def scalingUniform (s : Float) : Mat4 := scaling s s s

/-- Create a perspective projection matrix. -/
def perspective (fovY aspect near far : Float) : Mat4 :=
  let tanHalfFov := Float.tan (fovY / 2.0)
  let f := 1.0 / tanHalfFov
  let nf := 1.0 / (near - far)
  { data := #[
    f / aspect, 0, 0, 0,
    0, f, 0, 0,
    0, 0, (far + near) * nf, -1,
    0, 0, 2 * far * near * nf, 0
  ]}

/-- Create an orthographic projection matrix. -/
def orthographic (left right bottom top near far : Float) : Mat4 :=
  let rl := 1.0 / (right - left)
  let tb := 1.0 / (top - bottom)
  let fn := 1.0 / (far - near)
  { data := #[
    2 * rl, 0, 0, 0,
    0, 2 * tb, 0, 0,
    0, 0, -2 * fn, 0,
    -(right + left) * rl, -(top + bottom) * tb, -(far + near) * fn, 1
  ]}

/-- Create a look-at view matrix. -/
def lookAt (eye center up : Vec3) : Mat4 :=
  let f := (center.sub eye).normalize
  let s := f.cross up |>.normalize
  let u := s.cross f
  { data := #[
    s.x, u.x, -f.x, 0,
    s.y, u.y, -f.y, 0,
    s.z, u.z, -f.z, 0,
    -(s.dot eye), -(u.dot eye), f.dot eye, 1
  ]}

/-- Get the translation component. -/
def getTranslation (m : Mat4) : Vec3 := ⟨m.get 0 3, m.get 1 3, m.get 2 3⟩

/-- Get the scale component (assumes no shear). -/
def getScale (m : Mat4) : Vec3 :=
  let sx := Vec3.mk (m.get 0 0) (m.get 1 0) (m.get 2 0) |>.length
  let sy := Vec3.mk (m.get 0 1) (m.get 1 1) (m.get 2 1) |>.length
  let sz := Vec3.mk (m.get 0 2) (m.get 1 2) (m.get 2 2) |>.length
  ⟨sx, sy, sz⟩

/-- Check if two matrices are approximately equal. -/
def approxEq (a b : Mat4) (eps : Float := Float.epsilon) : Bool :=
  a.data.size == b.data.size &&
  (Array.zip a.data b.data).all fun (x, y) => Float.approxEq x y eps

/-- Convert to flat array (column-major). -/
def toArray (m : Mat4) : Array Float := m.data

/-- Extract upper-left 3x3 submatrix. -/
def toMat3 (m : Mat4) : Mat3 :=
  { data := #[
    m.get 0 0, m.get 1 0, m.get 2 0,
    m.get 0 1, m.get 1 1, m.get 2 1,
    m.get 0 2, m.get 1 2, m.get 2 2
  ]}

/-- Create Mat4 from Mat3 (embed in upper-left, 1 on w diagonal). -/
def fromMat3 (m : Mat3) : Mat4 :=
  { data := #[
    m.get 0 0, m.get 1 0, m.get 2 0, 0.0,
    m.get 0 1, m.get 1 1, m.get 2 1, 0.0,
    m.get 0 2, m.get 1 2, m.get 2 2, 0.0,
    0.0,       0.0,       0.0,       1.0
  ]}

instance : Add Mat4 := ⟨add⟩
instance : Sub Mat4 := ⟨sub⟩
instance : HMul Mat4 Mat4 Mat4 := ⟨multiply⟩
instance : HMul Mat4 Vec4 Vec4 := ⟨transformVec4⟩
instance : HMul Mat4 Float Mat4 := ⟨scale⟩

/-- Coerce Mat3 to Mat4 (embed in upper-left, 1 on w diagonal). -/
instance : Coe Mat3 Mat4 := ⟨fromMat3⟩

end Mat4

end Linalg
