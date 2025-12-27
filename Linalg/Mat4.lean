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

instance : Add Mat4 := ⟨add⟩
instance : Sub Mat4 := ⟨sub⟩
instance : HMul Mat4 Mat4 Mat4 := ⟨multiply⟩
instance : HMul Mat4 Vec4 Vec4 := ⟨transformVec4⟩
instance : HMul Mat4 Float Mat4 := ⟨scale⟩

end Mat4

end Linalg
