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
