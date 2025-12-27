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
