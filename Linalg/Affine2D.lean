/-
  Affine2D - 2D Affine Transform Matrix

  A 2x3 matrix representing 2D affine transformations including:
  - Translation
  - Rotation
  - Scaling (uniform and non-uniform)
  - Shearing

  Matrix layout (column-major storage):
  [a  c  tx]   transforms point (x,y) to:
  [b  d  ty]   (a*x + c*y + tx, b*x + d*y + ty)

  Storage: #[a, b, c, d, tx, ty] (6 elements, column-major)
-/

import Linalg.Core
import Linalg.Vec2
import Linalg.Mat2

namespace Linalg

/-- 2D affine transform matrix (2x3, column-major). -/
structure Affine2D where
  data : Array Float
deriving Inhabited, Repr

namespace Affine2D

/-- Identity transform. -/
def identity : Affine2D := { data := #[1, 0, 0, 1, 0, 0] }

/-- Create from column vectors. -/
def fromColumns (c0 c1 t : Vec2) : Affine2D :=
  { data := #[c0.x, c0.y, c1.x, c1.y, t.x, t.y] }

/-- Create from individual components. -/
def mk' (a b c d tx ty : Float) : Affine2D :=
  { data := #[a, b, c, d, tx, ty] }

/-- Get element at (row, col). Row 0-1, Col 0-2. -/
@[inline]
def get (m : Affine2D) (row col : Nat) : Float :=
  m.data.getD (col * 2 + row) 0.0

/-- Set element at (row, col). -/
def set (m : Affine2D) (row col : Nat) (v : Float) : Affine2D :=
  { data := m.data.set! (col * 2 + row) v }

/-- Get column as Vec2. -/
def getColumn (m : Affine2D) (col : Nat) : Vec2 :=
  ⟨m.get 0 col, m.get 1 col⟩

/-- First basis vector (x-axis after transform). -/
def basisX (m : Affine2D) : Vec2 := m.getColumn 0

/-- Second basis vector (y-axis after transform). -/
def basisY (m : Affine2D) : Vec2 := m.getColumn 1

/-- Translation component. -/
def translation (m : Affine2D) : Vec2 := m.getColumn 2

-- ============================================================================
-- Constructors
-- ============================================================================

/-- Pure translation transform. -/
def translationMat (tx ty : Float) : Affine2D :=
  { data := #[1, 0, 0, 1, tx, ty] }

/-- Pure translation from vector. -/
def translationV (v : Vec2) : Affine2D :=
  translationMat v.x v.y

/-- Pure rotation transform. -/
def rotation (angle : Float) : Affine2D :=
  let c := Float.cos angle
  let s := Float.sin angle
  { data := #[c, s, -s, c, 0, 0] }

/-- Pure scaling transform. -/
def scaling (sx sy : Float) : Affine2D :=
  { data := #[sx, 0, 0, sy, 0, 0] }

/-- Uniform scaling transform. -/
def uniformScaling (s : Float) : Affine2D :=
  scaling s s

/-- Shear transform. sx shears along x, sy shears along y. -/
def shear (sx sy : Float) : Affine2D :=
  { data := #[1, sy, sx, 1, 0, 0] }

/-- Create rotation around a pivot point. -/
def rotationAround (pivot : Vec2) (angle : Float) : Affine2D :=
  -- Translate to origin, rotate, translate back
  let c := Float.cos angle
  let s := Float.sin angle
  let tx := pivot.x - c * pivot.x + s * pivot.y
  let ty := pivot.y - s * pivot.x - c * pivot.y
  { data := #[c, s, -s, c, tx, ty] }

/-- Create scaling around a pivot point. -/
def scalingAround (pivot : Vec2) (sx sy : Float) : Affine2D :=
  let tx := pivot.x * (1 - sx)
  let ty := pivot.y * (1 - sy)
  { data := #[sx, 0, 0, sy, tx, ty] }

-- ============================================================================
-- Transform Operations
-- ============================================================================

/-- Transform a point (applies translation). -/
@[inline]
def transformPoint (m : Affine2D) (p : Vec2) : Vec2 :=
  ⟨m.get 0 0 * p.x + m.get 0 1 * p.y + m.get 0 2,
   m.get 1 0 * p.x + m.get 1 1 * p.y + m.get 1 2⟩

/-- Transform a vector (no translation). -/
@[inline]
def transformVector (m : Affine2D) (v : Vec2) : Vec2 :=
  ⟨m.get 0 0 * v.x + m.get 0 1 * v.y,
   m.get 1 0 * v.x + m.get 1 1 * v.y⟩

/-- Compose two affine transforms (this * other). -/
def compose (a b : Affine2D) : Affine2D :=
  -- [a0 a2 a4]   [b0 b2 b4]   [a0*b0+a2*b1  a0*b2+a2*b3  a0*b4+a2*b5+a4]
  -- [a1 a3 a5] * [b1 b3 b5] = [a1*b0+a3*b1  a1*b2+a3*b3  a1*b4+a3*b5+a5]
  --              [0  0  1 ]
  let a0 := a.get 0 0; let a1 := a.get 1 0
  let a2 := a.get 0 1; let a3 := a.get 1 1
  let a4 := a.get 0 2; let a5 := a.get 1 2
  let b0 := b.get 0 0; let b1 := b.get 1 0
  let b2 := b.get 0 1; let b3 := b.get 1 1
  let b4 := b.get 0 2; let b5 := b.get 1 2
  { data := #[
    a0*b0 + a2*b1, a1*b0 + a3*b1,
    a0*b2 + a2*b3, a1*b2 + a3*b3,
    a0*b4 + a2*b5 + a4, a1*b4 + a3*b5 + a5
  ]}

/-- Determinant of the linear part. -/
def determinant (m : Affine2D) : Float :=
  m.get 0 0 * m.get 1 1 - m.get 0 1 * m.get 1 0

/-- Inverse transform. Returns None if singular. -/
def inverse (m : Affine2D) : Option Affine2D :=
  let det := m.determinant
  if Float.abs' det < Float.epsilon then none
  else
    let invDet := 1.0 / det
    let a := m.get 0 0; let b := m.get 1 0
    let c := m.get 0 1; let d := m.get 1 1
    let tx := m.get 0 2; let ty := m.get 1 2
    -- Inverse of linear part
    let a' := d * invDet
    let b' := -b * invDet
    let c' := -c * invDet
    let d' := a * invDet
    -- Inverse translation: -M^(-1) * t
    let tx' := -(a' * tx + c' * ty)
    let ty' := -(b' * tx + d' * ty)
    some { data := #[a', b', c', d', tx', ty'] }

-- ============================================================================
-- Decomposition
-- ============================================================================

/-- Extract the linear part as Mat2. -/
def toMat2 (m : Affine2D) : Mat2 :=
  { data := #[m.get 0 0, m.get 1 0, m.get 0 1, m.get 1 1] }

/-- Extract rotation angle (assumes no shear). -/
def extractRotation (m : Affine2D) : Float :=
  Float.atan2 (m.get 1 0) (m.get 0 0)

/-- Extract scale factors (assumes no shear). -/
def extractScale (m : Affine2D) : Vec2 :=
  let sx := (m.basisX).length
  let sy := (m.basisY).length
  ⟨sx, sy⟩

/-- Check if this is approximately the identity transform. -/
def isIdentity (m : Affine2D) (epsilon : Float := Float.epsilon) : Bool :=
  Float.approxEq (m.get 0 0) 1.0 epsilon &&
  Float.approxEq (m.get 1 0) 0.0 epsilon &&
  Float.approxEq (m.get 0 1) 0.0 epsilon &&
  Float.approxEq (m.get 1 1) 1.0 epsilon &&
  Float.approxEq (m.get 0 2) 0.0 epsilon &&
  Float.approxEq (m.get 1 2) 0.0 epsilon

-- ============================================================================
-- Additional Operations
-- ============================================================================

/-- Translate this transform. -/
def translate (m : Affine2D) (tx ty : Float) : Affine2D :=
  m.compose (translationMat tx ty)

/-- Pre-translate this transform. -/
def preTranslate (m : Affine2D) (tx ty : Float) : Affine2D :=
  (translationMat tx ty).compose m

/-- Rotate this transform. -/
def rotate (m : Affine2D) (angle : Float) : Affine2D :=
  m.compose (rotation angle)

/-- Scale this transform. -/
def scale (m : Affine2D) (sx sy : Float) : Affine2D :=
  m.compose (scaling sx sy)

/-- Lerp between two transforms (component-wise). -/
def lerp (a b : Affine2D) (t : Float) : Affine2D :=
  let result := Id.run do
    let mut arr : Array Float := #[]
    for i in [:6] do
      let va := a.data.getD i 0.0
      let vb := b.data.getD i 0.0
      arr := arr.push (Float.lerp va vb t)
    return arr
  { data := result }

-- ============================================================================
-- Typeclass Instances
-- ============================================================================

instance : HMul Affine2D Affine2D Affine2D := ⟨compose⟩
instance : HMul Affine2D Vec2 Vec2 := ⟨transformPoint⟩

instance : BEq Affine2D where
  beq a b := a.data == b.data

end Affine2D

end Linalg
