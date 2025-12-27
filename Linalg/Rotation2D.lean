/-
  Rotation2D - 2D Rotation Type

  A simple angle-based 2D rotation type. More efficient than a full matrix
  for pure rotations, and provides a clear semantic type for rotation values.
-/

import Linalg.Core
import Linalg.Vec2
import Linalg.Mat2
import Linalg.Affine2D

namespace Linalg

/-- 2D rotation represented as an angle in radians. -/
structure Rotation2D where
  angle : Float
deriving Inhabited, Repr, BEq

namespace Rotation2D

-- ============================================================================
-- Constructors
-- ============================================================================

/-- Identity (zero) rotation. -/
def identity : Rotation2D := ⟨0.0⟩

/-- Create from radians. -/
def fromRadians (rad : Float) : Rotation2D := ⟨rad⟩

/-- Create from degrees. -/
def fromDegrees (deg : Float) : Rotation2D := ⟨Float.toRadians deg⟩

/-- Create rotation to point in a given direction. -/
def fromDirection (dir : Vec2) : Rotation2D :=
  ⟨Float.atan2 dir.y dir.x⟩

/-- Create rotation from one direction to another. -/
def fromTo (fromDir toDir : Vec2) : Rotation2D :=
  let a1 := Float.atan2 fromDir.y fromDir.x
  let a2 := Float.atan2 toDir.y toDir.x
  ⟨a2 - a1⟩

-- ============================================================================
-- Accessors
-- ============================================================================

/-- Get angle in degrees. -/
def degrees (r : Rotation2D) : Float :=
  Float.toDegrees r.angle

/-- Get cosine of the angle. -/
@[inline]
def cos (r : Rotation2D) : Float :=
  Float.cos r.angle

/-- Get sine of the angle. -/
@[inline]
def sin (r : Rotation2D) : Float :=
  Float.sin r.angle

/-- Get unit vector in the direction of this rotation. -/
def direction (r : Rotation2D) : Vec2 :=
  ⟨r.cos, r.sin⟩

-- ============================================================================
-- Operations
-- ============================================================================

/-- Rotate a vector. -/
@[inline]
def rotate (r : Rotation2D) (v : Vec2) : Vec2 :=
  let c := r.cos
  let s := r.sin
  ⟨c * v.x - s * v.y, s * v.x + c * v.y⟩

/-- Inverse rotation (negate angle). -/
def inverse (r : Rotation2D) : Rotation2D :=
  ⟨-r.angle⟩

/-- Compose two rotations (add angles). -/
def compose (a b : Rotation2D) : Rotation2D :=
  ⟨a.angle + b.angle⟩

/-- Normalize angle to [-pi, pi]. -/
def normalize (r : Rotation2D) : Rotation2D :=
  let a := r.angle
  let twoPi := 2.0 * Float.pi
  let normalized := a - twoPi * Float.floor ((a + Float.pi) / twoPi)
  ⟨normalized⟩

/-- Lerp between two rotations. -/
def lerp (a b : Rotation2D) (t : Float) : Rotation2D :=
  -- Handle wraparound by finding shortest path
  let diff := (b.angle - a.angle)
  let twoPi := 2.0 * Float.pi
  -- Normalize diff to [-pi, pi]
  let normalizedDiff := diff - twoPi * Float.floor ((diff + Float.pi) / twoPi)
  ⟨a.angle + normalizedDiff * t⟩

/-- Spherical lerp (same as lerp for 2D). -/
def slerp (a b : Rotation2D) (t : Float) : Rotation2D :=
  lerp a b t

-- ============================================================================
-- Conversion
-- ============================================================================

/-- Convert to 2x2 rotation matrix. -/
def toMat2 (r : Rotation2D) : Mat2 :=
  let c := r.cos
  let s := r.sin
  { data := #[c, s, -s, c] }

/-- Convert to 2D affine transform. -/
def toAffine2D (r : Rotation2D) : Affine2D :=
  Affine2D.rotation r.angle

-- ============================================================================
-- Approximate Equality
-- ============================================================================

/-- Check if two rotations are approximately equal. -/
def approxEq (a b : Rotation2D) (epsilon : Float := Float.epsilon) : Bool :=
  let diff := (a.normalize).angle - (b.normalize).angle
  Float.abs' diff < epsilon

-- ============================================================================
-- Typeclass Instances
-- ============================================================================

instance : HMul Rotation2D Rotation2D Rotation2D := ⟨compose⟩
instance : HMul Rotation2D Vec2 Vec2 := ⟨rotate⟩

end Rotation2D

end Linalg
