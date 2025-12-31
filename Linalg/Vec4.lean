/-
  4D Vector type and operations (homogeneous coordinates).
-/

import Linalg.Core
import Linalg.Vec3

namespace Linalg

/-- 4D vector with x, y, z, and w components. -/
structure Vec4 where
  x : Float
  y : Float
  z : Float
  w : Float
deriving Repr, BEq, Inhabited

namespace Vec4

/-- Zero vector. -/
def zero : Vec4 := ⟨0.0, 0.0, 0.0, 0.0⟩

/-- One vector (all components 1). -/
def one : Vec4 := ⟨1.0, 1.0, 1.0, 1.0⟩

/-- Unit vector along X axis. -/
def unitX : Vec4 := ⟨1.0, 0.0, 0.0, 0.0⟩

/-- Unit vector along Y axis. -/
def unitY : Vec4 := ⟨0.0, 1.0, 0.0, 0.0⟩

/-- Unit vector along Z axis. -/
def unitZ : Vec4 := ⟨0.0, 0.0, 1.0, 0.0⟩

/-- Unit vector along W axis. -/
def unitW : Vec4 := ⟨0.0, 0.0, 0.0, 1.0⟩

/-- Component-wise addition. -/
@[inline]
def add (a b : Vec4) : Vec4 := ⟨a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w⟩

/-- Component-wise subtraction. -/
@[inline]
def sub (a b : Vec4) : Vec4 := ⟨a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w⟩

/-- Component-wise multiplication. -/
@[inline]
def mul (a b : Vec4) : Vec4 := ⟨a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w⟩

/-- Component-wise division. -/
@[inline]
def div (a b : Vec4) : Vec4 := ⟨a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w⟩

/-- Negate a vector. -/
@[inline]
def neg (v : Vec4) : Vec4 := ⟨-v.x, -v.y, -v.z, -v.w⟩

/-- Scale a vector by a scalar. -/
@[inline]
def scale (v : Vec4) (s : Float) : Vec4 := ⟨v.x * s, v.y * s, v.z * s, v.w * s⟩

/-- Dot product of two vectors. -/
@[inline]
def dot (a b : Vec4) : Float := a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w

/-- Squared length of a vector. -/
@[inline]
def lengthSquared (v : Vec4) : Float := v.dot v

/-- Length (magnitude) of a vector. -/
@[inline]
def length (v : Vec4) : Float := Float.sqrt v.lengthSquared

/-- Normalize a vector to unit length. Returns zero if length is too small. -/
def normalize (v : Vec4) : Vec4 :=
  let len := v.length
  if len > Float.epsilon then v.scale (1.0 / len) else zero

/-- Linear interpolation between two vectors. -/
def lerp (a b : Vec4) (t : Float) : Vec4 :=
  ⟨Float.lerp a.x b.x t, Float.lerp a.y b.y t, Float.lerp a.z b.z t, Float.lerp a.w b.w t⟩

/-- Create from a point in 3D space (w=1 for homogeneous coordinates). -/
def fromPoint (v : Vec3) : Vec4 := ⟨v.x, v.y, v.z, 1.0⟩

/-- Create from a direction in 3D space (w=0 for homogeneous coordinates). -/
def fromDirection (v : Vec3) : Vec4 := ⟨v.x, v.y, v.z, 0.0⟩

/-- Convert to Vec3 (drop w component). -/
def toVec3 (v : Vec4) : Vec3 := ⟨v.x, v.y, v.z⟩

/-- Convert to Vec3 with perspective divide (divide xyz by w). -/
def toVec3Normalized (v : Vec4) : Vec3 :=
  if Float.abs' v.w > Float.epsilon then ⟨v.x / v.w, v.y / v.w, v.z / v.w⟩
  else ⟨v.x, v.y, v.z⟩

/-- Check if two vectors are approximately equal. -/
def approxEq (a b : Vec4) (eps : Float := Float.epsilon) : Bool :=
  Float.approxEq a.x b.x eps && Float.approxEq a.y b.y eps &&
  Float.approxEq a.z b.z eps && Float.approxEq a.w b.w eps

/-- Minimum of each component. -/
def min (a b : Vec4) : Vec4 :=
  ⟨Float.min a.x b.x, Float.min a.y b.y, Float.min a.z b.z, Float.min a.w b.w⟩

/-- Maximum of each component. -/
def max (a b : Vec4) : Vec4 :=
  ⟨Float.max a.x b.x, Float.max a.y b.y, Float.max a.z b.z, Float.max a.w b.w⟩

/-- Clamp each component to a range. -/
def clamp (v lo hi : Vec4) : Vec4 :=
  ⟨Float.clamp v.x lo.x hi.x, Float.clamp v.y lo.y hi.y,
   Float.clamp v.z lo.z hi.z, Float.clamp v.w lo.w hi.w⟩

/-- Absolute value of each component. -/
def abs (v : Vec4) : Vec4 := ⟨Float.abs' v.x, Float.abs' v.y, Float.abs' v.z, Float.abs' v.w⟩

instance : Add Vec4 := ⟨add⟩
instance : Sub Vec4 := ⟨sub⟩
instance : Neg Vec4 := ⟨neg⟩
instance : HMul Vec4 Float Vec4 := ⟨scale⟩
instance : HMul Float Vec4 Vec4 := ⟨fun s v => scale v s⟩
instance : HDiv Vec4 Float Vec4 := ⟨fun v s => scale v (1.0 / s)⟩

/-- Coerce Vec3 to Vec4 with w=0 (direction). -/
instance : Coe Vec3 Vec4 := ⟨fromDirection⟩

/-- Coerce a Float to a uniform Vec4. -/
instance : Coe Float Vec4 := ⟨fun f => ⟨f, f, f, f⟩⟩

end Vec4

end Linalg
