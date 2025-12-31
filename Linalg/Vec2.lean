/-
  2D Vector type and operations.
-/

import Linalg.Core

namespace Linalg

/-- 2D vector with x and y components. -/
structure Vec2 where
  x : Float
  y : Float
deriving Repr, BEq, Inhabited

namespace Vec2

/-- Zero vector. -/
def zero : Vec2 := ⟨0.0, 0.0⟩

/-- One vector (all components 1). -/
def one : Vec2 := ⟨1.0, 1.0⟩

/-- Unit vector along X axis. -/
def unitX : Vec2 := ⟨1.0, 0.0⟩

/-- Unit vector along Y axis. -/
def unitY : Vec2 := ⟨0.0, 1.0⟩

/-- Component-wise addition. -/
@[inline]
def add (a b : Vec2) : Vec2 := ⟨a.x + b.x, a.y + b.y⟩

/-- Component-wise subtraction. -/
@[inline]
def sub (a b : Vec2) : Vec2 := ⟨a.x - b.x, a.y - b.y⟩

/-- Component-wise multiplication. -/
@[inline]
def mul (a b : Vec2) : Vec2 := ⟨a.x * b.x, a.y * b.y⟩

/-- Component-wise division. -/
@[inline]
def div (a b : Vec2) : Vec2 := ⟨a.x / b.x, a.y / b.y⟩

/-- Negate a vector. -/
@[inline]
def neg (v : Vec2) : Vec2 := ⟨-v.x, -v.y⟩

/-- Scale a vector by a scalar. -/
@[inline]
def scale (v : Vec2) (s : Float) : Vec2 := ⟨v.x * s, v.y * s⟩

/-- Dot product of two vectors. -/
@[inline]
def dot (a b : Vec2) : Float := a.x * b.x + a.y * b.y

/-- Squared length of a vector. -/
@[inline]
def lengthSquared (v : Vec2) : Float := v.dot v

/-- Length (magnitude) of a vector. -/
@[inline]
def length (v : Vec2) : Float := Float.sqrt v.lengthSquared

/-- Normalize a vector to unit length. Returns zero if length is too small. -/
def normalize (v : Vec2) : Vec2 :=
  let len := v.length
  if len > Float.epsilon then v.scale (1.0 / len) else zero

/-- Distance between two points. -/
def distance (a b : Vec2) : Float := (b.sub a).length

/-- Squared distance between two points. -/
def distanceSquared (a b : Vec2) : Float := (b.sub a).lengthSquared

/-- Linear interpolation between two vectors. -/
def lerp (a b : Vec2) (t : Float) : Vec2 :=
  ⟨Float.lerp a.x b.x t, Float.lerp a.y b.y t⟩

/-- Project vector v onto vector onto. -/
def project (v onto : Vec2) : Vec2 :=
  let d := onto.dot onto
  if d > Float.epsilon then onto.scale (v.dot onto / d) else zero

/-- Reflect vector v off a surface with the given normal. -/
def reflect (v normal : Vec2) : Vec2 :=
  v.sub (normal.scale (2.0 * v.dot normal))

/-- Angle between two vectors in radians. -/
def angle (a b : Vec2) : Float :=
  let d := a.dot b
  let la := a.length
  let lb := b.length
  if la > Float.epsilon && lb > Float.epsilon then
    Float.acos (Float.clamp (d / (la * lb)) (-1.0) 1.0)
  else 0.0

/-- Get perpendicular vector (90 degrees counter-clockwise). -/
def perpendicular (v : Vec2) : Vec2 := ⟨-v.y, v.x⟩

/-- 2D cross product (returns the z component of the 3D cross product). -/
def cross (a b : Vec2) : Float := a.x * b.y - a.y * b.x

/-- Check if two vectors are approximately equal. -/
def approxEq (a b : Vec2) (eps : Float := Float.epsilon) : Bool :=
  Float.approxEq a.x b.x eps && Float.approxEq a.y b.y eps

/-- Minimum of each component. -/
def min (a b : Vec2) : Vec2 := ⟨Float.min a.x b.x, Float.min a.y b.y⟩

/-- Maximum of each component. -/
def max (a b : Vec2) : Vec2 := ⟨Float.max a.x b.x, Float.max a.y b.y⟩

/-- Clamp each component to a range. -/
def clamp (v lo hi : Vec2) : Vec2 :=
  ⟨Float.clamp v.x lo.x hi.x, Float.clamp v.y lo.y hi.y⟩

/-- Absolute value of each component. -/
def abs (v : Vec2) : Vec2 := ⟨Float.abs' v.x, Float.abs' v.y⟩

instance : Add Vec2 := ⟨add⟩
instance : Sub Vec2 := ⟨sub⟩
instance : Neg Vec2 := ⟨neg⟩
instance : HMul Vec2 Float Vec2 := ⟨scale⟩
instance : HMul Float Vec2 Vec2 := ⟨fun s v => scale v s⟩
instance : HDiv Vec2 Float Vec2 := ⟨fun v s => scale v (1.0 / s)⟩

/-- Coerce a Float to a uniform Vec2. -/
instance : Coe Float Vec2 := ⟨fun f => ⟨f, f⟩⟩

end Vec2

end Linalg
