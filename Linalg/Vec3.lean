/-
  3D Vector type and operations.
-/

import Linalg.Core
import Linalg.Vec2

namespace Linalg

/-- 3D vector with x, y, and z components. -/
structure Vec3 where
  x : Float
  y : Float
  z : Float
deriving Repr, BEq, Inhabited

namespace Vec3

/-- Zero vector. -/
def zero : Vec3 := ⟨0.0, 0.0, 0.0⟩

/-- One vector (all components 1). -/
def one : Vec3 := ⟨1.0, 1.0, 1.0⟩

/-- Unit vector along X axis. -/
def unitX : Vec3 := ⟨1.0, 0.0, 0.0⟩

/-- Unit vector along Y axis. -/
def unitY : Vec3 := ⟨0.0, 1.0, 0.0⟩

/-- Unit vector along Z axis. -/
def unitZ : Vec3 := ⟨0.0, 0.0, 1.0⟩

/-- Up direction (positive Y). -/
def up : Vec3 := unitY

/-- Down direction (negative Y). -/
def down : Vec3 := ⟨0.0, -1.0, 0.0⟩

/-- Forward direction (negative Z, right-handed). -/
def forward : Vec3 := ⟨0.0, 0.0, -1.0⟩

/-- Back direction (positive Z). -/
def back : Vec3 := unitZ

/-- Right direction (positive X). -/
def right : Vec3 := unitX

/-- Left direction (negative X). -/
def left : Vec3 := ⟨-1.0, 0.0, 0.0⟩

/-- Component-wise addition. -/
@[inline]
def add (a b : Vec3) : Vec3 := ⟨a.x + b.x, a.y + b.y, a.z + b.z⟩

/-- Component-wise subtraction. -/
@[inline]
def sub (a b : Vec3) : Vec3 := ⟨a.x - b.x, a.y - b.y, a.z - b.z⟩

/-- Component-wise multiplication. -/
@[inline]
def mul (a b : Vec3) : Vec3 := ⟨a.x * b.x, a.y * b.y, a.z * b.z⟩

/-- Component-wise division. -/
@[inline]
def div (a b : Vec3) : Vec3 := ⟨a.x / b.x, a.y / b.y, a.z / b.z⟩

/-- Negate a vector. -/
@[inline]
def neg (v : Vec3) : Vec3 := ⟨-v.x, -v.y, -v.z⟩

/-- Scale a vector by a scalar. -/
@[inline]
def scale (v : Vec3) (s : Float) : Vec3 := ⟨v.x * s, v.y * s, v.z * s⟩

/-- Dot product of two vectors. -/
@[inline]
def dot (a b : Vec3) : Float := a.x * b.x + a.y * b.y + a.z * b.z

/-- Cross product of two vectors. -/
def cross (a b : Vec3) : Vec3 :=
  ⟨a.y * b.z - a.z * b.y,
   a.z * b.x - a.x * b.z,
   a.x * b.y - a.y * b.x⟩

/-- Squared length of a vector. -/
@[inline]
def lengthSquared (v : Vec3) : Float := v.dot v

/-- Length (magnitude) of a vector. -/
@[inline]
def length (v : Vec3) : Float := Float.sqrt v.lengthSquared

/-- Normalize a vector to unit length. Returns zero if length is too small. -/
def normalize (v : Vec3) : Vec3 :=
  let len := v.length
  if len > Float.epsilon then v.scale (1.0 / len) else zero

/-- Distance between two points. -/
def distance (a b : Vec3) : Float := (b.sub a).length

/-- Squared distance between two points. -/
def distanceSquared (a b : Vec3) : Float := (b.sub a).lengthSquared

/-- Linear interpolation between two vectors. -/
def lerp (a b : Vec3) (t : Float) : Vec3 :=
  ⟨Float.lerp a.x b.x t, Float.lerp a.y b.y t, Float.lerp a.z b.z t⟩

/-- Project vector v onto vector onto. -/
def project (v onto : Vec3) : Vec3 :=
  let d := onto.dot onto
  if d > Float.epsilon then onto.scale (v.dot onto / d) else zero

/-- Reflect vector v off a surface with the given normal. -/
def reflect (v normal : Vec3) : Vec3 :=
  v.sub (normal.scale (2.0 * v.dot normal))

/-- Angle between two vectors in radians. -/
def angle (a b : Vec3) : Float :=
  let d := a.dot b
  let la := a.length
  let lb := b.length
  if la > Float.epsilon && lb > Float.epsilon then
    Float.acos (Float.clamp (d / (la * lb)) (-1.0) 1.0)
  else 0.0

/-- Convert to Vec2 (drop z component). -/
def toVec2 (v : Vec3) : Vec2 := ⟨v.x, v.y⟩

/-- Create from Vec2 with specified z. -/
def fromVec2 (v : Vec2) (z : Float := 0.0) : Vec3 := ⟨v.x, v.y, z⟩

/-- Check if two vectors are approximately equal. -/
def approxEq (a b : Vec3) (eps : Float := Float.epsilon) : Bool :=
  Float.approxEq a.x b.x eps && Float.approxEq a.y b.y eps && Float.approxEq a.z b.z eps

/-- Minimum of each component. -/
def min (a b : Vec3) : Vec3 := ⟨Float.min a.x b.x, Float.min a.y b.y, Float.min a.z b.z⟩

/-- Maximum of each component. -/
def max (a b : Vec3) : Vec3 := ⟨Float.max a.x b.x, Float.max a.y b.y, Float.max a.z b.z⟩

/-- Clamp each component to a range. -/
def clamp (v lo hi : Vec3) : Vec3 :=
  ⟨Float.clamp v.x lo.x hi.x, Float.clamp v.y lo.y hi.y, Float.clamp v.z lo.z hi.z⟩

/-- Absolute value of each component. -/
def abs (v : Vec3) : Vec3 := ⟨Float.abs' v.x, Float.abs' v.y, Float.abs' v.z⟩

instance : Add Vec3 := ⟨add⟩
instance : Sub Vec3 := ⟨sub⟩
instance : Neg Vec3 := ⟨neg⟩
instance : HMul Vec3 Float Vec3 := ⟨scale⟩
instance : HMul Float Vec3 Vec3 := ⟨fun s v => scale v s⟩
instance : HDiv Vec3 Float Vec3 := ⟨fun v s => scale v (1.0 / s)⟩

/-- Coerce Vec2 to Vec3 with z=0. -/
instance : Coe Vec2 Vec3 := ⟨fromVec2⟩

/-- Coerce a Float to a uniform Vec3. -/
instance : Coe Float Vec3 := ⟨fun f => ⟨f, f, f⟩⟩

end Vec3

end Linalg
