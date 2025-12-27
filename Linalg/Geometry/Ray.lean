/-
  Ray primitive for raycasting and intersection tests.
-/

import Linalg.Vec3

namespace Linalg

/-- A ray with origin and direction. -/
structure Ray where
  origin : Vec3
  direction : Vec3  -- Should be normalized
deriving Repr, BEq, Inhabited

namespace Ray

/-- Create a ray with automatic direction normalization. -/
def mk' (origin direction : Vec3) : Ray :=
  ⟨origin, direction.normalize⟩

/-- Get point along the ray at parameter t. -/
def pointAt (r : Ray) (t : Float) : Vec3 :=
  r.origin.add (r.direction.scale t)

/-- Check if two rays are approximately equal. -/
def approxEq (a b : Ray) (eps : Float := Float.epsilon) : Bool :=
  a.origin.approxEq b.origin eps && a.direction.approxEq b.direction eps

end Ray

end Linalg
