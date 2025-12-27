/-
  Axis-Aligned Bounding Box primitive.
-/

import Linalg.Vec3

namespace Linalg

/-- Axis-Aligned Bounding Box defined by min and max corners. -/
structure AABB where
  min : Vec3
  max : Vec3
deriving Repr, BEq, Inhabited

namespace AABB

/-- Create AABB from min and max points. -/
def fromMinMax (min max : Vec3) : AABB := ⟨min, max⟩

/-- Create AABB from center and half-extents. -/
def fromCenterExtents (center extents : Vec3) : AABB :=
  ⟨center.sub extents, center.add extents⟩

/-- Create AABB from center and full size. -/
def fromCenterSize (center size : Vec3) : AABB :=
  let halfSize := size.scale 0.5
  ⟨center.sub halfSize, center.add halfSize⟩

/-- Create AABB containing a single point. -/
def fromPoint (p : Vec3) : AABB := ⟨p, p⟩

/-- Get the center point of the AABB. -/
def center (b : AABB) : Vec3 :=
  Vec3.lerp b.min b.max 0.5

/-- Get the half-extents (half-size along each axis). -/
def extents (b : AABB) : Vec3 :=
  (b.max.sub b.min).scale 0.5

/-- Get the full size of the AABB. -/
def size (b : AABB) : Vec3 := b.max.sub b.min

/-- Get the volume of the AABB. -/
def volume (b : AABB) : Float :=
  let s := b.size
  s.x * s.y * s.z

/-- Get the surface area of the AABB. -/
def surfaceArea (b : AABB) : Float :=
  let s := b.size
  2.0 * (s.x * s.y + s.y * s.z + s.z * s.x)

/-- Check if a point is contained within the AABB. -/
def containsPoint (b : AABB) (p : Vec3) : Bool :=
  p.x >= b.min.x && p.x <= b.max.x &&
  p.y >= b.min.y && p.y <= b.max.y &&
  p.z >= b.min.z && p.z <= b.max.z

/-- Check if this AABB fully contains another AABB. -/
def containsAABB (outer inner : AABB) : Bool :=
  outer.containsPoint inner.min && outer.containsPoint inner.max

/-- Expand the AABB to include a point. -/
def expand (b : AABB) (p : Vec3) : AABB :=
  ⟨Vec3.min b.min p, Vec3.max b.max p⟩

/-- Merge two AABBs into one that contains both. -/
def merge (a b : AABB) : AABB :=
  ⟨Vec3.min a.min b.min, Vec3.max a.max b.max⟩

/-- Get the closest point on the AABB to a given point. -/
def closestPoint (b : AABB) (p : Vec3) : Vec3 :=
  Vec3.clamp p b.min b.max

/-- Get the squared distance from a point to the AABB. -/
def distanceSquared (b : AABB) (p : Vec3) : Float :=
  let closest := b.closestPoint p
  p.distanceSquared closest

/-- Get the distance from a point to the AABB. -/
def distance (b : AABB) (p : Vec3) : Float :=
  Float.sqrt (b.distanceSquared p)

/-- Check if two AABBs are approximately equal. -/
def approxEq (a b : AABB) (eps : Float := Float.epsilon) : Bool :=
  a.min.approxEq b.min eps && a.max.approxEq b.max eps

/-- Check if the AABB is valid (min <= max for all components). -/
def isValid (b : AABB) : Bool :=
  b.min.x <= b.max.x && b.min.y <= b.max.y && b.min.z <= b.max.z

end AABB

end Linalg
