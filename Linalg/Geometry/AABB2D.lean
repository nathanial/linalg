/-
  2D Axis-Aligned Bounding Box primitive.
-/

import Linalg.Vec2

namespace Linalg

/-- 2D Axis-Aligned Bounding Box defined by min and max corners. -/
structure AABB2D where
  min : Vec2
  max : Vec2
deriving Repr, BEq, Inhabited

namespace AABB2D

/-- Create AABB2D from min and max points. -/
def fromMinMax (min max : Vec2) : AABB2D := ⟨min, max⟩

/-- Create AABB2D from center and half-extents. -/
def fromCenterExtents (center extents : Vec2) : AABB2D :=
  ⟨center.sub extents, center.add extents⟩

/-- Create AABB2D from center and full size. -/
def fromCenterSize (center size : Vec2) : AABB2D :=
  let halfSize := size.scale 0.5
  ⟨center.sub halfSize, center.add halfSize⟩

/-- Create AABB2D containing a single point. -/
def fromPoint (p : Vec2) : AABB2D := ⟨p, p⟩

/-- Create AABB2D from an array of points. -/
def fromPoints (points : Array Vec2) : Option AABB2D :=
  if h : points.size > 0 then
    let first := points[0]
    some (points.foldl (fun acc p => ⟨Vec2.min acc.min p, Vec2.max acc.max p⟩) (fromPoint first))
  else
    none

/-- Get the center point of the AABB2D. -/
def center (b : AABB2D) : Vec2 :=
  Vec2.lerp b.min b.max 0.5

/-- Get the half-extents (half-size along each axis). -/
def extents (b : AABB2D) : Vec2 :=
  (b.max.sub b.min).scale 0.5

/-- Get the full size of the AABB2D. -/
def size (b : AABB2D) : Vec2 := b.max.sub b.min

/-- Get the area of the AABB2D. -/
def area (b : AABB2D) : Float :=
  let s := b.size
  s.x * s.y

/-- Get the perimeter of the AABB2D. -/
def perimeter (b : AABB2D) : Float :=
  let s := b.size
  2.0 * (s.x + s.y)

/-- Get the width (x dimension) of the AABB2D. -/
def width (b : AABB2D) : Float := b.max.x - b.min.x

/-- Get the height (y dimension) of the AABB2D. -/
def height (b : AABB2D) : Float := b.max.y - b.min.y

/-- Check if a point is contained within the AABB2D. -/
def containsPoint (b : AABB2D) (p : Vec2) : Bool :=
  p.x >= b.min.x && p.x <= b.max.x &&
  p.y >= b.min.y && p.y <= b.max.y

/-- Check if this AABB2D fully contains another AABB2D. -/
def containsAABB (outer inner : AABB2D) : Bool :=
  outer.containsPoint inner.min && outer.containsPoint inner.max

/-- Check if two AABB2Ds intersect (overlap). -/
def intersects (a b : AABB2D) : Bool :=
  a.min.x <= b.max.x && a.max.x >= b.min.x &&
  a.min.y <= b.max.y && a.max.y >= b.min.y

/-- Compute the intersection of two AABB2Ds (returns none if no overlap). -/
def intersection (a b : AABB2D) : Option AABB2D :=
  if a.intersects b then
    some ⟨Vec2.max a.min b.min, Vec2.min a.max b.max⟩
  else
    none

/-- Expand the AABB2D to include a point. -/
def expand (b : AABB2D) (p : Vec2) : AABB2D :=
  ⟨Vec2.min b.min p, Vec2.max b.max p⟩

/-- Merge two AABB2Ds into one that contains both. -/
def merge (a b : AABB2D) : AABB2D :=
  ⟨Vec2.min a.min b.min, Vec2.max a.max b.max⟩

/-- Get the closest point on the AABB2D to a given point. -/
def closestPoint (b : AABB2D) (p : Vec2) : Vec2 :=
  Vec2.clamp p b.min b.max

/-- Get the squared distance from a point to the AABB2D. -/
def distanceSquared (b : AABB2D) (p : Vec2) : Float :=
  let closest := b.closestPoint p
  p.distanceSquared closest

/-- Get the distance from a point to the AABB2D. -/
def distance (b : AABB2D) (p : Vec2) : Float :=
  Float.sqrt (b.distanceSquared p)

/-- Check if two AABB2Ds are approximately equal. -/
def approxEq (a b : AABB2D) (eps : Float := Float.epsilon) : Bool :=
  a.min.approxEq b.min eps && a.max.approxEq b.max eps

/-- Check if the AABB2D is valid (min <= max for all components). -/
def isValid (b : AABB2D) : Bool :=
  b.min.x <= b.max.x && b.min.y <= b.max.y

/-- Translate the AABB2D by a vector. -/
def translate (b : AABB2D) (v : Vec2) : AABB2D :=
  ⟨b.min.add v, b.max.add v⟩

/-- Scale the AABB2D from its center. -/
def scale (b : AABB2D) (s : Float) : AABB2D :=
  let c := b.center
  let e := b.extents.scale s
  ⟨c.sub e, c.add e⟩

/-- Get the four corner points of the AABB2D. -/
def corners (b : AABB2D) : Array Vec2 :=
  #[b.min,
    Vec2.mk b.max.x b.min.y,
    b.max,
    Vec2.mk b.min.x b.max.y]

/-- Get corner by index (0=bottomLeft, 1=bottomRight, 2=topRight, 3=topLeft). -/
def corner (b : AABB2D) (idx : Fin 4) : Vec2 :=
  match idx with
  | 0 => b.min
  | 1 => Vec2.mk b.max.x b.min.y
  | 2 => b.max
  | 3 => Vec2.mk b.min.x b.max.y

/-- Check if a circle intersects this AABB2D. -/
def intersectsCircle (b : AABB2D) (center : Vec2) (radius : Float) : Bool :=
  b.distanceSquared center <= radius * radius

/-- Subdivide into four quadrants (for quadtree). -/
def subdivide (b : AABB2D) : Array AABB2D :=
  let c := b.center
  #[-- Bottom-left (SW)
    AABB2D.fromMinMax b.min c,
    -- Bottom-right (SE)
    AABB2D.fromMinMax (Vec2.mk c.x b.min.y) (Vec2.mk b.max.x c.y),
    -- Top-left (NW)
    AABB2D.fromMinMax (Vec2.mk b.min.x c.y) (Vec2.mk c.x b.max.y),
    -- Top-right (NE)
    AABB2D.fromMinMax c b.max]

/-- Get the quadrant index for a point (0=SW, 1=SE, 2=NW, 3=NE). -/
def quadrantFor (b : AABB2D) (p : Vec2) : Fin 4 :=
  let c := b.center
  let xBit : Fin 2 := if p.x >= c.x then 1 else 0
  let yBit : Fin 2 := if p.y >= c.y then 1 else 0
  ⟨xBit.val + yBit.val * 2, by
    have hx := xBit.isLt
    have hy := yBit.isLt
    omega⟩

end AABB2D

end Linalg
