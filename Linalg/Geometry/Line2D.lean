/-
  Line2D - 2D Line and Line Segment Primitives

  Line2D: An infinite line defined by a point and direction.
  Segment2D: A finite line segment defined by two endpoints.
-/

import Linalg.Core
import Linalg.Vec2

namespace Linalg

-- ============================================================================
-- Line2D - Infinite Line
-- ============================================================================

/-- An infinite 2D line defined by a point on the line and a direction. -/
structure Line2D where
  point : Vec2
  direction : Vec2
deriving Repr, BEq, Inhabited

namespace Line2D

-- ============================================================================
-- Constructors
-- ============================================================================

/-- Create a line from a point and direction (normalizes direction). -/
def mk' (point direction : Vec2) : Line2D :=
  { point := point, direction := direction.normalize }

/-- Create a line from two points. -/
def fromPoints (p1 p2 : Vec2) : Line2D :=
  mk' p1 (p2 - p1)

/-- Create a line from point-slope form (y = mx + b). -/
def fromSlope (slope : Float) (yIntercept : Float) : Line2D :=
  let point := Vec2.mk 0.0 yIntercept
  let direction := Vec2.mk 1.0 slope
  mk' point direction

/-- Create horizontal line at given y. -/
def horizontal (y : Float) : Line2D :=
  { point := Vec2.mk 0.0 y, direction := Vec2.mk 1.0 0.0 }

/-- Create vertical line at given x. -/
def vertical (x : Float) : Line2D :=
  { point := Vec2.mk x 0.0, direction := Vec2.mk 0.0 1.0 }

-- ============================================================================
-- Properties
-- ============================================================================

/-- Get the normal vector (perpendicular to direction). -/
def normal (l : Line2D) : Vec2 := l.direction.perpendicular

/-- Check if the line is approximately horizontal. -/
def isHorizontal (l : Line2D) (tolerance : Float := Float.epsilon) : Bool :=
  Float.abs' l.direction.y < tolerance

/-- Check if the line is approximately vertical. -/
def isVertical (l : Line2D) (tolerance : Float := Float.epsilon) : Bool :=
  Float.abs' l.direction.x < tolerance

/-- Get the slope of the line (returns none for vertical lines). -/
def slope (l : Line2D) : Option Float :=
  if Float.abs' l.direction.x < Float.epsilon then none
  else some (l.direction.y / l.direction.x)

-- ============================================================================
-- Geometric Queries
-- ============================================================================

/-- Get a point on the line at parameter t. -/
def pointAt (l : Line2D) (t : Float) : Vec2 :=
  l.point + l.direction.scale t

/-- Signed distance from a point to the line (positive on left side). -/
def signedDistance (l : Line2D) (p : Vec2) : Float :=
  let toPoint := p - l.point
  toPoint.cross l.direction

/-- Absolute distance from a point to the line. -/
def distance (l : Line2D) (p : Vec2) : Float :=
  Float.abs' (l.signedDistance p)

/-- Project a point onto the line, returning the closest point. -/
def closestPoint (l : Line2D) (p : Vec2) : Vec2 :=
  let toPoint := p - l.point
  let t := toPoint.dot l.direction
  l.point + l.direction.scale t

/-- Get the parameter t for the closest point to p. -/
def closestPointParameter (l : Line2D) (p : Vec2) : Float :=
  (p - l.point).dot l.direction

/-- Check which side of the line a point is on (-1 = right, 0 = on, 1 = left). -/
def side (l : Line2D) (p : Vec2) : Int :=
  let d := l.signedDistance p
  if d > Float.epsilon then 1
  else if d < -Float.epsilon then -1
  else 0

-- ============================================================================
-- Line-Line Intersection
-- ============================================================================

/-- Check if two lines are parallel. -/
def isParallel (l1 l2 : Line2D) (tolerance : Float := Float.epsilon) : Bool :=
  let cross := l1.direction.cross l2.direction
  Float.abs' cross < tolerance

/-- Find intersection point of two lines (None if parallel). -/
def intersection (l1 l2 : Line2D) : Option Vec2 :=
  let cross := l1.direction.cross l2.direction
  if Float.abs' cross < Float.epsilon then none
  else
    let diff := l2.point - l1.point
    let t := diff.cross l2.direction / cross
    some (l1.pointAt t)

end Line2D

-- ============================================================================
-- Segment2D - Finite Line Segment
-- ============================================================================

/-- A finite 2D line segment defined by two endpoints. -/
structure Segment2D where
  a : Vec2  -- Start point
  b : Vec2  -- End point
deriving Repr, BEq, Inhabited

namespace Segment2D

-- ============================================================================
-- Constructors
-- ============================================================================

/-- Create a segment from coordinates. -/
def mk' (ax ay bx by_ : Float) : Segment2D :=
  { a := Vec2.mk ax ay, b := Vec2.mk bx by_ }

/-- Create a segment from a point, direction, and length. -/
def fromPointDirection (point direction : Vec2) (len : Float) : Segment2D :=
  { a := point, b := point + direction.normalize.scale len }

-- ============================================================================
-- Properties
-- ============================================================================

/-- Length of the segment. -/
def length (s : Segment2D) : Float := (s.b - s.a).length

/-- Squared length of the segment. -/
def lengthSquared (s : Segment2D) : Float := (s.b - s.a).lengthSquared

/-- Direction vector (not normalized). -/
def delta (s : Segment2D) : Vec2 := s.b - s.a

/-- Normalized direction vector. -/
def direction (s : Segment2D) : Vec2 := s.delta.normalize

/-- Midpoint of the segment. -/
def midpoint (s : Segment2D) : Vec2 := (s.a + s.b).scale 0.5

/-- Normal vector (perpendicular to segment, pointing left). -/
def normal (s : Segment2D) : Vec2 := s.direction.perpendicular

/-- Convert to infinite line. -/
def toLine (s : Segment2D) : Line2D :=
  { point := s.a, direction := s.direction }

-- ============================================================================
-- Point Queries
-- ============================================================================

/-- Get a point on the segment at parameter t (0 = a, 1 = b). -/
def pointAt (s : Segment2D) (t : Float) : Vec2 :=
  Vec2.lerp s.a s.b t

/-- Find the closest point on the segment to a given point. -/
def closestPoint (s : Segment2D) (p : Vec2) : Vec2 :=
  let ab := s.b - s.a
  let ap := p - s.a
  let lenSq := ab.lengthSquared
  if lenSq < Float.epsilon then s.a
  else
    let t := Float.clamp (ap.dot ab / lenSq) 0.0 1.0
    s.a + ab.scale t

/-- Get the parameter t (clamped to [0,1]) for the closest point. -/
def closestPointParameter (s : Segment2D) (p : Vec2) : Float :=
  let ab := s.b - s.a
  let ap := p - s.a
  let lenSq := ab.lengthSquared
  if lenSq < Float.epsilon then 0.0
  else Float.clamp (ap.dot ab / lenSq) 0.0 1.0

/-- Distance from a point to the segment. -/
def distance (s : Segment2D) (p : Vec2) : Float :=
  (p - s.closestPoint p).length

/-- Squared distance from a point to the segment. -/
def distanceSquared (s : Segment2D) (p : Vec2) : Float :=
  (p - s.closestPoint p).lengthSquared

/-- Check if a point lies on the segment. -/
def containsPoint (s : Segment2D) (p : Vec2) (tolerance : Float := Float.epsilon) : Bool :=
  s.distance p < tolerance

-- ============================================================================
-- Segment-Segment Intersection
-- ============================================================================

/-- Check if two segments intersect. -/
def intersects (s1 s2 : Segment2D) : Bool :=
  let d1 := s1.delta
  let d2 := s2.delta
  let cross := d1.cross d2

  if Float.abs' cross < Float.epsilon then
    -- Parallel segments - check if collinear and overlapping
    let d := s2.a - s1.a
    if Float.abs' (d.cross d1) < Float.epsilon then
      -- Collinear - check overlap
      let t0 := d.dot d1 / d1.lengthSquared
      let t1 := t0 + d2.dot d1 / d1.lengthSquared
      let tMin := Float.min t0 t1
      let tMax := Float.max t0 t1
      tMax >= 0.0 && tMin <= 1.0
    else false
  else
    let d := s2.a - s1.a
    let t := d.cross d2 / cross
    let u := d.cross d1 / cross
    t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0

/-- Find intersection point of two segments (None if no intersection). -/
def intersection (s1 s2 : Segment2D) : Option Vec2 :=
  let d1 := s1.delta
  let d2 := s2.delta
  let cross := d1.cross d2

  if Float.abs' cross < Float.epsilon then none
  else
    let d := s2.a - s1.a
    let t := d.cross d2 / cross
    let u := d.cross d1 / cross
    if t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0 then
      some (s1.pointAt t)
    else none

-- ============================================================================
-- Segment-Line Intersection
-- ============================================================================

/-- Find intersection of segment with infinite line. -/
def intersectLine (s : Segment2D) (l : Line2D) : Option Vec2 :=
  let d := s.delta
  let cross := d.cross l.direction
  if Float.abs' cross < Float.epsilon then none
  else
    let diff := l.point - s.a
    let t := diff.cross l.direction / cross
    if t >= 0.0 && t <= 1.0 then
      some (s.pointAt t)
    else none

-- ============================================================================
-- Utilities
-- ============================================================================

/-- Reverse the segment (swap endpoints). -/
def reverse (s : Segment2D) : Segment2D :=
  { a := s.b, b := s.a }

/-- Extend the segment by a given amount on both ends. -/
def extend (s : Segment2D) (amount : Float) : Segment2D :=
  let dir := s.direction
  { a := s.a - dir.scale amount, b := s.b + dir.scale amount }

/-- Shrink the segment by a given amount on both ends. -/
def shrink (s : Segment2D) (amount : Float) : Segment2D :=
  s.extend (-amount)

end Segment2D

end Linalg
