/-
  Capsule primitive.

  A capsule is a cylinder with hemispherical caps, defined by:
  - a: first endpoint of the central line segment
  - b: second endpoint of the central line segment
  - radius: the radius of the capsule

  Useful for character collision detection and swept sphere tests.
-/

import Linalg.Vec3
import Linalg.Geometry.AABB
import Linalg.Geometry.Sphere

namespace Linalg

/-- Capsule defined by two endpoints and a radius. -/
structure Capsule where
  a : Vec3      -- First endpoint
  b : Vec3      -- Second endpoint
  radius : Float
deriving Repr, BEq, Inhabited

namespace Capsule

/-- Create a capsule from two endpoints and radius. -/
def fromEndpoints (a b : Vec3) (radius : Float) : Capsule :=
  { a := a, b := b, radius := radius }

/-- Create a capsule from center, direction, half-height, and radius. -/
def fromCenterDir (center direction : Vec3) (halfHeight radius : Float) : Capsule :=
  let offset := direction.normalize.scale halfHeight
  { a := center.sub offset, b := center.add offset, radius := radius }

/-- Create a vertical capsule (along Y axis). -/
def vertical (center : Vec3) (halfHeight radius : Float) : Capsule :=
  fromCenterDir center Vec3.unitY halfHeight radius

/-- Get the center of the capsule. -/
def center (c : Capsule) : Vec3 :=
  Vec3.lerp c.a c.b 0.5

/-- Get the direction from a to b (not normalized). -/
def segment (c : Capsule) : Vec3 :=
  c.b.sub c.a

/-- Get the normalized direction from a to b. -/
def direction (c : Capsule) : Vec3 :=
  c.segment.normalize

/-- Get the length of the central line segment. -/
def segmentLength (c : Capsule) : Float :=
  c.a.distance c.b

/-- Get the half-length of the central line segment. -/
def halfHeight (c : Capsule) : Float :=
  c.segmentLength / 2.0

/-- Get the total height (segment length + 2 * radius for the caps). -/
def totalHeight (c : Capsule) : Float :=
  c.segmentLength + 2.0 * c.radius

/-- Get the closest point on the central line segment to a given point. -/
def closestPointOnSegment (c : Capsule) (p : Vec3) : Vec3 :=
  let ab := c.segment
  let ap := p.sub c.a
  let lenSq := ab.lengthSquared
  if lenSq < Float.epsilon then
    c.a  -- Degenerate capsule (a == b)
  else
    let t := Float.clamp (ap.dot ab / lenSq) 0.0 1.0
    c.a.add (ab.scale t)

/-- Get the closest point on the capsule surface to a given point. -/
def closestPoint (c : Capsule) (p : Vec3) : Vec3 :=
  let onSegment := c.closestPointOnSegment p
  let dir := p.sub onSegment
  let dist := dir.length
  if dist < Float.epsilon then
    -- Point is on the segment axis, pick arbitrary direction
    onSegment.add (Vec3.unitX.scale c.radius)
  else
    onSegment.add (dir.scale (c.radius / dist))

/-- Get the squared distance from a point to the capsule surface. -/
def distanceSquared (c : Capsule) (p : Vec3) : Float :=
  let onSegment := c.closestPointOnSegment p
  let dist := p.distance onSegment
  let surfaceDist := dist - c.radius
  if surfaceDist <= 0.0 then 0.0
  else surfaceDist * surfaceDist

/-- Get the distance from a point to the capsule surface. -/
def distance (c : Capsule) (p : Vec3) : Float :=
  let onSegment := c.closestPointOnSegment p
  let dist := p.distance onSegment
  Float.max (dist - c.radius) 0.0

/-- Get the signed distance from a point to the capsule (negative if inside). -/
def signedDistance (c : Capsule) (p : Vec3) : Float :=
  let onSegment := c.closestPointOnSegment p
  p.distance onSegment - c.radius

/-- Check if a point is inside the capsule. -/
def containsPoint (c : Capsule) (p : Vec3) : Bool :=
  c.signedDistance p <= 0.0

/-- Get the AABB that encloses this capsule. -/
def toAABB (c : Capsule) : AABB :=
  let minPt := Vec3.min c.a c.b
  let maxPt := Vec3.max c.a c.b
  let radiusVec := Vec3.mk c.radius c.radius c.radius
  AABB.fromMinMax (minPt.sub radiusVec) (maxPt.add radiusVec)

/-- Get the bounding sphere that encloses this capsule. -/
def boundingSphere (c : Capsule) : Sphere :=
  let center := c.center
  let halfLen := c.segmentLength / 2.0
  { center := center, radius := halfLen + c.radius }

/-- Get the volume of the capsule. -/
def volume (c : Capsule) : Float :=
  let r := c.radius
  let h := c.segmentLength
  -- Volume = cylinder + sphere
  -- Cylinder: π * r² * h
  -- Sphere: (4/3) * π * r³
  Float.pi * r * r * h + (4.0 / 3.0) * Float.pi * r * r * r

/-- Get the surface area of the capsule. -/
def surfaceArea (c : Capsule) : Float :=
  let r := c.radius
  let h := c.segmentLength
  -- Surface = cylinder lateral area + sphere
  -- Cylinder lateral: 2 * π * r * h
  -- Sphere: 4 * π * r²
  2.0 * Float.pi * r * h + 4.0 * Float.pi * r * r

/-- Translate the capsule by a vector. -/
def translate (c : Capsule) (v : Vec3) : Capsule :=
  { c with a := c.a.add v, b := c.b.add v }

/-- Scale the capsule uniformly from origin. -/
def scale (c : Capsule) (s : Float) : Capsule :=
  { a := c.a.scale s, b := c.b.scale s, radius := c.radius * s }

/-- Check if two capsules are approximately equal. -/
def approxEq (c1 c2 : Capsule) (eps : Float := Float.epsilon) : Bool :=
  c1.a.approxEq c2.a eps && c1.b.approxEq c2.b eps &&
  Float.approxEq c1.radius c2.radius eps

/-- Get the closest points between two line segments.
    Returns (point on segment 1, point on segment 2). -/
private def closestPointsBetweenSegments (a1 b1 a2 b2 : Vec3) : Vec3 × Vec3 := Id.run do
  let d1 := b1.sub a1  -- Direction of segment 1
  let d2 := b2.sub a2  -- Direction of segment 2
  let r := a1.sub a2

  let a := d1.dot d1  -- Squared length of segment 1
  let e := d2.dot d2  -- Squared length of segment 2
  let f := d2.dot r

  -- Check for degenerate cases (one or both segments are points)
  if a < Float.epsilon && e < Float.epsilon then
    return (a1, a2)

  if a < Float.epsilon then
    -- First segment is a point
    let t := Float.clamp (f / e) 0.0 1.0
    return (a1, a2.add (d2.scale t))

  let c := d1.dot r
  if e < Float.epsilon then
    -- Second segment is a point
    let s := Float.clamp (-c / a) 0.0 1.0
    return (a1.add (d1.scale s), a2)

  -- General case: both segments have length
  let b := d1.dot d2
  let denom := a * e - b * b

  -- If segments are parallel, pick a point on segment 1
  let mut s := 0.0
  if denom > Float.epsilon then
    s := Float.clamp ((b * f - c * e) / denom) 0.0 1.0

  -- Compute closest point on segment 2
  let mut t := (b * s + f) / e

  -- Clamp t and recompute s if needed
  if t < 0.0 then
    t := 0.0
    s := Float.clamp (-c / a) 0.0 1.0
  else if t > 1.0 then
    t := 1.0
    s := Float.clamp ((b - c) / a) 0.0 1.0

  let p1 := a1.add (d1.scale s)
  let p2 := a2.add (d2.scale t)
  return (p1, p2)

/-- Test if two capsules intersect. -/
def intersectsCapsule (c1 c2 : Capsule) : Bool :=
  let (p1, p2) := closestPointsBetweenSegments c1.a c1.b c2.a c2.b
  let dist := p1.distance p2
  dist <= c1.radius + c2.radius

/-- Test if the capsule intersects a sphere. -/
def intersectsSphere (c : Capsule) (sphere : Sphere) : Bool :=
  let closest := c.closestPointOnSegment sphere.center
  let dist := closest.distance sphere.center
  dist <= c.radius + sphere.radius

/-- Test if the capsule intersects an AABB. -/
def intersectsAABB (c : Capsule) (aabb : AABB) : Bool :=
  -- Conservative test: check if AABB of capsule intersects, then do precise test
  let capsuleAABB := c.toAABB
  if capsuleAABB.min.x > aabb.max.x || capsuleAABB.max.x < aabb.min.x ||
     capsuleAABB.min.y > aabb.max.y || capsuleAABB.max.y < aabb.min.y ||
     capsuleAABB.min.z > aabb.max.z || capsuleAABB.max.z < aabb.min.z then
    false
  else
    -- More precise: check distance from segment to AABB
    let closest := aabb.closestPoint (c.closestPointOnSegment aabb.center)
    let onSegment := c.closestPointOnSegment closest
    onSegment.distance closest <= c.radius

end Capsule

end Linalg
