/-
  Circle - 2D Circle Primitive

  A circle defined by center and radius, with containment tests,
  intersection queries, and utility functions.
-/

import Linalg.Core
import Linalg.Vec2

namespace Linalg

/-- A 2D circle defined by center and radius. -/
structure Circle where
  center : Vec2
  radius : Float
deriving Repr, BEq, Inhabited

namespace Circle

-- ============================================================================
-- Constructors
-- ============================================================================

/-- Create a circle at the origin with given radius. -/
def atOrigin (radius : Float) : Circle :=
  { center := Vec2.zero, radius := radius }

/-- Create a unit circle at the origin. -/
def unit : Circle := atOrigin 1.0

/-- Create a circle from center coordinates and radius. -/
def mk' (cx cy radius : Float) : Circle :=
  { center := Vec2.mk cx cy, radius := radius }

/-- Create circle from three points on its circumference. -/
def fromThreePoints (p1 p2 p3 : Vec2) : Option Circle :=
  -- Using circumcenter formula
  let ax := p1.x; let ay := p1.y
  let bx := p2.x; let by_ := p2.y
  let cx := p3.x; let cy := p3.y
  let d := 2.0 * (ax * (by_ - cy) + bx * (cy - ay) + cx * (ay - by_))
  if Float.abs' d < Float.epsilon then
    none  -- Points are collinear
  else
    let aSq := ax * ax + ay * ay
    let bSq := bx * bx + by_ * by_
    let cSq := cx * cx + cy * cy
    let ux := (aSq * (by_ - cy) + bSq * (cy - ay) + cSq * (ay - by_)) / d
    let uy := (aSq * (cx - bx) + bSq * (ax - cx) + cSq * (bx - ax)) / d
    let center := Vec2.mk ux uy
    let radius := (p1 - center).length
    some { center := center, radius := radius }

-- ============================================================================
-- Properties
-- ============================================================================

/-- Diameter of the circle. -/
def diameter (c : Circle) : Float := c.radius * 2.0

/-- Circumference of the circle. -/
def circumference (c : Circle) : Float := 2.0 * Float.pi * c.radius

/-- Area of the circle. -/
def area (c : Circle) : Float := Float.pi * c.radius * c.radius

-- ============================================================================
-- Containment Tests
-- ============================================================================

/-- Check if a point is inside the circle (exclusive of boundary). -/
def containsPoint (c : Circle) (p : Vec2) : Bool :=
  (p - c.center).lengthSquared < c.radius * c.radius

/-- Check if a point is on or inside the circle. -/
def containsPointInclusive (c : Circle) (p : Vec2) : Bool :=
  (p - c.center).lengthSquared <= c.radius * c.radius

/-- Check if a point is on the boundary of the circle. -/
def isOnBoundary (c : Circle) (p : Vec2) (tolerance : Float := Float.epsilon) : Bool :=
  let dist := (p - c.center).length
  Float.abs' (dist - c.radius) < tolerance

/-- Signed distance from point to circle boundary (negative if inside). -/
def signedDistance (c : Circle) (p : Vec2) : Float :=
  (p - c.center).length - c.radius

/-- Distance from point to nearest point on circle. -/
def distance (c : Circle) (p : Vec2) : Float :=
  Float.abs' (c.signedDistance p)

-- ============================================================================
-- Geometric Queries
-- ============================================================================

/-- Closest point on circle boundary to a given point. -/
def closestPoint (c : Circle) (p : Vec2) : Vec2 :=
  let dir := p - c.center
  let len := dir.length
  if len < Float.epsilon then
    -- Point is at center, return arbitrary point on boundary
    c.center + Vec2.mk c.radius 0.0
  else
    c.center + dir.scale (c.radius / len)

/-- Get a point on the circle at a given angle (radians from positive x-axis). -/
def pointAtAngle (c : Circle) (angle : Float) : Vec2 :=
  c.center + Vec2.mk (c.radius * Float.cos angle) (c.radius * Float.sin angle)

/-- Get the tangent vector at a point on the circle (counter-clockwise). -/
def tangentAt (c : Circle) (p : Vec2) : Vec2 :=
  let radial := (p - c.center).normalize
  radial.perpendicular

-- ============================================================================
-- Circle-Circle Relationships
-- ============================================================================

/-- Check if two circles intersect. -/
def intersects (c1 c2 : Circle) : Bool :=
  let distSq := (c2.center - c1.center).lengthSquared
  let sumRadii := c1.radius + c2.radius
  distSq <= sumRadii * sumRadii

/-- Check if two circles are separate (not touching or overlapping). -/
def isSeparate (c1 c2 : Circle) : Bool :=
  let dist := (c2.center - c1.center).length
  dist > c1.radius + c2.radius

/-- Check if this circle contains another circle entirely. -/
def containsCircle (outer inner : Circle) : Bool :=
  let dist := (inner.center - outer.center).length
  dist + inner.radius <= outer.radius

/-- Check if two circles are tangent (touching at exactly one point). -/
def isTangent (c1 c2 : Circle) (tolerance : Float := 0.0001) : Bool :=
  let dist := (c2.center - c1.center).length
  let sumRadii := c1.radius + c2.radius
  let diffRadii := Float.abs' (c1.radius - c2.radius)
  Float.abs' (dist - sumRadii) < tolerance ||  -- External tangent
  Float.abs' (dist - diffRadii) < tolerance    -- Internal tangent

/-- Get intersection points of two circles (0, 1, or 2 points). -/
def intersectionPoints (c1 c2 : Circle) : Array Vec2 :=
  let d := (c2.center - c1.center).length
  let r1 := c1.radius
  let r2 := c2.radius

  -- No intersection cases
  if d > r1 + r2 then #[]  -- Too far apart
  else if d < Float.abs' (r1 - r2) then #[]  -- One inside other
  else if d < Float.epsilon then #[]  -- Concentric
  else
    -- Calculate intersection points
    let a := (r1 * r1 - r2 * r2 + d * d) / (2.0 * d)
    let hSq := r1 * r1 - a * a
    if hSq < 0.0 then #[]
    else
      let h := Float.sqrt hSq
      let dir := (c2.center - c1.center).normalize
      let midpoint := c1.center + dir.scale a
      let perp := dir.perpendicular
      if h < Float.epsilon then
        -- Single point (tangent)
        #[midpoint]
      else
        -- Two points
        #[midpoint + perp.scale h, midpoint - perp.scale h]

-- ============================================================================
-- Bounding
-- ============================================================================

/-- Get the axis-aligned bounding box of the circle. -/
def boundingBox (c : Circle) : Vec2 × Vec2 :=
  let min := Vec2.mk (c.center.x - c.radius) (c.center.y - c.radius)
  let max := Vec2.mk (c.center.x + c.radius) (c.center.y + c.radius)
  (min, max)

-- ============================================================================
-- Transformations
-- ============================================================================

/-- Translate the circle. -/
def translate (c : Circle) (offset : Vec2) : Circle :=
  { c with center := c.center + offset }

/-- Scale the circle uniformly from its center. -/
def scale (c : Circle) (factor : Float) : Circle :=
  { c with radius := c.radius * Float.abs' factor }

/-- Scale the circle uniformly from a given origin. -/
def scaleFrom (c : Circle) (origin : Vec2) (factor : Float) : Circle :=
  let newCenter := origin + (c.center - origin).scale factor
  { center := newCenter, radius := c.radius * Float.abs' factor }

end Circle

end Linalg
