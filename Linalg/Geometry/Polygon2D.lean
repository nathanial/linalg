/-
  Polygon2D - 2D Polygon Primitive

  A polygon defined by an array of vertices, with support for:
  - Winding order detection and normalization
  - Area calculation (signed and unsigned)
  - Point-in-polygon tests
  - Centroid calculation
  - Bounding box
-/

import Linalg.Core
import Linalg.Vec2
import Linalg.Geometry.Line2D

namespace Linalg

/-- Winding order of a polygon. -/
inductive WindingOrder
  | clockwise
  | counterClockwise
  | degenerate  -- Area is zero (collinear points)
deriving Repr, BEq, Inhabited

/-- A 2D polygon defined by an array of vertices in order. -/
structure Polygon2D where
  vertices : Array Vec2
deriving Repr, Inhabited

namespace Polygon2D

-- ============================================================================
-- Constructors
-- ============================================================================

/-- Create a polygon from vertices. -/
def fromVertices (vertices : Array Vec2) : Polygon2D := { vertices := vertices }

/-- Create a triangle. -/
def triangle (v0 v1 v2 : Vec2) : Polygon2D :=
  { vertices := #[v0, v1, v2] }

/-- Create a rectangle from min and max corners. -/
def rectangle (min max : Vec2) : Polygon2D :=
  { vertices := #[
    min,
    Vec2.mk max.x min.y,
    max,
    Vec2.mk min.x max.y
  ]}

/-- Create a rectangle centered at origin with given width and height. -/
def centeredRectangle (width height : Float) : Polygon2D :=
  let hw := width / 2.0
  let hh := height / 2.0
  rectangle (Vec2.mk (-hw) (-hh)) (Vec2.mk hw hh)

/-- Create a regular polygon with n sides, centered at origin. -/
def regular (n : Nat) (radius : Float) : Polygon2D :=
  if n < 3 then { vertices := #[] }
  else
    let angleStep := 2.0 * Float.pi / n.toFloat
    let verts := Id.run do
      let mut arr : Array Vec2 := #[]
      for i in [:n] do
        let angle := i.toFloat * angleStep - Float.pi / 2.0  -- Start at top
        arr := arr.push (Vec2.mk (radius * Float.cos angle) (radius * Float.sin angle))
      return arr
    { vertices := verts }

-- ============================================================================
-- Basic Properties
-- ============================================================================

/-- Number of vertices. -/
def vertexCount (p : Polygon2D) : Nat := p.vertices.size

/-- Number of edges (same as vertex count for closed polygon). -/
def edgeCount (p : Polygon2D) : Nat := p.vertices.size

/-- Check if the polygon is valid (has at least 3 vertices). -/
def isValid (p : Polygon2D) : Bool := p.vertices.size >= 3

/-- Get vertex at index (with wraparound). -/
def vertex (p : Polygon2D) (i : Nat) : Vec2 :=
  if p.vertices.isEmpty then Vec2.zero
  else p.vertices[i % p.vertices.size]!

/-- Get edge as a segment (from vertex i to vertex i+1). -/
def edge (p : Polygon2D) (i : Nat) : Segment2D :=
  { a := p.vertex i, b := p.vertex (i + 1) }

-- ============================================================================
-- Area and Winding
-- ============================================================================

/-- Signed area (positive for counter-clockwise, negative for clockwise).
    Uses the shoelace formula. -/
def signedArea (p : Polygon2D) : Float :=
  if p.vertices.size < 3 then 0.0
  else Id.run do
    let mut sum := 0.0
    for i in [:p.vertices.size] do
      let v0 := p.vertex i
      let v1 := p.vertex (i + 1)
      sum := sum + v0.cross v1
    return sum / 2.0

/-- Absolute area of the polygon. -/
def area (p : Polygon2D) : Float := Float.abs' p.signedArea

/-- Determine the winding order of the polygon. -/
def windingOrder (p : Polygon2D) : WindingOrder :=
  let sa := p.signedArea
  if Float.abs' sa < Float.epsilon then WindingOrder.degenerate
  else if sa > 0.0 then WindingOrder.counterClockwise
  else WindingOrder.clockwise

/-- Check if polygon is counter-clockwise. -/
def isCounterClockwise (p : Polygon2D) : Bool :=
  p.signedArea > Float.epsilon

/-- Check if polygon is clockwise. -/
def isClockwise (p : Polygon2D) : Bool :=
  p.signedArea < -Float.epsilon

-- ============================================================================
-- Geometric Properties
-- ============================================================================

/-- Centroid (center of mass) of the polygon. -/
def centroid (p : Polygon2D) : Vec2 :=
  if p.vertices.size < 3 then Vec2.zero
  else
    let sa := p.signedArea
    if Float.abs' sa < Float.epsilon then
      -- Degenerate: return average of vertices
      let sum := p.vertices.foldl (· + ·) Vec2.zero
      sum.scale (1.0 / p.vertices.size.toFloat)
    else Id.run do
      let mut cx := 0.0
      let mut cy := 0.0
      for i in [:p.vertices.size] do
        let v0 := p.vertex i
        let v1 := p.vertex (i + 1)
        let cross := v0.cross v1
        cx := cx + (v0.x + v1.x) * cross
        cy := cy + (v0.y + v1.y) * cross
      let factor := 1.0 / (6.0 * sa)
      return Vec2.mk (cx * factor) (cy * factor)

/-- Perimeter of the polygon. -/
def perimeter (p : Polygon2D) : Float := Id.run do
  let mut sum := 0.0
  for i in [:p.vertices.size] do
    sum := sum + (p.edge i).length
  return sum

/-- Axis-aligned bounding box (returns min, max). -/
def boundingBox (p : Polygon2D) : Vec2 × Vec2 :=
  if p.vertices.isEmpty then (Vec2.zero, Vec2.zero)
  else
    let first := p.vertices[0]!
    let (minV, maxV) := p.vertices.foldl (fun (minV, maxV) v =>
      (Vec2.mk (Float.min minV.x v.x) (Float.min minV.y v.y),
       Vec2.mk (Float.max maxV.x v.x) (Float.max maxV.y v.y))
    ) (first, first)
    (minV, maxV)

-- ============================================================================
-- Point Tests
-- ============================================================================

/-- Check if a point is inside the polygon using ray casting algorithm. -/
def containsPoint (p : Polygon2D) (point : Vec2) : Bool :=
  if p.vertices.size < 3 then false
  else Id.run do
    let mut inside := false
    for i in [:p.vertices.size] do
      let v0 := p.vertex i
      let v1 := p.vertex (i + 1)
      -- Ray casting: count intersections with horizontal ray to the right
      if ((v0.y > point.y) != (v1.y > point.y)) then
        let slope := (v1.x - v0.x) / (v1.y - v0.y)
        let xIntersect := v0.x + slope * (point.y - v0.y)
        if point.x < xIntersect then
          inside := !inside
    return inside

/-- Check if a point is on the boundary of the polygon. -/
def isOnBoundary (p : Polygon2D) (point : Vec2) (tolerance : Float := Float.epsilon) : Bool := Id.run do
  for i in [:p.vertices.size] do
    if (p.edge i).containsPoint point tolerance then
      return true
  return false

/-- Check if a point is inside or on the boundary. -/
def containsPointInclusive (p : Polygon2D) (point : Vec2) (tolerance : Float := Float.epsilon) : Bool :=
  p.containsPoint point || p.isOnBoundary point tolerance

-- ============================================================================
-- Transformations
-- ============================================================================

/-- Reverse the winding order of the polygon. -/
def reverse (p : Polygon2D) : Polygon2D :=
  { vertices := p.vertices.reverse }

/-- Ensure the polygon is counter-clockwise. -/
def makeCounterClockwise (p : Polygon2D) : Polygon2D :=
  if p.isClockwise then p.reverse else p

/-- Ensure the polygon is clockwise. -/
def makeClockwise (p : Polygon2D) : Polygon2D :=
  if p.isCounterClockwise then p.reverse else p

/-- Translate the polygon. -/
def translate (p : Polygon2D) (offset : Vec2) : Polygon2D :=
  { vertices := p.vertices.map (· + offset) }

/-- Scale the polygon uniformly from origin. -/
def scale (p : Polygon2D) (factor : Float) : Polygon2D :=
  { vertices := p.vertices.map (·.scale factor) }

/-- Scale the polygon from its centroid. -/
def scaleFromCenter (p : Polygon2D) (factor : Float) : Polygon2D :=
  let c := p.centroid
  { vertices := p.vertices.map (fun v => c + (v - c).scale factor) }

/-- Rotate the polygon around origin. -/
def rotate (p : Polygon2D) (angle : Float) : Polygon2D :=
  let cosA := Float.cos angle
  let sinA := Float.sin angle
  { vertices := p.vertices.map (fun v =>
    Vec2.mk (v.x * cosA - v.y * sinA) (v.x * sinA + v.y * cosA)) }

/-- Rotate the polygon around its centroid. -/
def rotateAroundCenter (p : Polygon2D) (angle : Float) : Polygon2D :=
  let c := p.centroid
  let translated := p.translate (Vec2.zero - c)
  let rotated := translated.rotate angle
  rotated.translate c

-- ============================================================================
-- Convexity
-- ============================================================================

/-- Check if the polygon is convex. -/
def isConvex (p : Polygon2D) : Bool :=
  if p.vertices.size < 3 then false
  else Id.run do
    let mut sign : Option Bool := none
    for i in [:p.vertices.size] do
      let v0 := p.vertex i
      let v1 := p.vertex (i + 1)
      let v2 := p.vertex (i + 2)
      let cross := (v1 - v0).cross (v2 - v1)
      if Float.abs' cross > Float.epsilon then
        let positive := cross > 0.0
        match sign with
        | none => sign := some positive
        | some s => if s != positive then return false
    return true

-- ============================================================================
-- Edge Queries
-- ============================================================================

/-- Get all edges as segments. -/
def edges (p : Polygon2D) : Array Segment2D := Id.run do
  let mut result : Array Segment2D := #[]
  for i in [:p.vertices.size] do
    result := result.push (p.edge i)
  return result

/-- Find the closest point on the polygon boundary to a given point. -/
def closestPointOnBoundary (p : Polygon2D) (point : Vec2) : Vec2 :=
  if p.vertices.isEmpty then Vec2.zero
  else Id.run do
    let mut closest := p.vertices[0]!
    let mut minDistSq := Float.infinity
    for i in [:p.vertices.size] do
      let edgeClosest := (p.edge i).closestPoint point
      let distSq := (point - edgeClosest).lengthSquared
      if distSq < minDistSq then
        minDistSq := distSq
        closest := edgeClosest
    return closest

/-- Distance from a point to the polygon boundary. -/
def distanceToBoundary (p : Polygon2D) (point : Vec2) : Float :=
  (point - p.closestPointOnBoundary point).length

end Polygon2D

end Linalg
