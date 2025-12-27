/-
  2D Collision Detection

  Implements:
  - SAT (Separating Axis Theorem) for convex polygon collision
  - GJK (Gilbert-Johnson-Keerthi) for general convex shape collision

  Both algorithms work with convex shapes only.
-/

import Linalg.Core
import Linalg.Vec2
import Linalg.Geometry.Polygon2D
import Linalg.Geometry.Circle
import Linalg.Geometry.AABB2D

namespace Linalg

-- ============================================================================
-- Collision Result Types
-- ============================================================================

/-- Result of a collision test with penetration information. -/
structure CollisionResult where
  /-- Whether the shapes are colliding. -/
  colliding : Bool
  /-- Minimum Translation Vector - direction to push shapes apart.
      Points from shape A toward shape B. -/
  mtv : Vec2
  /-- Penetration depth (magnitude of MTV). -/
  depth : Float
deriving Repr, BEq, Inhabited

namespace CollisionResult

/-- No collision result. -/
def none : CollisionResult := ⟨false, Vec2.zero, 0.0⟩

/-- Create a collision result with MTV. -/
def collision (mtv : Vec2) (depth : Float) : CollisionResult :=
  ⟨true, mtv, depth⟩

end CollisionResult

-- ============================================================================
-- Support Function Typeclass (for GJK)
-- ============================================================================

/-- Typeclass for shapes that can provide a support point.
    The support function returns the point on the shape's boundary
    that is furthest in the given direction. -/
class Support2D (α : Type) where
  /-- Get the support point in the given direction. -/
  support : α → Vec2 → Vec2

-- ============================================================================
-- SAT (Separating Axis Theorem)
-- ============================================================================

namespace SAT

/-- Project a polygon onto an axis and return the min/max interval. -/
def projectPolygon (poly : Polygon2D) (axis : Vec2) : Float × Float :=
  if poly.vertices.isEmpty then (0.0, 0.0)
  else
    let first := axis.dot poly.vertices[0]!
    poly.vertices.foldl (fun (minP, maxP) v =>
      let proj := axis.dot v
      (Float.min minP proj, Float.max maxP proj)
    ) (first, first)

/-- Project a circle onto an axis and return the min/max interval. -/
def projectCircle (circle : Circle) (axis : Vec2) : Float × Float :=
  let centerProj := axis.dot circle.center
  (centerProj - circle.radius, centerProj + circle.radius)

/-- Check if two intervals overlap and return the overlap amount. -/
def intervalOverlap (min1 max1 min2 max2 : Float) : Option Float :=
  let overlap := Float.min max1 max2 - Float.max min1 min2
  if overlap > 0.0 then some overlap else none

/-- Get edge normal (perpendicular to edge, pointing outward for CCW polygon). -/
def edgeNormal (a b : Vec2) : Vec2 :=
  let edge := b - a
  Vec2.mk edge.y (-edge.x) |>.normalize

/-- Test collision between two convex polygons using SAT.
    Returns collision result with MTV pointing from poly1 toward poly2. -/
def polygonPolygon (poly1 poly2 : Polygon2D) : CollisionResult :=
  if poly1.vertices.size < 3 || poly2.vertices.size < 3 then
    CollisionResult.none
  else
    -- Ensure both polygons are CCW for consistent normal direction
    let p1 := poly1.makeCounterClockwise
    let p2 := poly2.makeCounterClockwise

    -- Collect all potential separating axes (edge normals from both polygons)
    let axes1 := collectAxes p1
    let axes2 := collectAxes p2

    -- Find minimum overlap across all axes
    match findMinOverlap p1 p2 (axes1 ++ axes2) with
    | none => CollisionResult.none
    | some (axis, depth) =>
      -- Ensure MTV points from poly1 center toward poly2 center
      let center1 := p1.centroid
      let center2 := p2.centroid
      let dir := center2 - center1
      let mtv := if axis.dot dir >= 0.0 then axis else axis.scale (-1.0)
      CollisionResult.collision mtv depth
where
  collectAxes (poly : Polygon2D) : Array Vec2 := Id.run do
    let mut axes : Array Vec2 := #[]
    for i in [:poly.vertices.size] do
      let v1 := poly.vertex i
      let v2 := poly.vertex (i + 1)
      axes := axes.push (edgeNormal v1 v2)
    return axes

  findMinOverlap (p1 p2 : Polygon2D) (axes : Array Vec2) : Option (Vec2 × Float) :=
    axes.foldl (fun best axis =>
      let (min1, max1) := projectPolygon p1 axis
      let (min2, max2) := projectPolygon p2 axis
      match intervalOverlap min1 max1 min2 max2 with
      | none => none  -- Found separating axis
      | some overlap =>
        match best with
        | none => none  -- Already found a separating axis
        | some (_, minOverlap) =>
          if overlap < minOverlap then some (axis, overlap)
          else best
    ) (some (Vec2.zero, Float.infinity))

/-- Test collision between a convex polygon and a circle using SAT. -/
def polygonCircle (poly : Polygon2D) (circle : Circle) : CollisionResult :=
  if poly.vertices.size < 3 then CollisionResult.none
  else
    let p := poly.makeCounterClockwise

    -- Axes: edge normals + axis from closest vertex to circle center
    let edgeAxes := collectEdgeAxes p
    let closestVertex := findClosestVertex p circle.center
    let vertexAxis := (circle.center - closestVertex).normalize
    let axes := edgeAxes.push vertexAxis

    match findMinOverlapCircle p circle axes with
    | none => CollisionResult.none
    | some (axis, depth) =>
      -- Ensure MTV points from polygon toward circle
      let polyCenter := p.centroid
      let dir := circle.center - polyCenter
      let mtv := if axis.dot dir >= 0.0 then axis else axis.scale (-1.0)
      CollisionResult.collision mtv depth
where
  collectEdgeAxes (poly : Polygon2D) : Array Vec2 := Id.run do
    let mut axes : Array Vec2 := #[]
    for i in [:poly.vertices.size] do
      let v1 := poly.vertex i
      let v2 := poly.vertex (i + 1)
      axes := axes.push (edgeNormal v1 v2)
    return axes

  findClosestVertex (poly : Polygon2D) (point : Vec2) : Vec2 :=
    poly.vertices.foldl (fun (closest, minDist) v =>
      let dist := (v - point).lengthSquared
      if dist < minDist then (v, dist) else (closest, minDist)
    ) (poly.vertices[0]!, Float.infinity) |>.1

  findMinOverlapCircle (poly : Polygon2D) (circle : Circle) (axes : Array Vec2) : Option (Vec2 × Float) :=
    axes.foldl (fun best axis =>
      let (min1, max1) := projectPolygon poly axis
      let (min2, max2) := projectCircle circle axis
      match intervalOverlap min1 max1 min2 max2 with
      | none => none
      | some overlap =>
        match best with
        | none => none
        | some (_, minOverlap) =>
          if overlap < minOverlap then some (axis, overlap)
          else best
    ) (some (Vec2.zero, Float.infinity))

/-- Test collision between two circles. -/
def circleCircle (c1 c2 : Circle) : CollisionResult :=
  let diff := c2.center - c1.center
  let distSq := diff.lengthSquared
  let radiusSum := c1.radius + c2.radius

  if distSq >= radiusSum * radiusSum then
    CollisionResult.none
  else
    let dist := Float.sqrt distSq
    if dist < Float.epsilon then
      -- Circles are concentric, pick arbitrary direction
      CollisionResult.collision Vec2.unitX radiusSum
    else
      let normal := diff.scale (1.0 / dist)
      let depth := radiusSum - dist
      CollisionResult.collision normal depth

/-- Quick boolean test: do two convex polygons intersect? -/
def polygonsIntersect (poly1 poly2 : Polygon2D) : Bool :=
  (polygonPolygon poly1 poly2).colliding

/-- Quick boolean test: does a polygon intersect a circle? -/
def polygonIntersectsCircle (poly : Polygon2D) (circle : Circle) : Bool :=
  (polygonCircle poly circle).colliding

/-- Quick boolean test: do two circles intersect? -/
def circlesIntersect (c1 c2 : Circle) : Bool :=
  (circleCircle c1 c2).colliding

end SAT

-- ============================================================================
-- GJK (Gilbert-Johnson-Keerthi Algorithm)
-- ============================================================================

namespace GJK

/-- Support point in Minkowski difference with original points tracked. -/
structure MinkowskiPoint where
  /-- Point in Minkowski difference space (a - b). -/
  point : Vec2
  /-- Support point from shape A. -/
  supportA : Vec2
  /-- Support point from shape B. -/
  supportB : Vec2
deriving Repr, Inhabited

/-- Get support point in Minkowski difference A - B. -/
def minkowskiSupport {α β : Type} [Support2D α] [Support2D β] (a : α) (b : β) (dir : Vec2) : MinkowskiPoint :=
  let sa := Support2D.support a dir
  let sb := Support2D.support b (dir.scale (-1.0))
  ⟨sa - sb, sa, sb⟩

/-- Triple product: (A × B) × C, useful for finding perpendicular in 2D. -/
def tripleProduct (a b c : Vec2) : Vec2 :=
  -- In 2D: (A × B) × C = B * (A · C) - A * (B · C)
  let ac := a.dot c
  let bc := b.dot c
  Vec2.mk (b.x * ac - a.x * bc) (b.y * ac - a.y * bc)

/-- Process a 2-point simplex (line). -/
private def processLine (simplex : Array MinkowskiPoint) : Bool × Array MinkowskiPoint × Vec2 :=
  let a := simplex[1]!.point  -- Most recent point
  let b := simplex[0]!.point
  let ab := b - a
  let ao := Vec2.zero - a

  -- Get perpendicular toward origin
  let abPerp := tripleProduct ab ao ab
  if abPerp.lengthSquared < Float.epsilon then
    -- Origin is on the line, use perpendicular
    (false, simplex, Vec2.mk (-ab.y) ab.x)
  else
    (false, simplex, abPerp)

/-- Process a 3-point simplex (triangle). -/
private def processTriangle (simplex : Array MinkowskiPoint) : Bool × Array MinkowskiPoint × Vec2 :=
  let a := simplex[2]!.point  -- Most recent point
  let b := simplex[1]!.point
  let c := simplex[0]!.point
  let ao := Vec2.zero - a
  let ab := b - a
  let ac := c - a

  -- Check which edge region contains the origin
  let abPerp := tripleProduct ac ab ab
  let acPerp := tripleProduct ab ac ac

  if abPerp.dot ao > 0.0 then
    -- Origin is outside AB edge
    (false, #[simplex[1]!, simplex[2]!], abPerp)
  else if acPerp.dot ao > 0.0 then
    -- Origin is outside AC edge
    (false, #[simplex[0]!, simplex[2]!], acPerp)
  else
    -- Origin is inside the triangle
    (true, simplex, Vec2.zero)

/-- Check if simplex contains the origin and update search direction.
    Returns (containsOrigin, newSimplex, newDirection). -/
def processSimplex (simplex : Array MinkowskiPoint) : Bool × Array MinkowskiPoint × Vec2 :=
  match simplex.size with
  | 2 => processLine simplex
  | 3 => processTriangle simplex
  | _ => (false, simplex, Vec2.zero)

/-- Maximum GJK iterations to prevent infinite loops. -/
def maxIterations : Nat := 32

/-- GJK loop for intersection test. -/
private partial def gjkLoop {α β : Type} [Support2D α] [Support2D β]
    (a : α) (b : β) (simplex : Array MinkowskiPoint) (dir : Vec2) (iter : Nat) : Bool :=
  if iter >= maxIterations then false
  else if dir.lengthSquared < Float.epsilon then true  -- Origin is at support point
  else
    -- Get new support point in direction toward origin
    let support := minkowskiSupport a b dir
    -- Check if we passed the origin
    if support.point.dot dir <= 0.0 then
      false  -- Didn't pass origin, no intersection
    else
      -- Add to simplex and check if it contains origin
      let newSimplex := simplex.push support
      let (containsOrigin, simplex', newDir) := processSimplex newSimplex
      if containsOrigin then true
      else gjkLoop a b simplex' newDir (iter + 1)

/-- Run GJK algorithm to test if two convex shapes intersect.
    Returns true if shapes intersect. -/
def intersects {α β : Type} [Support2D α] [Support2D β] (a : α) (b : β) : Bool :=
  -- Initial direction (arbitrary, toward B from A would be ideal but we use X axis)
  let initialDir := Vec2.unitX

  -- Get first support point
  let support1 := minkowskiSupport a b initialDir
  if support1.point.dot initialDir <= 0.0 then
    false  -- No intersection possible in this direction
  else
    -- Start simplex with first point, search toward origin
    let dir := Vec2.zero - support1.point
    gjkLoop a b #[support1] dir 0

/-- GJK loop that returns the final simplex. -/
private partial def gjkLoopWithSimplex {α β : Type} [Support2D α] [Support2D β]
    (a : α) (b : β) (simplex : Array MinkowskiPoint) (dir : Vec2) (iter : Nat) : Option (Array MinkowskiPoint) :=
  if iter >= maxIterations then none
  else if dir.lengthSquared < Float.epsilon then some simplex
  else
    let support := minkowskiSupport a b dir
    if support.point.dot dir <= 0.0 then none
    else
      let newSimplex := simplex.push support
      let (containsOrigin, simplex', newDir) := processSimplex newSimplex
      if containsOrigin then some simplex'
      else gjkLoopWithSimplex a b simplex' newDir (iter + 1)

/-- Find closest edge to origin in polytope. -/
private def findClosestEdge (polytope : Array MinkowskiPoint) : Vec2 × Float × Nat :=
  let n := polytope.size
  (List.range n).foldl (fun (bestNormal, bestDist, bestIdx) i =>
    let j := (i + 1) % n
    let a := polytope[i]!.point
    let b := polytope[j]!.point
    let edge := b - a
    -- Normal pointing outward (away from origin)
    let normal := Vec2.mk edge.y (-edge.x) |>.normalize
    let dist := normal.dot a
    -- Ensure normal points away from origin
    let (normal, dist) := if dist < 0.0 then (normal.scale (-1.0), -dist) else (normal, dist)
    if dist < bestDist then (normal, dist, i)
    else (bestNormal, bestDist, bestIdx)
  ) (Vec2.zero, Float.infinity, 0)

/-- Simplified EPA to find penetration depth from final simplex. -/
private def epaLite (simplex : Array MinkowskiPoint) : CollisionResult :=
  if simplex.size < 3 then
    -- Degenerate case - shapes are just touching
    CollisionResult.collision Vec2.unitX Float.epsilon
  else
    -- Find closest edge to origin
    let (edgeNormal, edgeDist, _) := findClosestEdge simplex
    -- Use edge normal as MTV, distance as depth
    CollisionResult.collision edgeNormal edgeDist

/-- GJK with collision result (includes approximate MTV via EPA-lite). -/
def collision {α β : Type} [Support2D α] [Support2D β] (a : α) (b : β) : CollisionResult :=
  let initialDir := Vec2.unitX
  let support1 := minkowskiSupport a b initialDir

  if support1.point.dot initialDir <= 0.0 then
    CollisionResult.none
  else
    let dir := Vec2.zero - support1.point
    match gjkLoopWithSimplex a b #[support1] dir 0 with
    | none => CollisionResult.none
    | some simplex => epaLite simplex

end GJK

-- ============================================================================
-- Support2D Instances
-- ============================================================================

instance : Support2D Polygon2D where
  support poly dir :=
    if poly.vertices.isEmpty then Vec2.zero
    else
      poly.vertices.foldl (fun (best, bestDot) v =>
        let d := dir.dot v
        if d > bestDot then (v, d) else (best, bestDot)
      ) (poly.vertices[0]!, dir.dot poly.vertices[0]!) |>.1

instance : Support2D Circle where
  support circle dir :=
    let normalizedDir := if dir.lengthSquared > Float.epsilon
                         then dir.normalize
                         else Vec2.unitX
    circle.center + normalizedDir.scale circle.radius

instance : Support2D AABB2D where
  support aabb dir :=
    Vec2.mk
      (if dir.x >= 0.0 then aabb.max.x else aabb.min.x)
      (if dir.y >= 0.0 then aabb.max.y else aabb.min.y)

-- ============================================================================
-- Convenience Functions
-- ============================================================================

/-- Test if two convex polygons collide using SAT. -/
def collidePolygons (p1 p2 : Polygon2D) : CollisionResult :=
  SAT.polygonPolygon p1 p2

/-- Test if a polygon and circle collide using SAT. -/
def collidePolygonCircle (poly : Polygon2D) (circle : Circle) : CollisionResult :=
  SAT.polygonCircle poly circle

/-- Test if two circles collide. -/
def collideCircles (c1 c2 : Circle) : CollisionResult :=
  SAT.circleCircle c1 c2

/-- Test if two convex shapes collide using GJK. -/
def collideGJK {α β : Type} [Support2D α] [Support2D β] (a : α) (b : β) : CollisionResult :=
  GJK.collision a b

/-- Quick intersection test using GJK. -/
def intersectsGJK {α β : Type} [Support2D α] [Support2D β] (a : α) (b : β) : Bool :=
  GJK.intersects a b

end Linalg
