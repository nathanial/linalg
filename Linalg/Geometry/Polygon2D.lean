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

-- ============================================================================
-- Boolean Operations (Greiner-Hormann, simple polygons)
-- ============================================================================

inductive BooleanOp
  | intersection
  | union
  | difference
deriving Repr, BEq, Inhabited

private structure PolyNode where
  point : Vec2
  next : Nat
  prev : Nat
  intersection : Bool
  alpha : Float
  neighbor : Option Nat
  entry : Bool
  visited : Bool
deriving Repr, Inhabited

private def buildNodes (verts : Array Vec2) : Array PolyNode := Id.run do
  if verts.isEmpty then
    return #[]
  let n := verts.size
  let mut nodes : Array PolyNode := Array.mkEmpty n
  for i in [:n] do
    let next := (i + 1) % n
    let prev := (i + n - 1) % n
    nodes := nodes.push {
      point := verts[i]!
      next := next
      prev := prev
      intersection := false
      alpha := 0.0
      neighbor := none
      entry := false
      visited := false
    }
  return nodes

private def segmentIntersectionParam (p1 p2 q1 q2 : Vec2) : Option (Vec2 × Float × Float) :=
  let r := p2 - p1
  let s := q2 - q1
  let denom := r.cross s
  if Float.abs' denom < Float.epsilon then none
  else
    let qp := q1 - p1
    let t := qp.cross s / denom
    let u := qp.cross r / denom
    let eps := Float.epsilon * 10.0
    if t >= -eps && t <= 1.0 + eps && u >= -eps && u <= 1.0 + eps then
      some (p1 + r.scale t, t, u)
    else none

private def hasIntersectionAt (nodes : Array PolyNode) (start endIdx : Nat) (pt : Vec2) : Bool := Id.run do
  if nodes.isEmpty then
    return false
  let eps := Float.epsilon * 10.0
  let epsSq := eps * eps
  let mut idx := start
  let mut first := true
  let mut found := false
  while (first || idx != endIdx) && !found do
    let node := nodes[idx]!
    if node.intersection then
      if (node.point - pt).lengthSquared <= epsSq then
        found := true
    idx := node.next
    first := false
  return found

private def insertIntersection (nodes : Array PolyNode) (start endIdx newIdx : Nat) (alpha : Float) :
    Array PolyNode := Id.run do
  let mut result := nodes
  let mut curr := start
  let mut next := result[curr]!.next
  while next != endIdx && result[next]!.intersection && result[next]!.alpha < alpha do
    curr := next
    next := result[curr]!.next
  let currNode := result[curr]!
  let nextNode := result[next]!
  result := result.set! newIdx { result[newIdx]! with prev := curr, next := next }
  result := result.set! curr { currNode with next := newIdx }
  result := result.set! next { nextNode with prev := newIdx }
  return result

private def markEntryFlags (nodes : Array PolyNode) (other : Polygon2D) : Array PolyNode := Id.run do
  if nodes.isEmpty then
    return nodes
  let mut result := nodes
  let start := 0
  let mut idx := start
  let mut inside := other.containsPointInclusive result[idx]!.point
  let mut first := true
  while first || idx != start do
    let node := result[idx]!
    if node.intersection then
      let entry := !inside
      result := result.set! idx { node with entry := entry }
      inside := !inside
    idx := result[idx]!.next
    first := false
  return result

private def markVisited (nodes : Array PolyNode) (idx : Nat) : Array PolyNode :=
  let node := nodes[idx]!
  if node.visited then nodes else nodes.set! idx { node with visited := true }

private def findStart (nodes : Array PolyNode) (wantEntry : Bool) : Option Nat := Id.run do
  for i in [:nodes.size] do
    let node := nodes[i]!
    if node.intersection && !node.visited && node.entry == wantEntry then
      return some i
  return none

private def cleanPoints (points : Array Vec2) : Array Vec2 := Id.run do
  let eps := Float.epsilon * 10.0
  let epsSq := eps * eps
  let mut result : Array Vec2 := #[]
  for p in points do
    if result.isEmpty then
      result := result.push p
    else
      let last := result[result.size - 1]!
      if (p - last).lengthSquared > epsSq then
        result := result.push p
  if result.size > 1 then
    let first := result[0]!
    let last := result[result.size - 1]!
    if (first - last).lengthSquared <= epsSq then
      result := result.pop
  return result

private partial def tracePolygon (nodesA nodesB : Array PolyNode) (startIdx : Nat) (op : BooleanOp) :
    Array Vec2 × Array PolyNode × Array PolyNode := Id.run do
  let mut result : Array Vec2 := #[]
  let mut aNodes := nodesA
  let mut bNodes := nodesB
  let mut idx := startIdx
  let mut inA := true
  let startInA := true
  let mut first := true
  while first || !(idx == startIdx && inA == startInA) do
    let node := if inA then aNodes[idx]! else bNodes[idx]!
    if !node.intersection || !node.visited then
      result := result.push node.point
    if node.intersection then
      if inA then
        aNodes := markVisited aNodes idx
        let otherIdx := node.neighbor.getD idx
        bNodes := markVisited bNodes otherIdx
        inA := false
        idx := otherIdx
      else
        bNodes := markVisited bNodes idx
        let otherIdx := node.neighbor.getD idx
        aNodes := markVisited aNodes otherIdx
        inA := true
        idx := otherIdx
      if inA then
        idx := aNodes[idx]!.next
      else
        idx := if op == BooleanOp.difference then bNodes[idx]!.prev else bNodes[idx]!.next
    else
      if inA then
        idx := aNodes[idx]!.next
      else
        idx := if op == BooleanOp.difference then bNodes[idx]!.prev else bNodes[idx]!.next
    first := false
  return (result, aNodes, bNodes)

private def booleanOp (a b : Polygon2D) (op : BooleanOp) : Array Polygon2D := Id.run do
  if !a.isValid || !b.isValid then
    match op with
    | BooleanOp.union => return if a.isValid then #[a] else if b.isValid then #[b] else #[]
    | _ => return #[]

  let mut nodesA := buildNodes a.vertices
  let mut nodesB := buildNodes b.vertices
  let nA := a.vertices.size
  let nB := b.vertices.size
  let mut hasIntersections := false

  for i in [:nA] do
    let a0 := nodesA[i]!.point
    let a1 := nodesA[(i + 1) % nA]!.point
    for j in [:nB] do
      let b0 := nodesB[j]!.point
      let b1 := nodesB[(j + 1) % nB]!.point
      match segmentIntersectionParam a0 a1 b0 b1 with
      | none => ()
      | some (pt, ta, tb) =>
        hasIntersections := true
        let endA := (i + 1) % nA
        let endB := (j + 1) % nB
        if !hasIntersectionAt nodesA i endA pt && !hasIntersectionAt nodesB j endB pt then
          let idxA := nodesA.size
          let idxB := nodesB.size
          nodesA := nodesA.push {
            point := pt
            next := idxA
            prev := idxA
            intersection := true
            alpha := ta
            neighbor := some idxB
            entry := false
            visited := false
          }
          nodesB := nodesB.push {
            point := pt
            next := idxB
            prev := idxB
            intersection := true
            alpha := tb
            neighbor := some idxA
            entry := false
            visited := false
          }
          nodesA := insertIntersection nodesA i endA idxA ta
          nodesB := insertIntersection nodesB j endB idxB tb

  if !hasIntersections then
    let aInB := b.containsPointInclusive a.vertices[0]!
    let bInA := a.containsPointInclusive b.vertices[0]!
    match op with
    | BooleanOp.intersection =>
      if aInB then return #[a]
      else if bInA then return #[b]
      else return #[]
    | BooleanOp.union =>
      if aInB then return #[b]
      else if bInA then return #[a]
      else return #[a, b]
    | BooleanOp.difference =>
      if aInB then return #[]
      else if bInA then return #[a, b.reverse]
      else return #[a]

  nodesA := markEntryFlags nodesA b
  nodesB := markEntryFlags nodesB a

  let wantEntry := match op with
    | BooleanOp.intersection => false
    | BooleanOp.union => true
    | BooleanOp.difference => true

  let mut result : Array Polygon2D := #[]
  let mut startOpt := findStart nodesA wantEntry
  while startOpt.isSome do
    let start := startOpt.get!
    let (points, newA, newB) := tracePolygon nodesA nodesB start op
    nodesA := newA
    nodesB := newB
    let cleaned := cleanPoints points
    if cleaned.size >= 3 then
      result := result.push { vertices := cleaned }
    startOpt := findStart nodesA wantEntry
  return result

/-- Polygon intersection (simple polygons). -/
def intersection (a b : Polygon2D) : Array Polygon2D :=
  booleanOp a b BooleanOp.intersection

/-- Polygon union (simple polygons). -/
def union (a b : Polygon2D) : Array Polygon2D :=
  booleanOp a b BooleanOp.union

/-- Polygon difference (A \ B). Holes may be represented as clockwise polygons. -/
def difference (a b : Polygon2D) : Array Polygon2D :=
  booleanOp a b BooleanOp.difference

-- ============================================================================
-- Offset / Buffer
-- ============================================================================

/-- Offset polygon by distance (outward for CCW, inward for CW). -/
def offset (p : Polygon2D) (distance : Float) : Polygon2D :=
  if p.vertices.size < 3 then p
  else
    let poly := p
    let outwardSign := if poly.isCounterClockwise then -1.0 else 1.0
    let d := distance * outwardSign
    let verts := Id.run do
      let mut out : Array Vec2 := #[]
      for i in [:poly.vertices.size] do
        let prev := poly.vertex (i + poly.vertices.size - 1)
        let curr := poly.vertex i
        let next := poly.vertex (i + 1)
        let dirPrev := (curr - prev).normalize
        let dirNext := (next - curr).normalize
        let nPrev := dirPrev.perpendicular.scale d
        let nNext := dirNext.perpendicular.scale d
        let line1 := Line2D.fromPoints (prev + nPrev) (curr + nPrev)
        let line2 := Line2D.fromPoints (curr + nNext) (next + nNext)
        match Line2D.intersection line1 line2 with
        | some p => out := out.push p
        | none =>
          let fallback := curr + (nPrev + nNext)
          out := out.push fallback
      return out
    { vertices := verts }

-- ============================================================================
-- Convex Hull
-- ============================================================================

/-- Compute the convex hull of a set of points using Andrew's monotone chain algorithm.
    Returns a new polygon with vertices in counter-clockwise order. -/
def convexHull (points : Array Vec2) : Polygon2D :=
  if points.size < 3 then { vertices := points }
  else
    -- Sort points by x, then by y
    let sorted := points.qsort (fun a b =>
      if a.x != b.x then a.x < b.x else a.y < b.y)

    -- Build lower hull
    let lower := sorted.foldl (fun hull p =>
      let hull := removeNonLeftTurns hull p
      hull.push p
    ) #[]

    -- Build upper hull
    let upper := sorted.foldr (fun p hull =>
      let hull := removeNonLeftTurns hull p
      hull.push p
    ) #[]

    -- Concatenate hulls (remove duplicate endpoints)
    let lowerWithoutLast := if lower.size > 0 then lower.pop else lower
    let upperWithoutLast := if upper.size > 0 then upper.pop else upper
    { vertices := lowerWithoutLast ++ upperWithoutLast }
where
  /-- Remove points from the end of hull that don't make a left turn with p. -/
  removeNonLeftTurns (hull : Array Vec2) (p : Vec2) : Array Vec2 :=
    if hull.size < 2 then hull
    else
      let last := hull[hull.size - 1]!
      let secondLast := hull[hull.size - 2]!
      let cross := (last - secondLast).cross (p - last)
      if cross <= 0.0 then
        removeNonLeftTurns hull.pop p
      else
        hull

/-- Compute convex hull of the polygon's vertices. -/
def toConvexHull (p : Polygon2D) : Polygon2D :=
  convexHull p.vertices

-- ============================================================================
-- Triangulation
-- ============================================================================

/-- A triangle represented by three vertex indices. -/
structure TriangleIndices where
  i0 : Nat
  i1 : Nat
  i2 : Nat
deriving Repr, BEq, Inhabited

/-- Check if vertex at index i is an ear (can be clipped).
    An ear is a convex vertex whose triangle doesn't contain any other vertices. -/
private def isEar (vertices : Array Vec2) (indices : Array Nat) (i : Nat) : Bool :=
  if indices.size < 3 then false
  else
    let n := indices.size
    let prevIdx := (i + n - 1) % n
    let nextIdx := (i + 1) % n

    let prev := vertices[indices[prevIdx]!]!
    let curr := vertices[indices[i]!]!
    let next := vertices[indices[nextIdx]!]!

    -- Check if the vertex is convex (makes a left turn)
    let cross := (curr - prev).cross (next - curr)
    if cross <= 0.0 then false  -- Not convex
    else
      -- Check if any other vertex is inside this triangle
      let tri := Polygon2D.triangle prev curr next
      -- Get the actual vertex indices we're using for the triangle
      let prevVertIdx := indices[prevIdx]!
      let currVertIdx := indices[i]!
      let nextVertIdx := indices[nextIdx]!
      -- j iterates over vertex indices (values in indices array)
      indices.foldl (init := true) fun canClip j =>
        if !canClip then false
        else if j == prevVertIdx || j == currVertIdx || j == nextVertIdx then true
        else
          let pt := vertices[j]!
          !tri.containsPoint pt

/-- Remove element at index from array. -/
private def eraseAt (arr : Array Nat) (idx : Nat) : Array Nat :=
  arr.foldl (init := (#[], 0)) (fun (result, i) elem =>
    if i == idx then (result, i + 1)
    else (result.push elem, i + 1)
  ) |>.1

/-- Find an ear in the polygon for triangulation. -/
private partial def findEarForTriangulation (vertices : Array Vec2) (indices : Array Nat) (start : Nat) : Option Nat :=
  if start >= indices.size then none
  else if isEar vertices indices start then some start
  else findEarForTriangulation vertices indices (start + 1)

/-- Main triangulation loop. -/
private partial def triangulateLoop (vertices : Array Vec2) (indices : Array Nat) (triangles : Array TriangleIndices) : Array TriangleIndices :=
  if indices.size < 3 then triangles
  else if indices.size == 3 then
    triangles.push ⟨indices[0]!, indices[1]!, indices[2]!⟩
  else
    -- Find an ear to clip
    match findEarForTriangulation vertices indices 0 with
    | none => triangles  -- Couldn't find ear (degenerate polygon)
    | some earIdx =>
      let n := indices.size
      let prevIdx := (earIdx + n - 1) % n
      let nextIdx := (earIdx + 1) % n
      let newTriangle := ⟨indices[prevIdx]!, indices[earIdx]!, indices[nextIdx]!⟩
      let newIndices := eraseAt indices earIdx
      triangulateLoop vertices newIndices (triangles.push newTriangle)

/-- Triangulate a simple polygon using the ear clipping algorithm.
    Returns an array of triangles represented by vertex indices.
    Works for convex and concave polygons, but not for self-intersecting ones. -/
def triangulate (p : Polygon2D) : Array TriangleIndices :=
  if p.vertices.size < 3 then #[]
  else if p.vertices.size == 3 then #[⟨0, 1, 2⟩]
  else
    -- Ensure counter-clockwise winding
    let poly := p.makeCounterClockwise
    -- Create index list
    let indices := (List.range poly.vertices.size).toArray
    triangulateLoop poly.vertices indices #[]

/-- Triangulate and return actual triangle polygons. -/
def triangulateToPolygons (p : Polygon2D) : Array Polygon2D :=
  let triangleIndices := p.triangulate
  triangleIndices.map fun tri =>
    Polygon2D.triangle p.vertices[tri.i0]! p.vertices[tri.i1]! p.vertices[tri.i2]!

/-- Get the triangles as arrays of Vec2 (each inner array has 3 vertices). -/
def triangulateToVertices (p : Polygon2D) : Array (Array Vec2) :=
  let triangleIndices := p.triangulate
  triangleIndices.map fun tri =>
    #[p.vertices[tri.i0]!, p.vertices[tri.i1]!, p.vertices[tri.i2]!]

-- ============================================================================
-- Triangulation with Holes
-- ============================================================================

/-- Polygon with holes (outer boundary plus zero or more holes). -/
structure WithHoles where
  outer : Polygon2D
  holes : Array Polygon2D
deriving Repr, Inhabited

namespace WithHoles

private def normalize (p : WithHoles) : WithHoles :=
  { outer := p.outer.makeCounterClockwise
    holes := p.holes.map (·.makeClockwise) }

private def rightmostIndex (verts : Array Vec2) : Nat := Id.run do
  if verts.isEmpty then
    return 0
  let mut best := 0
  for i in [1:verts.size] do
    let v := verts[i]!
    let b := verts[best]!
    if v.x > b.x || (Float.abs' (v.x - b.x) < Float.epsilon && v.y < b.y) then
      best := i
  return best

private def rotate (verts : Array Vec2) (start : Nat) : Array Vec2 :=
  if verts.isEmpty then verts
  else
    let idx := start % verts.size
    let a := verts.extract idx verts.size
    let b := verts.extract 0 idx
    a ++ b

private def insertAt (verts : Array Vec2) (idx : Nat) (v : Vec2) : Array Vec2 :=
  let front := verts.extract 0 idx
  let back := verts.extract idx verts.size
  front ++ #[v] ++ back

private def findBridgeIntersection (outerVerts : Array Vec2) (holePt : Vec2) : Option (Nat × Vec2) := Id.run do
  if outerVerts.size < 2 then
    return none
  let mut bestX := Float.infinity
  let mut bestEdge : Option Nat := none
  let mut bestPt := Vec2.zero
  for i in [:outerVerts.size] do
    let a := outerVerts[i]!
    let b := outerVerts[(i + 1) % outerVerts.size]!
    let dy := b.y - a.y
    if Float.abs' dy < Float.epsilon then
      continue
    let yMin := Float.min a.y b.y
    let yMax := Float.max a.y b.y
    if holePt.y > yMin && holePt.y <= yMax then
      let t := (holePt.y - a.y) / dy
      let x := a.x + t * (b.x - a.x)
      if x > holePt.x + Float.epsilon && x < bestX then
        bestX := x
        bestEdge := some i
        bestPt := Vec2.mk x holePt.y
  match bestEdge with
  | none => return none
  | some idx => return some (idx, bestPt)

private def spliceHole (outerVerts holeVerts : Array Vec2) : Array Vec2 := Id.run do
  if outerVerts.isEmpty || holeVerts.isEmpty then
    return outerVerts
  let holeStart := rightmostIndex holeVerts
  let holeRot := rotate holeVerts holeStart
  let holePt := holeRot[0]!
  match findBridgeIntersection outerVerts holePt with
  | none => return outerVerts
  | some (edgeIdx, bridgePt) =>
    let outerInserted := insertAt outerVerts (edgeIdx + 1) bridgePt
    let insertIdx := edgeIdx + 1
    let front := outerInserted.extract 0 (insertIdx + 1)
    let back := outerInserted.extract (insertIdx + 1) outerInserted.size
    let holeLoop := holeRot ++ #[holeRot[0]!]
    return front ++ holeLoop ++ #[bridgePt] ++ back

/-- Merge holes into a single simple polygon by adding bridge edges. -/
def mergeHoles (p : WithHoles) : Polygon2D := Id.run do
  let normalized := normalize p
  let mut verts := normalized.outer.vertices
  for hole in normalized.holes do
    verts := spliceHole verts hole.vertices
  return { vertices := verts }

/-- Triangulate a polygon with holes by bridging holes and ear clipping. -/
def triangulate (p : WithHoles) : Array TriangleIndices :=
  let merged := p.mergeHoles
  merged.triangulate

/-- Triangulate and return actual triangle polygons. -/
def triangulateToPolygons (p : WithHoles) : Array Polygon2D :=
  let merged := p.mergeHoles
  merged.triangulateToPolygons

/-- Get the triangles as arrays of Vec2 (each inner array has 3 vertices). -/
def triangulateToVertices (p : WithHoles) : Array (Array Vec2) :=
  let merged := p.mergeHoles
  merged.triangulateToVertices

end WithHoles

end Polygon2D

end Linalg
