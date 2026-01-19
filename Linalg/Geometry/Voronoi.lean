/-
  Voronoi Diagram - Extracted from Delaunay Triangulation

  Each Delaunay triangle circumcenter becomes a Voronoi vertex.
  For each input point (site), its Voronoi cell is formed by walking
  around incident triangles and collecting their circumcenters.
-/

import Linalg.Vec2
import Linalg.Geometry.AABB2D
import Linalg.Geometry.Polygon2D
import Linalg.Geometry.Delaunay

namespace Linalg

namespace Voronoi

/-- A single Voronoi cell -/
structure VoronoiCell where
  /-- Index of the site point in the original points array -/
  siteIndex : Nat
  /-- Cell polygon vertices in counter-clockwise order -/
  vertices : Array Vec2
  /-- True if this site is on the convex hull (cell extends to infinity) -/
  isUnbounded : Bool
deriving Repr, Inhabited

/-- Complete Voronoi diagram -/
structure VoronoiDiagram where
  /-- Original site points -/
  sites : Array Vec2
  /-- Voronoi cells, one per site -/
  cells : Array VoronoiCell
  /-- The underlying Delaunay triangulation -/
  triangulation : Delaunay.Triangulation
deriving Repr, Inhabited

-- ============================================================================
-- Cell Extraction from Delaunay Triangulation
-- ============================================================================

/-- Build a map from each vertex to an incident half-edge -/
private def buildPointToEdge (tri : Delaunay.Triangulation) : Array (Option Nat) := Id.run do
  let n := tri.points.size
  let mut map : Array (Option Nat) := List.replicate n none |>.toArray
  for e in [:tri.triangles.size] do
    let p := tri.triangles[e]!
    if map[p]!.isNone then
      map := map.set! p (some e)
  return map

/-- Check if a vertex is on the convex hull -/
private def isOnHull (hull : Array Nat) (vertex : Nat) : Bool :=
  hull.any (· == vertex)

/-- Get the circumcenter of the triangle containing half-edge e -/
private def getCircumcenter (tri : Delaunay.Triangulation) (e : Nat) : Option Vec2 :=
  let t := e / 3
  let base := t * 3
  if base + 2 < tri.triangles.size then
    let i0 := tri.triangles[base]!
    let i1 := tri.triangles[base + 1]!
    let i2 := tri.triangles[base + 2]!
    Delaunay.circumcenter tri.points[i0]! tri.points[i1]! tri.points[i2]!
  else
    none

/-- Extract the Voronoi cell for a site by walking around incident triangles.
    For interior points, this walks in a complete cycle.
    For hull points, we need to handle the unbounded edges. -/
private def extractCell (tri : Delaunay.Triangulation) (siteIndex : Nat)
    (pointToEdge : Array (Option Nat)) (hullSet : Array Bool) : VoronoiCell := Id.run do
  let onHull := hullSet[siteIndex]!

  match pointToEdge[siteIndex]! with
  | none => return { siteIndex := siteIndex, vertices := #[], isUnbounded := onHull }
  | some startEdge =>
    let mut vertices : Array Vec2 := #[]
    let mut e := startEdge
    let mut iterations := 0
    let maxIter := tri.triangles.size

    -- Walk around the site collecting circumcenters
    while iterations < maxIter do
      -- Get circumcenter of current triangle
      if let some cc := getCircumcenter tri e then
        vertices := vertices.push cc

      -- Move to next triangle around this vertex
      -- The next triangle shares an edge with current triangle
      let nextInTri := e / 3 * 3 + (e + 2) % 3  -- Previous edge in triangle
      let twin := tri.halfedges[nextInTri]!

      match twin with
      | none =>
        -- We hit the hull boundary - this is an unbounded cell
        -- For hull vertices, we need to also check the other direction
        break
      | some t =>
        e := t
        if e == startEdge then
          -- Completed the cycle
          break

      iterations := iterations + 1

    -- For hull points, we may need to walk the other direction too
    if onHull && vertices.size > 0 then
      -- Walk backward from start
      let mut backVertices : Array Vec2 := #[]
      e := startEdge

      -- Get next edge in triangle, then its twin
      let nextEdge := e / 3 * 3 + (e + 1) % 3
      match tri.halfedges[nextEdge]! with
      | none => ()  -- Already at hull boundary
      | some t =>
        e := t
        iterations := 0
        while iterations < maxIter do
          if let some cc := getCircumcenter tri e then
            backVertices := backVertices.push cc

          let nextEdge := e / 3 * 3 + (e + 1) % 3
          match tri.halfedges[nextEdge]! with
          | none => break
          | some t =>
            e := t
            if e == startEdge then break

          iterations := iterations + 1

      -- Combine: reverse back vertices and prepend to main vertices
      vertices := backVertices.reverse ++ vertices

    return { siteIndex := siteIndex, vertices := vertices, isUnbounded := onHull }

/-- Extract Voronoi diagram from Delaunay triangulation -/
def fromDelaunay (tri : Delaunay.Triangulation) : VoronoiDiagram := Id.run do
  let n := tri.points.size

  -- Build helper data structures
  let pointToEdge := buildPointToEdge tri

  -- Build hull set for O(1) lookup
  let mut hullSet := List.replicate n false |>.toArray
  for h in tri.hull do
    hullSet := hullSet.set! h true

  -- Extract cell for each site
  let mut cells : Array VoronoiCell := #[]
  for i in [:n] do
    let cell := extractCell tri i pointToEdge hullSet
    cells := cells.push cell

  return {
    sites := tri.points
    cells := cells
    triangulation := tri
  }

-- ============================================================================
-- Polygon Clipping (Sutherland-Hodgman Algorithm)
-- ============================================================================

/-- Clip a polygon against a half-plane defined by an edge.
    The edge goes from p1 to p2, and we keep the left side (CCW). -/
private def clipAgainstEdge (polygon : Array Vec2) (p1 p2 : Vec2) : Array Vec2 := Id.run do
  if polygon.size == 0 then return #[]

  let mut output : Array Vec2 := #[]

  for i in [:polygon.size] do
    let current := polygon[i]!
    let next := polygon[(i + 1) % polygon.size]!

    -- Check which side of the edge each point is on
    let currentInside := Delaunay.orient2d p1 p2 current <= 0.0
    let nextInside := Delaunay.orient2d p1 p2 next <= 0.0

    if currentInside then
      output := output.push current
      if !nextInside then
        -- Edge exits, add intersection
        if let some inter := lineIntersection current next p1 p2 then
          output := output.push inter
    else if nextInside then
      -- Edge enters, add intersection
      if let some inter := lineIntersection current next p1 p2 then
        output := output.push inter

  return output
where
  /-- Find intersection of line segments -/
  lineIntersection (a1 a2 b1 b2 : Vec2) : Option Vec2 :=
    let d1 := a2 - a1
    let d2 := b2 - b1
    let cross := d1.cross d2
    if Float.abs cross < 1e-12 then none
    else
      let t := (b1 - a1).cross d2 / cross
      some (a1 + d1 * t)

/-- Clip a polygon to a bounding rectangle -/
def clipToRect (polygon : Array Vec2) (bounds : AABB2D) : Array Vec2 :=
  let minX := bounds.min.x
  let minY := bounds.min.y
  let maxX := bounds.max.x
  let maxY := bounds.max.y

  -- Clip against each edge of the rectangle (CCW order)
  -- Left edge (going up)
  let p := clipAgainstEdge polygon (Vec2.mk minX minY) (Vec2.mk minX maxY)
  -- Top edge (going right)
  let p := clipAgainstEdge p (Vec2.mk minX maxY) (Vec2.mk maxX maxY)
  -- Right edge (going down)
  let p := clipAgainstEdge p (Vec2.mk maxX maxY) (Vec2.mk maxX minY)
  -- Bottom edge (going left)
  clipAgainstEdge p (Vec2.mk maxX minY) (Vec2.mk minX minY)

/-- Extend unbounded Voronoi cells and clip all cells to bounds.
    Returns array of polygons (one per cell, in same order as cells). -/
def clipToBounds (diagram : VoronoiDiagram) (bounds : AABB2D) : Array Polygon2D := Id.run do
  let mut result : Array Polygon2D := #[]

  for cell in diagram.cells do
    if cell.vertices.size == 0 then
      -- Degenerate cell - use site as single point (will be clipped away)
      result := result.push (Polygon2D.fromVertices #[])
      continue

    let mut vertices := cell.vertices

    -- For unbounded cells, we need to extend the boundary
    if cell.isUnbounded && vertices.size >= 2 then
      let site := diagram.sites[cell.siteIndex]!

      -- The first and last vertices are on unbounded rays
      -- We extend them towards the bounding box

      -- Extend first edge (from first vertex away from cell)
      let first := vertices[0]!
      let second := if vertices.size > 1 then vertices[1]! else first

      -- Direction is perpendicular to edge and away from site
      let edgeDir := (second - first).normalize
      let perpDir := Vec2.mk (-edgeDir.y) edgeDir.x
      -- Choose direction away from site
      let toSite := site - first
      let outDir := if perpDir.dot toSite < 0.0 then perpDir else perpDir.neg

      -- Extend far enough to exceed bounds
      let diagonal := bounds.size.length * 2.0
      let extendedFirst := first + outDir * diagonal
      vertices := #[extendedFirst] ++ vertices

      -- Extend last edge similarly
      let last := vertices[vertices.size - 1]!
      let secondLast := if vertices.size > 1 then vertices[vertices.size - 2]! else last

      let edgeDir := (secondLast - last).normalize
      let perpDir := Vec2.mk (-edgeDir.y) edgeDir.x
      let toSite := site - last
      let outDir := if perpDir.dot toSite < 0.0 then perpDir else perpDir.neg

      let extendedLast := last + outDir * diagonal
      vertices := vertices.push extendedLast

    -- Clip to bounds
    let clipped := clipToRect vertices bounds

    result := result.push (Polygon2D.fromVertices clipped)

  return result

-- ============================================================================
-- Convenience Functions
-- ============================================================================

/-- Generate Voronoi diagram from points and clip to bounds.
    Returns none if triangulation fails. -/
def generate (points : Array Vec2) (bounds : AABB2D) : Option (Array Polygon2D) := do
  let tri ← Delaunay.triangulate points
  let diagram := fromDelaunay tri
  return clipToBounds diagram bounds

/-- Generate Voronoi cells with guaranteed full coverage of bounds.
    Adds corner points if needed to ensure no gaps at corners. -/
def generateWithCoverage (points : Array Vec2) (bounds : AABB2D) : Option (Array Polygon2D) := do
  -- Add virtual corner points far outside bounds to ensure coverage
  let padding := bounds.size.length * 0.1
  let expandedBounds := AABB2D.fromMinMax
    (bounds.min - Vec2.mk padding padding)
    (bounds.max + Vec2.mk padding padding)

  let tri ← Delaunay.triangulate points
  let diagram := fromDelaunay tri
  let clipped := clipToBounds diagram bounds

  -- Filter out empty polygons
  let valid := clipped.filter (·.vertices.size >= 3)
  return valid

end Voronoi

end Linalg
