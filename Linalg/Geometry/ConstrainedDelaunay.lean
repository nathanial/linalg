/- 
  Constrained Delaunay Triangulation (polygons with holes)

  This module provides a constrained triangulation interface for polygons with
  holes by bridging holes and applying ear clipping. This produces triangles
  that respect the input boundaries and hole loops.
-/

import Std.Data.HashMap
import Std.Data.HashSet
import Linalg.Geometry.Polygon2D
import Linalg.Geometry.Line2D
import Linalg.Geometry.Delaunay

namespace Linalg

namespace ConstrainedDelaunay

/-- Result of constrained triangulation. -/
structure Triangulation where
  /-- Vertex positions used for triangulation. -/
  vertices : Array Vec2
  /-- Triangle indices in groups of 3. -/
  triangles : Array Nat
  /-- Twin half-edge index for each half-edge (none if on boundary). -/
  halfedges : Array (Option Nat)
deriving Repr, Inhabited

namespace Triangulation

/-- Number of triangles. -/
def triangleCount (t : Triangulation) : Nat := t.triangles.size / 3

/-- Get the half-edge index for the ith edge of a triangle. -/
def halfEdge (triangleIdx edgeIdx : Nat) : Nat := triangleIdx * 3 + edgeIdx

/-- Get triangle index from half-edge. -/
def triangleOfEdge (e : Nat) : Nat := e / 3

/-- Get next half-edge within the same triangle. -/
def nextHalfEdge (e : Nat) : Nat :=
  let t := e / 3
  t * 3 + (e + 1) % 3

/-- Get previous half-edge within the same triangle. -/
def prevHalfEdge (e : Nat) : Nat :=
  let t := e / 3
  t * 3 + (e + 2) % 3

/-- Get the three half-edges for a triangle. -/
def triangleEdges (t : Triangulation) (triangleIdx : Nat) : Option (Nat × Nat × Nat) :=
  let base := triangleIdx * 3
  if base + 2 < t.triangles.size then
    some (base, base + 1, base + 2)
  else
    none

/-- Number of half-edges. -/
def halfEdgeCount (t : Triangulation) : Nat := t.triangles.size

/-- Get the vertex index at half-edge e. -/
def edgeVertex (t : Triangulation) (e : Nat) : Option Nat :=
  if e < t.triangles.size then some t.triangles[e]! else none

/-- Get the start/end vertex indices for a half-edge. -/
def edgeVertices (t : Triangulation) (e : Nat) : Option (Nat × Nat) :=
  if e < t.triangles.size then
    let a := t.triangles[e]!
    let b := t.triangles[nextHalfEdge e]!
    some (a, b)
  else
    none

/-- Get the segment for a half-edge. -/
def edgeSegment (t : Triangulation) (e : Nat) : Option Segment2D :=
  match t.edgeVertices e with
  | none => none
  | some (a, b) =>
    if a < t.vertices.size && b < t.vertices.size then
      some { a := t.vertices[a]!, b := t.vertices[b]! }
    else
      none

/-- Get the length of a half-edge. -/
def edgeLength (t : Triangulation) (e : Nat) : Option Float :=
  match t.edgeSegment e with
  | none => none
  | some seg => some seg.length

/-- Get the twin half-edge if it exists. -/
def twinHalfEdge (t : Triangulation) (e : Nat) : Option Nat :=
  if e < t.halfedges.size then t.halfedges[e]! else none

/-- Check if a half-edge is on the boundary (no twin). -/
def isBoundaryEdge (t : Triangulation) (e : Nat) : Bool :=
  match t.twinHalfEdge e with
  | none => true
  | some _ => false

/-- Collect all boundary half-edges (no twin). -/
def boundaryEdges (t : Triangulation) : Array Nat := Id.run do
  let mut edges : Array Nat := #[]
  for e in [:t.halfedges.size] do
    if t.halfedges[e]!.isNone then
      edges := edges.push e
  return edges

private def buildBoundaryNext (t : Triangulation) (edges : Array Nat) : Std.HashMap Nat Nat := Id.run do
  let mut map : Std.HashMap Nat Nat := {}
  for e in edges do
    let a := t.triangles[e]!
    map := map.insert a e
  return map

/-- Walk boundary loops as arrays of vertex indices. -/
def boundaryLoops (t : Triangulation) : Array (Array Nat) := Id.run do
  let edges := t.boundaryEdges
  let mut visited : Std.HashSet Nat := {}
  let mut loops : Array (Array Nat) := #[]
  let nextMap := buildBoundaryNext t edges
  for start in edges do
    if visited.contains start then
      continue
    let mut loop : Array Nat := #[]
    let mut current := start
    let mut guard := 0
    while guard < edges.size do
      if visited.contains current then
        break
      visited := visited.insert current
      let a := t.triangles[current]!
      loop := loop.push a
      let b := t.triangles[nextHalfEdge current]!
      match Std.HashMap.get? nextMap b with
      | none => break
      | some next =>
        current := next
        if current == start then
          break
      guard := guard + 1
    if !loop.isEmpty then
      loops := loops.push loop
  return loops

/-- Get triangle indices for triangle at index i. -/
def getTriangle (t : Triangulation) (i : Nat) : Option (Nat × Nat × Nat) :=
  let base := i * 3
  if base + 2 < t.triangles.size then
    some (t.triangles[base]!, t.triangles[base + 1]!, t.triangles[base + 2]!)
  else
    none

/-- Get triangle vertices for triangle at index i. -/
def getTriangleVertices (t : Triangulation) (i : Nat) : Option (Vec2 × Vec2 × Vec2) :=
  match t.getTriangle i with
  | none => none
  | some (i0, i1, i2) =>
    if i0 < t.vertices.size && i1 < t.vertices.size && i2 < t.vertices.size then
      some (t.vertices[i0]!, t.vertices[i1]!, t.vertices[i2]!)
    else
      none

end Triangulation

private def constraintTolerance : Float := 1e-6

private def nextHalfEdge (e : Nat) : Nat :=
  let t := e / 3
  t * 3 + (e + 1) % 3

private def prevHalfEdge (e : Nat) : Nat :=
  let t := e / 3
  t * 3 + (e + 2) % 3

private def edgeKey (a b vertexCount : Nat) : Nat :=
  a * vertexCount + b

private def undirectedKey (a b vertexCount : Nat) : Nat :=
  if a < b then edgeKey a b vertexCount else edgeKey b a vertexCount

private def loopSegments (p : Polygon2D) : Array Segment2D := Id.run do
  let mut segs : Array Segment2D := #[]
  for i in [:p.vertices.size] do
    segs := segs.push { a := p.vertex i, b := p.vertex (i + 1) }
  return segs

private def boundarySegments (outer : Polygon2D) (holes : Array Polygon2D) : Array Segment2D :=
  holes.foldl (init := loopSegments outer) (fun acc hole => acc ++ loopSegments hole)

private def buildConstraintSet (merged : Polygon2D) (segments : Array Segment2D) : Std.HashSet Nat := Id.run do
  let n := merged.vertices.size
  let mut set : Std.HashSet Nat := {}
  if n < 2 then
    return set
  for i in [:n] do
    let a := merged.vertices[i]!
    let b := merged.vertices[(i + 1) % n]!
    let constrained := segments.any (fun s =>
      s.containsPoint a constraintTolerance && s.containsPoint b constraintTolerance)
    if constrained then
      let key := undirectedKey i ((i + 1) % n) n
      set := set.insert key
  return set

private def buildHalfEdges (triangles : Array Nat) (vertexCount : Nat) : Array (Option Nat) := Id.run do
  let edgeCount := triangles.size
  let mut map : Std.HashMap Nat Nat := {}
  for e in [:edgeCount] do
    let a := triangles[e]!
    let b := triangles[nextHalfEdge e]!
    map := map.insert (edgeKey a b vertexCount) e
  let mut halfedges := Array.replicate edgeCount none
  for e in [:edgeCount] do
    let a := triangles[e]!
    let b := triangles[nextHalfEdge e]!
    match Std.HashMap.get? map (edgeKey b a vertexCount) with
    | some t => halfedges := halfedges.set! e (some t)
    | none => ()
  return halfedges

private def isConstrainedEdge (triangles : Array Nat) (edge : Nat) (vertexCount : Nat) (constraints : Std.HashSet Nat) : Bool :=
  let a := triangles[edge]!
  let b := triangles[nextHalfEdge edge]!
  constraints.contains (undirectedKey a b vertexCount)

private def link (halfedges : Array (Option Nat)) (a : Nat) (b : Option Nat) : Array (Option Nat) :=
  let halfedges := halfedges.set! a b
  match b with
  | some b => halfedges.set! b (some a)
  | none => halfedges

structure TriState where
  triangles : Array Nat
  halfedges : Array (Option Nat)

private def isConvexQuad (points : Array Vec2) (pr pl p0 p1 : Nat) : Bool :=
  let a := points[pr]!
  let b := points[pl]!
  let o0 := Delaunay.orient2d a b points[p0]!
  let o1 := Delaunay.orient2d a b points[p1]!
  o0 * o1 < 0.0

private def legalize (state : TriState) (points : Array Vec2) (constraints : Std.HashSet Nat) (vertexCount : Nat) (startEdge : Nat) : TriState := Id.run do
  let mut st := state
  let mut stack : Array Nat := #[]
  let mut a := startEdge

  while true do
    if isConstrainedEdge st.triangles a vertexCount constraints then
      if stack.size == 0 then
        break
      a := stack[stack.size - 1]!
      stack := stack.pop
      continue

    let a0 := a - a % 3
    let ar := a0 + (a + 2) % 3

    match st.halfedges[a]! with
    | none =>
      if stack.size == 0 then
        break
      a := stack[stack.size - 1]!
      stack := stack.pop
    | some b =>
      let b0 := b - b % 3
      let al := a0 + (a + 1) % 3
      let bl := b0 + (b + 2) % 3

      let p0 := st.triangles[ar]!
      let pr := st.triangles[a]!
      let pl := st.triangles[al]!
      let p1 := st.triangles[bl]!

      let pt0 := points[p0]!
      let ptr := points[pr]!
      let ptl := points[pl]!
      let pt1 := points[p1]!

      let convex := isConvexQuad points pr pl p0 p1
      let illegal := convex && Delaunay.inCircle pt0 ptr ptl pt1 < 0.0

      if illegal then
        st := { st with triangles := st.triangles.set! a p1 |>.set! b p0 }

        let hbl := st.halfedges[bl]!
        let har := st.halfedges[ar]!

        let halfedges := link st.halfedges a hbl
        let halfedges := link halfedges b har
        let halfedges := link halfedges ar (some bl)

        st := { st with halfedges := halfedges }

        let br := b0 + (b + 1) % 3
        stack := stack.push br
      else
        if stack.size == 0 then
          break
        a := stack[stack.size - 1]!
        stack := stack.pop

  return st

private def legalizeAll (triangles : Array Nat) (halfedges : Array (Option Nat)) (points : Array Vec2) (constraints : Std.HashSet Nat) : TriState := Id.run do
  let vertexCount := points.size
  let mut st : TriState := { triangles := triangles, halfedges := halfedges }
  for e in [:triangles.size] do
    if !isConstrainedEdge st.triangles e vertexCount constraints then
      st := legalize st points constraints vertexCount e
  return st

private def flattenTriangles (tris : Array Polygon2D.TriangleIndices) : Array Nat :=
  tris.foldl (init := #[]) (fun acc tri =>
    acc.push tri.i0 |>.push tri.i1 |>.push tri.i2)

/-- Constrained triangulation for an outer polygon and optional holes.
    Uses hole bridging + ear clipping, ensuring triangles respect boundaries. -/
def triangulate (outer : Polygon2D) (holes : Array Polygon2D := #[]) : Triangulation :=
  let withHoles : Polygon2D.WithHoles := { outer := outer, holes := holes }
  let merged := withHoles.mergeHoles
  let tris := merged.triangulate
  let triangles := flattenTriangles tris
  let segments := boundarySegments outer holes
  let constraints := buildConstraintSet merged segments
  let halfedges := buildHalfEdges triangles merged.vertices.size
  let legalized := legalizeAll triangles halfedges merged.vertices constraints
  { vertices := merged.vertices, triangles := legalized.triangles, halfedges := legalized.halfedges }

/-- Constrained triangulation returning triangle polygons. -/
def triangulateToPolygons (outer : Polygon2D) (holes : Array Polygon2D := #[]) : Array Polygon2D :=
  Id.run do
    let tri := triangulate outer holes
    let triCount := tri.triangles.size / 3
    let mut polys : Array Polygon2D := #[]
    for i in [:triCount] do
      let base := i * 3
      let i0 := tri.triangles[base]!
      let i1 := tri.triangles[base + 1]!
      let i2 := tri.triangles[base + 2]!
      polys := polys.push (Polygon2D.triangle tri.vertices[i0]! tri.vertices[i1]! tri.vertices[i2]!)
    return polys

end ConstrainedDelaunay

end Linalg
