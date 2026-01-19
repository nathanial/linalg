/-
  Delaunay Triangulation - Port of Mapbox Delaunator

  Computes the Delaunay triangulation of 2D points using an incremental algorithm
  with half-edge data structure for efficient traversal.

  The half-edge data structure:
  - For triangle t (0-indexed): half-edges are 3*t, 3*t+1, 3*t+2
  - For half-edge e:
    - Triangle index: e / 3
    - Next half-edge in triangle: 3*(e/3) + (e+1)%3
    - Previous half-edge: 3*(e/3) + (e+2)%3
    - Twin half-edge: halfedges[e] (none if on hull)
-/

import Linalg.Vec2
import Linalg.Geometry.AABB2D

namespace Linalg

namespace Delaunay

/-- Tolerance for floating point comparisons -/
private def epsilon : Float := 1e-12

/-- Tolerance for near-duplicate point comparisons (matches JS EPSILON). -/
private def duplicateEpsilon : Float := Float.pow 2.0 (-52.0)

/-- Result of Delaunay triangulation using half-edge data structure -/
structure Triangulation where
  /-- Original input points -/
  points : Array Vec2
  /-- Vertex indices in groups of 3 (each group is a triangle) -/
  triangles : Array Nat
  /-- Twin half-edge index for each half-edge (none if on convex hull) -/
  halfedges : Array (Option Nat)
  /-- Convex hull vertex indices in counter-clockwise order -/
  hull : Array Nat
deriving Repr, Inhabited

namespace Triangulation

/-- Number of triangles in the triangulation -/
def triangleCount (t : Triangulation) : Nat := t.triangles.size / 3

/-- Get the vertex indices for triangle at index i -/
def getTriangle (t : Triangulation) (i : Nat) : Option (Nat × Nat × Nat) :=
  let base := i * 3
  if base + 2 < t.triangles.size then
    some (t.triangles[base]!, t.triangles[base + 1]!, t.triangles[base + 2]!)
  else
    none

/-- Get the half-edge index for the ith edge of a triangle -/
def halfEdge (triangleIdx edgeIdx : Nat) : Nat := triangleIdx * 3 + edgeIdx

/-- Get triangle index from half-edge -/
def triangleOfEdge (e : Nat) : Nat := e / 3

/-- Get next half-edge within the same triangle -/
def nextHalfEdge (e : Nat) : Nat :=
  let t := e / 3
  t * 3 + (e + 1) % 3

/-- Get previous half-edge within the same triangle -/
def prevHalfEdge (e : Nat) : Nat :=
  let t := e / 3
  t * 3 + (e + 2) % 3

end Triangulation

-- ============================================================================
-- Geometric Predicates
-- ============================================================================

/-- Orientation test: positive if c is right of a->b, negative if left, zero if collinear. -/
def orient2d (a b c : Vec2) : Float :=
  (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y)

/-- In-circle test: negative if d is strictly inside the circumcircle of a,b,c (CCW ordered).
    Returns positive if outside, zero if on circle. -/
def inCircle (a b c d : Vec2) : Float :=
  let dx := a.x - d.x
  let dy := a.y - d.y
  let ex := b.x - d.x
  let ey := b.y - d.y
  let fx := c.x - d.x
  let fy := c.y - d.y

  let ap := dx * dx + dy * dy
  let bp := ex * ex + ey * ey
  let cp := fx * fx + fy * fy

  dx * (ey * cp - bp * fy) -
  dy * (ex * cp - bp * fx) +
  ap * (ex * fy - ey * fx)

/-- Circumcenter of triangle a,b,c. Returns none if degenerate (collinear). -/
def circumcenter (a b c : Vec2) : Option Vec2 :=
  let dx := b.x - a.x
  let dy := b.y - a.y
  let ex := c.x - a.x
  let ey := c.y - a.y

  let bl := dx * dx + dy * dy
  let cl := ex * ex + ey * ey
  let d := 2.0 * (dx * ey - dy * ex)

  if Float.abs d < epsilon then none
  else
    let x := a.x + (ey * bl - dy * cl) / d
    let y := a.y + (dx * cl - ex * bl) / d
    some (Vec2.mk x y)

/-- Circumradius squared of triangle a,b,c -/
def circumradiusSq (a b c : Vec2) : Float :=
  let dx := b.x - a.x
  let dy := b.y - a.y
  let ex := c.x - a.x
  let ey := c.y - a.y

  let bl := dx * dx + dy * dy
  let cl := ex * ex + ey * ey
  let d := 2.0 * (dx * ey - dy * ex)

  if Float.abs d < epsilon then Float.infinity
  else
    let x := (ey * bl - dy * cl) / d
    let y := (dx * cl - ex * bl) / d
    x * x + y * y

/-- Fast pseudo-angle for sorting (monotonic in actual angle, but not linear).
    Maps angle to [0, 4) range. -/
def pseudoAngle (dx dy : Float) : Float :=
  let p := dx / (Float.abs dx + Float.abs dy)
  -- Maps to [0, 1) for positive y, [3, 4) for negative y
  if dy > 0.0 then 3.0 - p
  else 1.0 + p

-- ============================================================================
-- Internal State for Triangulation Algorithm
-- ============================================================================

private structure TriState where
  /-- Triangle vertex indices (groups of 3) -/
  triangles : Array Nat
  /-- Twin half-edge for each half-edge -/
  halfedges : Array (Option Nat)
  /-- Convex hull as linked list: hull[i] = next vertex on hull after vertex i -/
  hullNext : Array Nat
  /-- Previous vertex on hull -/
  hullPrev : Array Nat
  /-- Triangles adjacent to hull: hullTri[i] = starting half-edge for vertex i -/
  hullTri : Array Nat
  /-- Hash table for hull vertices -/
  hullHash : Array (Option Nat)
  /-- Starting vertex of hull -/
  hullStart : Nat
  /-- Number of triangles created -/
  trianglesLen : Nat
  /-- Hash table size -/
  hashSize : Nat

private def TriState.empty (n : Nat) : TriState :=
  let maxTriangles := if n > 2 then 2 * n - 5 else 0
  let hashSize := Float.sqrt n.toFloat |>.ceil.toUInt32.toNat
  {
    triangles := List.replicate (maxTriangles * 3) 0 |>.toArray
    halfedges := List.replicate (maxTriangles * 3) none |>.toArray
    hullNext := List.replicate n 0 |>.toArray
    hullPrev := List.replicate n 0 |>.toArray
    hullTri := List.replicate n 0 |>.toArray
    hullHash := List.replicate hashSize none |>.toArray
    hullStart := 0
    trianglesLen := 0
    hashSize := hashSize
  }

/-- Hash a point to a hull hash table index -/
private def hashKey (x y : Float) (cx cy : Float) (hashSize : Nat) : Nat :=
  let angle := pseudoAngle (x - cx) (y - cy)
  ((angle / 4.0 * hashSize.toFloat).floor.toUInt64.toNat) % hashSize

/-- Link two half-edges to each other. -/
private def link (state : TriState) (a : Nat) (b : Option Nat) : TriState :=
  let halfedges := state.halfedges.set! a b
  let halfedges := match b with
    | some b' => halfedges.set! b' (some a)
    | none => halfedges
  { state with halfedges := halfedges }

/-- Add a triangle and return its first half-edge index -/
private def addTriangle (state : TriState) (i0 i1 i2 : Nat)
    (a b c : Option Nat) : TriState × Nat :=
  let t := state.trianglesLen
  let triangles := state.triangles
    |>.set! t i0
    |>.set! (t + 1) i1
    |>.set! (t + 2) i2

  let halfedges := state.halfedges
    |>.set! t a
    |>.set! (t + 1) b
    |>.set! (t + 2) c

  -- Link opposite half-edges
  let halfedges := match a with
    | some ae => halfedges.set! ae (some t)
    | none => halfedges

  let halfedges := match b with
    | some be => halfedges.set! be (some (t + 1))
    | none => halfedges

  let halfedges := match c with
    | some ce => halfedges.set! ce (some (t + 2))
    | none => halfedges

  ({ state with
     triangles := triangles
     halfedges := halfedges
     trianglesLen := t + 3 }, t)

/-- Legalize a half-edge by flipping if the Delaunay condition is violated.
    Uses stack-based iteration instead of recursion. -/
private def legalize (state : TriState) (points : Array Vec2) (startEdge : Nat) : TriState × Nat := Id.run do
  let mut st := state
  let mut stack : Array Nat := #[]
  let mut a := startEdge
  let mut ar := startEdge

  while true do
    let b := st.halfedges[a]!
    let a0 := a - a % 3
    ar := a0 + (a + 2) % 3

    match b with
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

      let illegal := inCircle pt0 ptr ptl pt1 < 0.0

      if illegal then
        -- Flip the edge
        st := { st with triangles := st.triangles.set! a p1 |>.set! b p0 }

        let hbl := st.halfedges[bl]!
        let har := st.halfedges[ar]!

        -- Edge swapped on the other side of the hull (rare); fix the half-edge reference
        if hbl.isNone then
          let mut e := st.hullStart
          let mut count := 0
          let mut updated := false
          while count < st.hullPrev.size && !updated do
            if st.hullTri[e]! == bl then
              st := { st with hullTri := st.hullTri.set! e a }
              updated := true
            e := st.hullPrev[e]!
            if e == st.hullStart then
              count := st.hullPrev.size
            count := count + 1

        st := link st a hbl
        st := link st b har
        st := link st ar (some bl)

        let br := b0 + (b + 1) % 3
        stack := stack.push br
      else
        if stack.size == 0 then
          break
        a := stack[stack.size - 1]!
        stack := stack.pop

  return (st, ar)

-- ============================================================================
-- Main Triangulation Function
-- ============================================================================

/-- Compute Delaunay triangulation of the given points.
    Returns none if fewer than 3 points or all points are collinear. -/
def triangulate (points : Array Vec2) : Option Triangulation :=
  if points.size < 3 then none
  else Id.run do
    let n := points.size

    -- Find bounding box
    let mut minX := Float.infinity
    let mut minY := Float.infinity
    let mut maxX := Float.negInfinity
    let mut maxY := Float.negInfinity

    for p in points do
      minX := Float.min minX p.x
      minY := Float.min minY p.y
      maxX := Float.max maxX p.x
      maxY := Float.max maxY p.y

    -- Center of bounding box
    let cx := (minX + maxX) / 2.0
    let cy := (minY + maxY) / 2.0

    -- Find seed point closest to center
    let mut i0 := 0
    let mut minDist := Float.infinity
    for i in [:n] do
      let d := (points[i]!.x - cx) * (points[i]!.x - cx) +
               (points[i]!.y - cy) * (points[i]!.y - cy)
      if d < minDist then
        i0 := i
        minDist := d

    let p0 := points[i0]!

    -- Find point closest to seed
    let mut i1 := 0
    minDist := Float.infinity
    for i in [:n] do
      if i != i0 then
        let d := p0.distanceSquared points[i]!
        if d < minDist && d > 0.0 then
          i1 := i
          minDist := d

    if i0 == i1 then return none

    let mut p1 := points[i1]!

    -- Find point that creates smallest circumcircle with i0, i1
    let mut i2 := 0
    let mut minRadius := Float.infinity
    for i in [:n] do
      if i != i0 && i != i1 then
        let r := circumradiusSq p0 p1 points[i]!
        if r < minRadius then
          i2 := i
          minRadius := r

    if minRadius == Float.infinity then return none

    let mut p2 := points[i2]!

    -- Ensure counter-clockwise orientation
    if orient2d p0 p1 p2 < 0.0 then
      let tmp := i1; i1 := i2; i2 := tmp
      let ptmp := p1; p1 := p2; p2 := ptmp

    -- Circumcenter of seed triangle
    let cc := circumcenter p0 p1 p2
    match cc with
    | none => return none
    | some center =>
      -- Sort points by distance from circumcenter
      let dists := points.mapIdx fun idx p =>
        (idx, (p.x - center.x) * (p.x - center.x) + (p.y - center.y) * (p.y - center.y))

      let sorted := dists.qsort fun (_, d1) (_, d2) => d1 < d2
      let ids := sorted.map (·.1)

      -- Initialize state
      let mut state := TriState.empty n

      -- Add first triangle
      state := { state with hullStart := i0 }

      -- Initialize hull linked lists
      state := { state with hullNext := state.hullNext.set! i0 i1 |>.set! i1 i2 |>.set! i2 i0 }
      state := { state with hullPrev := state.hullPrev.set! i0 i2 |>.set! i1 i0 |>.set! i2 i1 }

      let (state', _) := addTriangle state i0 i1 i2 none none none
      state := state'

      state := { state with hullTri := state.hullTri.set! i0 0 |>.set! i1 1 |>.set! i2 2 }

      -- Hash hull vertices
      let hashKey0 := hashKey p0.x p0.y cx cy state.hashSize
      let hashKey1 := hashKey p1.x p1.y cx cy state.hashSize
      let hashKey2 := hashKey p2.x p2.y cx cy state.hashSize
      state := { state with hullHash := state.hullHash.set! hashKey0 (some i0)
                                          |>.set! hashKey1 (some i1)
                                          |>.set! hashKey2 (some i2) }

      -- Incrementally add remaining points
      let mut xp := 0.0
      let mut yp := 0.0
      for k in [:ids.size] do
        let i := ids[k]!
        let p := points[i]!

        -- Skip near-duplicate points
        if k > 0 &&
            Float.abs (p.x - xp) <= duplicateEpsilon &&
            Float.abs (p.y - yp) <= duplicateEpsilon then
          continue
        xp := p.x
        yp := p.y

        -- Skip if this is one of the seed points
        if i == i0 || i == i1 || i == i2 then continue

        -- Find a visible edge on the convex hull using edge hash
        let key := hashKey p.x p.y cx cy state.hashSize
        let mut start := state.hullStart
        let mut foundStart := false
        let mut j := 0
        while j < state.hashSize && !foundStart do
          let idx := (key + j) % state.hashSize
          match state.hullHash[idx]! with
          | some s =>
            if s != state.hullNext[s]! then
              start := s
              foundStart := true
          | none => ()
          j := j + 1
        if !foundStart then
          start := state.hullStart

        start := state.hullPrev[start]!
        let startEdge := start
        let mut e := start
        let mut valid := true
        while true do
          let q := state.hullNext[e]!
          if orient2d p (points[e]!) (points[q]!) >= 0.0 then
            e := q
            if e == startEdge then
              valid := false
              break
          else
            break
        if !valid then continue

        -- Add the first triangle from the point
        let q := state.hullNext[e]!
        let (state', t) := addTriangle state e i q none none (some state.hullTri[e]!)
        state := state'
        let (state', legalEdge) := legalize state points (t + 2)
        state := state'
        state := { state with hullTri := state.hullTri.set! i legalEdge |>.set! e t }

        -- Walk forward through the hull, adding more triangles and flipping
        let mut nEdge := state.hullNext[e]!
        while true do
          let n_q := state.hullNext[nEdge]!
          if orient2d p (points[nEdge]!) (points[n_q]!) < 0.0 then
            let (state', nt) := addTriangle state nEdge i n_q (some state.hullTri[i]!) none (some state.hullTri[nEdge]!)
            state := state'
            let (state', legalEdge') := legalize state points (nt + 2)
            state := state'
            state := { state with hullTri := state.hullTri.set! i legalEdge' }
            state := { state with hullNext := state.hullNext.set! nEdge nEdge }
            nEdge := n_q
          else
            break

        -- Walk backward from the other side, adding more triangles and flipping
        if e == startEdge then
          while true do
            let q := state.hullPrev[e]!
            if orient2d p (points[q]!) (points[e]!) < 0.0 then
              let (state', nt) := addTriangle state q i e none (some state.hullTri[e]!) (some state.hullTri[q]!)
              state := state'
              let (state', _) := legalize state points (nt + 2)
              state := state'
              state := { state with hullTri := state.hullTri.set! q nt }
              state := { state with hullNext := state.hullNext.set! e e }
              e := q
            else
              break

        -- Update the hull indices
        state := { state with hullStart := e }
        state := { state with hullPrev := state.hullPrev.set! i e |>.set! nEdge i }
        state := { state with hullNext := state.hullNext.set! e i |>.set! i nEdge }

        -- Save the two new edges in the hash table
        let hashKeyI := hashKey p.x p.y cx cy state.hashSize
        state := { state with hullHash := state.hullHash.set! hashKeyI (some i) }
        let hashKeyE := hashKey (points[e]!.x) (points[e]!.y) cx cy state.hashSize
        state := { state with hullHash := state.hullHash.set! hashKeyE (some e) }

      -- Extract convex hull
      let mut hull : Array Nat := #[]
      let mut e := state.hullStart
      let mut hullCount := 0
      while hullCount < n do
        hull := hull.push e
        e := state.hullNext[e]!
        if e == state.hullStart then break
        hullCount := hullCount + 1

      -- Trim arrays to actual size
      let actualLen := state.trianglesLen
      let triangles := state.triangles.toList.take actualLen |>.toArray
      let halfedges := state.halfedges.toList.take actualLen |>.toArray

      return some {
        points := points
        triangles := triangles
        halfedges := halfedges
        hull := hull
      }

end Delaunay

end Linalg
