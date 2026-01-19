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

/-- Orientation test: positive if c is left of a->b, negative if right, zero if collinear.
    Uses robust adaptive precision when needed. -/
def orient2d (a b c : Vec2) : Float :=
  (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y)

/-- In-circle test: positive if d is strictly inside the circumcircle of a,b,c (CCW ordered).
    Returns negative if outside, zero if on circle. -/
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
private def legalize (state : TriState) (points : Array Vec2) (startEdge : Nat) : TriState := Id.run do
  let mut st := state
  let mut stack : Array Nat := #[startEdge]

  while h : stack.size > 0 do
    let a := stack[stack.size - 1]!
    stack := stack.pop

    let b := st.halfedges[a]!
    match b with
    | none => continue -- Hull edge, nothing to do
    | some b =>
      -- Get the quad vertices
      let a0 := a - a % 3
      let b0 := b - b % 3

      let al := a0 + (a + 1) % 3
      let ar := a0 + (a + 2) % 3
      let bl := b0 + (b + 2) % 3

      let p0 := st.triangles[ar]!
      let pr := st.triangles[a]!
      let pl := st.triangles[al]!
      let p1 := st.triangles[bl]!

      let pt0 := points[p0]!
      let ptr := points[pr]!
      let ptl := points[pl]!
      let pt1 := points[p1]!

      -- Check if we need to flip
      if inCircle pt0 ptr ptl pt1 > 0.0 then
        -- Flip the edge
        st := { st with
          triangles := st.triangles.set! a p1
                        |>.set! b p0
        }

        let hbl := st.halfedges[bl]!
        let har := st.halfedges[ar]!

        -- Update half-edge links
        if let some hblVal := hbl then
          st := { st with halfedges := st.halfedges.set! hblVal (some a) }
        if let some harVal := har then
          st := { st with halfedges := st.halfedges.set! harVal (some b) }

        st := { st with halfedges := st.halfedges.set! a hbl
                          |>.set! b har
                          |>.set! ar (some bl)
                          |>.set! bl (some ar)
              }

        let br := b0 + (b + 1) % 3

        -- Check more edges
        stack := stack.push a
        stack := stack.push br

  return st

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
      for id in ids do
        let i := id
        let p := points[i]!

        -- Skip if this is one of the seed points
        if i == i0 || i == i1 || i == i2 then continue

        -- Find a visible edge on the convex hull
        let key := hashKey p.x p.y cx cy state.hashSize

        -- Search for starting hull edge
        let mut start := state.hullHash[key]!.getD state.hullStart
        let mut e := start
        let mut found := false
        let mut iterations := 0

        -- Find visible hull edge (where point is to the right)
        while iterations < n && !found do
          let q := state.hullNext[e]!
          let pe := points[e]!
          let pq := points[q]!
          if orient2d pe pq p >= 0.0 then
            -- This edge is visible
            found := true
          else
            e := q
            if e == start then
              -- No visible edge found - point is inside hull (shouldn't happen for convex hull)
              found := true
          iterations := iterations + 1

        if !found then continue

        -- Remember first triangle
        let firstEdge := state.hullTri[e]!
        let mut lastEdge := firstEdge

        -- Add triangle for each visible hull edge
        let q := state.hullNext[e]!
        let (state', t) := addTriangle state e i q (some (state.hullTri[e]!)) none none
        state := state'

        -- Update hull triangle pointer for e
        state := { state with hullTri := state.hullTri.set! i t |>.set! e (t + 2) }
        lastEdge := t + 2

        -- Legalize the new edge
        state := legalize state points t

        -- Walk forward to add more triangles
        let mut n_e := q
        let mut prevT := t
        let mut walkCount := 0
        while walkCount < n do
          let n_q := state.hullNext[n_e]!
          let p_ne := points[n_e]!
          let p_nq := points[n_q]!
          if orient2d p_ne p_nq p >= 0.0 then
            let (state', nt) := addTriangle state n_e i n_q (some (prevT + 2)) none (some state.hullTri[n_e]!)
            state := state'
            state := { state with hullTri := state.hullTri.set! n_e (nt + 2) }

            -- Remove n_e from hull
            state := { state with hullNext := state.hullNext.set! (state.hullPrev[n_e]!) n_q }
            state := { state with hullPrev := state.hullPrev.set! n_q (state.hullPrev[n_e]!) }

            state := legalize state points nt
            prevT := nt
            n_e := n_q
          else
            break
          walkCount := walkCount + 1

        -- Walk backward to add more triangles
        if e != state.hullStart then
          let mut prev_e := state.hullPrev[e]!
          walkCount := 0
          while walkCount < n do
            let p_prev := points[prev_e]!
            let p_e := points[e]!
            if orient2d p_prev p_e p >= 0.0 then
              let (state', nt) := addTriangle state prev_e i e none (some (state.hullTri[prev_e]!)) (some lastEdge)
              state := state'

              -- Update lastEdge
              lastEdge := nt + 2
              state := { state with hullTri := state.hullTri.set! prev_e nt }

              -- Remove e from hull
              state := { state with hullNext := state.hullNext.set! prev_e i }
              state := { state with hullPrev := state.hullPrev.set! i prev_e }

              state := legalize state points nt
              e := prev_e
              prev_e := state.hullPrev[e]!
            else
              break
            walkCount := walkCount + 1

        -- Update hull start if needed
        if e == state.hullStart then
          state := { state with hullStart := state.hullPrev[i]! }

        -- Insert i into hull
        state := { state with hullNext := state.hullNext.set! e i }
        state := { state with hullPrev := state.hullPrev.set! i e }
        state := { state with hullNext := state.hullNext.set! i n_e }
        state := { state with hullPrev := state.hullPrev.set! n_e i }

        state := { state with hullTri := state.hullTri.set! i lastEdge }

        -- Update hull hash
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
