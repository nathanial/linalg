/-
  ConvexHull3D - Quickhull for 3D point clouds.

  Produces a triangulated convex hull for use in collision and physics tooling.
  Degenerate inputs (fewer than 4 non-coplanar points) return an empty hull.
-/

import Linalg.Vec3
import Linalg.Geometry.Triangle

namespace Linalg

/-- A triangulated 3D convex hull. -/
structure ConvexHull3D where
  /-- Original input points. -/
  points : Array Vec3
  /-- Triangle faces as point indices. -/
  faces : Array (Nat × Nat × Nat)
deriving Repr, Inhabited

namespace ConvexHull3D

/-- Number of faces in the hull. -/
def faceCount (h : ConvexHull3D) : Nat := h.faces.size

/-- Triangles as geometry. -/
def triangles (h : ConvexHull3D) : Array Triangle :=
  h.faces.map (fun (a, b, c) =>
    Triangle.mk' h.points[a]! h.points[b]! h.points[c]!)

/-- Unique vertex indices used by the hull. -/
def vertexIndices (h : ConvexHull3D) : Array Nat := Id.run do
  if h.points.isEmpty then
    return #[]
  let mut used := List.replicate h.points.size false |>.toArray
  for (a, b, c) in h.faces do
    used := used.set! a true
    used := used.set! b true
    used := used.set! c true
  let mut indices : Array Nat := #[]
  for i in [:used.size] do
    if used[i]! then
      indices := indices.push i
  return indices

-- ============================================================================
-- Quickhull Implementation
-- ============================================================================

private def epsilon : Float := 1e-6

private structure Face where
  a : Nat
  b : Nat
  c : Nat
  normal : Vec3
  offset : Float
  outside : Array Nat
deriving Repr, Inhabited

private def planeDistance (face : Face) (p : Vec3) : Float :=
  face.normal.dot p - face.offset

private def makeFace (points : Array Vec3) (a b c : Nat) (inside : Vec3) : Face :=
  let pa := points[a]!
  let pb := points[b]!
  let pc := points[c]!
  let n := (pb - pa).cross (pc - pa)
  if n.dot (inside - pa) > 0.0 then
    let n' := -n
    { a := a, b := c, c := b, normal := n', offset := n'.dot pa, outside := #[] }
  else
    { a := a, b := b, c := c, normal := n, offset := n.dot pa, outside := #[] }

private def isOutside (face : Face) (p : Vec3) : Bool :=
  planeDistance face p > epsilon

private def edgeKey (edge : Nat × Nat) : Nat × Nat :=
  if edge.1 <= edge.2 then edge else (edge.2, edge.1)

private def edgeEqualUnordered (e1 e2 : Nat × Nat) : Bool :=
  let k1 := edgeKey e1
  let k2 := edgeKey e2
  k1.1 == k2.1 && k1.2 == k2.2

private def findEdgeIdx (edges : Array (Nat × Nat)) (edge : Nat × Nat) : Option Nat := Id.run do
  let mut idx := 0
  for e in edges do
    if edgeEqualUnordered e edge then
      return some idx
    idx := idx + 1
  return none

private def eraseIdx {α : Type} (arr : Array α) (idx : Nat) : Array α :=
  (arr.extract 0 idx).append (arr.extract (idx + 1) arr.size)

private def addOrRemoveEdge (edges : Array (Nat × Nat)) (edge : Nat × Nat) : Array (Nat × Nat) :=
  match findEdgeIdx edges edge with
  | some idx => eraseIdx edges idx
  | none => edges.push edge

private def assignOutsidePoints (points : Array Vec3) (faces : Array Face) (candidates : Array Nat) : Array Face := Id.run do
  let mut updated := faces
  for idx in candidates do
    let p := points[idx]!
    let mut bestFace : Option Nat := none
    let mut bestDist := epsilon
    for fi in [:updated.size] do
      let dist := planeDistance updated[fi]! p
      if dist > bestDist then
        bestDist := dist
        bestFace := some fi
    match bestFace with
    | some fi =>
        let face := updated[fi]!
        updated := updated.set! fi { face with outside := face.outside.push idx }
    | none => ()
  return updated

private def findExtremesX (points : Array Vec3) : Nat × Nat := Id.run do
  let mut minIdx := 0
  let mut maxIdx := 0
  for i in [1:points.size] do
    let p := points[i]!
    if p.x < points[minIdx]!.x then
      minIdx := i
    if p.x > points[maxIdx]!.x then
      maxIdx := i
  return (minIdx, maxIdx)

private def farthestFromLine (points : Array Vec3) (aIdx bIdx : Nat) : Option Nat := Id.run do
  let a := points[aIdx]!
  let b := points[bIdx]!
  let ab := b - a
  if ab.lengthSquared <= epsilon then
    return none
  let mut bestIdx : Option Nat := none
  let mut bestDist := 0.0
  for i in [:points.size] do
    if i != aIdx && i != bIdx then
      let ap := points[i]! - a
      let dist := (ab.cross ap).lengthSquared
      if dist > bestDist then
        bestDist := dist
        bestIdx := some i
  if bestDist <= epsilon then none else bestIdx

private def farthestFromPlane (points : Array Vec3) (aIdx bIdx cIdx : Nat) : Option Nat := Id.run do
  let a := points[aIdx]!
  let b := points[bIdx]!
  let c := points[cIdx]!
  let n := (b - a).cross (c - a)
  if n.lengthSquared <= epsilon then
    return none
  let mut bestIdx : Option Nat := none
  let mut bestDist := 0.0
  for i in [:points.size] do
    if i != aIdx && i != bIdx && i != cIdx then
      let dist := Float.abs' (n.dot (points[i]! - a))
      if dist > bestDist then
        bestDist := dist
        bestIdx := some i
  if bestDist <= epsilon then none else bestIdx

private partial def buildHull (points : Array Vec3) (inside : Vec3) (faces : Array Face) : Array Face := Id.run do
  let mut current := faces
  let mut searching := true
  while searching do
    let mut bestFace : Option Nat := none
    let mut bestPoint := 0
    let mut bestDist := 0.0
    for fi in [:current.size] do
      let face := current[fi]!
      if !face.outside.isEmpty then
        for idx in face.outside do
          let dist := planeDistance face points[idx]!
          if dist > bestDist then
            bestDist := dist
            bestFace := some fi
            bestPoint := idx
    if bestFace.isNone then
      searching := false
    else
      let pIdx := bestPoint
      let p := points[pIdx]!
      let mut visibleMask := List.replicate current.size false |>.toArray
      let mut visibleFaces : Array Nat := #[]
      for fi in [:current.size] do
        if isOutside current[fi]! p then
          visibleMask := visibleMask.set! fi true
          visibleFaces := visibleFaces.push fi
      let mut horizon : Array (Nat × Nat) := #[]
      for fi in visibleFaces do
        let f := current[fi]!
        let edges : Array (Nat × Nat) := #[(f.a, f.b), (f.b, f.c), (f.c, f.a)]
        for e in edges do
          horizon := addOrRemoveEdge horizon e
      let mut reassigned : Array Nat := #[]
      for fi in visibleFaces do
        let f := current[fi]!
        for idx in f.outside do
          if idx != pIdx then
            reassigned := reassigned.push idx
      let mut nextFaces : Array Face := #[]
      for fi in [:current.size] do
        if !visibleMask[fi]! then
          nextFaces := nextFaces.push current[fi]!
      for e in horizon do
        nextFaces := nextFaces.push (makeFace points e.1 e.2 pIdx inside)
      current := assignOutsidePoints points nextFaces reassigned
  return current

/-- Compute the convex hull of a set of 3D points using Quickhull. -/
def quickHull (points : Array Vec3) : ConvexHull3D := Id.run do
  if points.size < 4 then
    return { points := points, faces := #[] }
  let (i0, i1) := findExtremesX points
  if i0 == i1 then
    return { points := points, faces := #[] }
  match farthestFromLine points i0 i1 with
  | none => return { points := points, faces := #[] }
  | some i2 =>
      match farthestFromPlane points i0 i1 i2 with
      | none => return { points := points, faces := #[] }
      | some i3 =>
          let p0 := points[i0]!
          let p1 := points[i1]!
          let p2 := points[i2]!
          let p3 := points[i3]!
          let inside := (p0 + p1 + p2 + p3) / 4.0
          let f0 := makeFace points i0 i1 i2 inside
          let f1 := makeFace points i0 i3 i1 inside
          let f2 := makeFace points i1 i3 i2 inside
          let f3 := makeFace points i2 i3 i0 inside
          let mut faces : Array Face := #[f0, f1, f2, f3]
          let mut remaining : Array Nat := #[]
          for i in [:points.size] do
            if i != i0 && i != i1 && i != i2 && i != i3 then
              remaining := remaining.push i
          faces := assignOutsidePoints points faces remaining
          let finalFaces := buildHull points inside faces
          let faceIndices := finalFaces.map (fun f => (f.a, f.b, f.c))
          return { points := points, faces := faceIndices }

end ConvexHull3D

end Linalg
