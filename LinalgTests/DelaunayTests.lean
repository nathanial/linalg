/-
  Tests for Delaunay triangulation.
-/

import Linalg
import Crucible

namespace LinalgTests.DelaunayTests

open Crucible
open Linalg

private def orient (a b c : Vec2) : Float :=
  (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)

private def polygonArea2 (points : Array Vec2) (hull : Array Nat) : Float := Id.run do
  if hull.isEmpty then
    return 0.0
  let mut area := 0.0
  for i in [:hull.size] do
    let j := (i + 1) % hull.size
    let pi := points[hull[i]!]!
    let pj := points[hull[j]!]!
    area := area + (pi.x * pj.y - pj.x * pi.y)
  return Float.abs area

private def trianglesArea2 (points : Array Vec2) (triangles : Array Nat) : Float := Id.run do
  let mut area := 0.0
  for t in [:triangles.size / 3] do
    let base := t * 3
    let a := points[triangles[base]!]!
    let b := points[triangles[base + 1]!]!
    let c := points[triangles[base + 2]!]!
    area := area + Float.abs (orient a b c)
  return area

private def findIdx? (arr : Array Nat) (value : Nat) : Option Nat := Id.run do
  let mut i := 0
  while i < arr.size do
    if arr[i]! == value then
      return some i
    i := i + 1
  return none

private def sameCycle (a b : Array Nat) : Bool := Id.run do
  if a.size != b.size then
    return false
  if a.isEmpty then
    return true
  match findIdx? b a[0]! with
  | none => return false
  | some start =>
    let mut ok := true
    for i in [:a.size] do
      if a[i]! != b[(start + i) % b.size]! then
        ok := false
    return ok

private def isPermutation (arr : Array Nat) (n : Nat) : Bool := Id.run do
  if arr.size != n then
    return false
  let mut seen : Array Bool := List.replicate n false |>.toArray
  for v in arr do
    if v >= n then
      return false
    if seen[v]! then
      return false
    seen := seen.set! v true
  return true

private def countHullEdges (halfedges : Array (Option Nat)) : Nat := Id.run do
  let mut count := 0
  for h in halfedges do
    if h.isNone then
      count := count + 1
  return count

private def regularPolygon (n : Nat) (radius : Float) : Array Vec2 := Id.run do
  let mut verts : Array Vec2 := #[]
  for i in [:n] do
    let angle := 2.0 * Float.pi * (i.toFloat / n.toFloat)
    verts := verts.push (Vec2.mk (radius * Float.cos angle) (radius * Float.sin angle))
  return verts

private def validateTriangulation (points : Array Vec2) (tri : Delaunay.Triangulation) : IO Unit := do
  ensure (tri.triangles.size == tri.halfedges.size) "triangles and halfedges sizes must match"
  ensure (tri.triangles.size % 3 == 0) "triangles length must be a multiple of 3"

  for e in [:tri.halfedges.size] do
    match tri.halfedges[e]! with
    | none => pure ()
    | some twin =>
      ensure (twin < tri.halfedges.size) "halfedge twin out of range"
      match tri.halfedges[twin]! with
      | some back => ensure (back == e) "halfedge connection mismatch"
      | none => ensure false "halfedge twin missing"

  for t in [:tri.triangles.size / 3] do
    let base := t * 3
    let i0 := tri.triangles[base]!
    let i1 := tri.triangles[base + 1]!
    let i2 := tri.triangles[base + 2]!
    ensure (i0 < points.size && i1 < points.size && i2 < points.size) "triangle index out of range"
    ensure (i0 != i1 && i1 != i2 && i0 != i2) "triangle has duplicate vertices"

  if tri.hull.size >= 3 then
    let mut sign := 0.0
    for i in [:tri.hull.size] do
      let prev := tri.hull[(i + tri.hull.size - 1) % tri.hull.size]!
      let curr := tri.hull[i]!
      let next := tri.hull[(i + 1) % tri.hull.size]!
      let o := orient points[prev]! points[curr]! points[next]!
      if Float.abs o > 1e-9 && sign == 0.0 then
        sign := if o > 0.0 then 1.0 else -1.0
      if sign != 0.0 then
        ensure (o * sign >= -1e-9) "hull should be convex"

  let hullArea := polygonArea2 points tri.hull
  if hullArea > 0.0 then
    let triArea := trianglesArea2 points tri.triangles
    let err := Float.abs (hullArea - triArea) / hullArea
    ensure (err <= 1e-6) "triangulation area should match hull area"

private def validateDelaunay (points : Array Vec2) (tri : Delaunay.Triangulation) : IO Unit := do
  let tol := 1e-9
  for t in [:tri.triangles.size / 3] do
    let base := t * 3
    let i0 := tri.triangles[base]!
    let i1 := tri.triangles[base + 1]!
    let i2 := tri.triangles[base + 2]!
    let a := points[i0]!
    let b := points[i1]!
    let c := points[i2]!
    let o := Delaunay.orient2d a b c
    ensure (Float.abs o > 1e-12) "triangle should be non-degenerate"
    let sign := if o < 0.0 then -1.0 else 1.0
    for pIdx in [:points.size] do
      if pIdx == i0 || pIdx == i1 || pIdx == i2 then continue
      let p := points[pIdx]!
      let val := Delaunay.inCircle a b c p * sign
      ensure (val >= -tol) "point should not be inside circumcircle"

testSuite "Delaunay"

test "square with center triangulates to star" := do
  let points := #[(Vec2.mk 0.0 0.0), (Vec2.mk 2.0 0.0), (Vec2.mk 2.0 2.0), (Vec2.mk 0.0 2.0), (Vec2.mk 1.0 1.1)]
  match Delaunay.triangulate points with
  | none => ensure false "expected triangulation"
  | some tri =>
    ensure (tri.triangles.size == 12) "expected 4 triangles"
    let mut triSign := 0.0
    for t in [:tri.triangles.size / 3] do
      let base := t * 3
      let i0 := tri.triangles[base]!
      let i1 := tri.triangles[base + 1]!
      let i2 := tri.triangles[base + 2]!
      ensure (i0 == 4 || i1 == 4 || i2 == 4) "triangle should include center"
      let a := points[i0]!
      let b := points[i1]!
      let c := points[i2]!
      let o := orient a b c
      ensure (Float.abs o > 1e-12) "triangle should be non-degenerate"
      if triSign == 0.0 then
        triSign := if o > 0.0 then 1.0 else -1.0
      ensure (o * triSign > 0.0) "triangle winding should be consistent"
    ensure (tri.hull.size == 4) "expected 4 hull points"
    let expected := #[0, 1, 2, 3]
    let expectedRev := #[0, 3, 2, 1]
    ensure (sameCycle tri.hull expected || sameCycle tri.hull expectedRev) "hull should match square"

test "triangulation passes validation" := do
  let points := #[(Vec2.mk 0.0 0.0), (Vec2.mk 1.0 0.2), (Vec2.mk 2.0 0.0), (Vec2.mk 2.5 1.2),
    (Vec2.mk 1.8 2.4), (Vec2.mk 0.6 2.0), (Vec2.mk (-0.4) 1.2), (Vec2.mk 0.4 1.0),
    (Vec2.mk 1.2 1.5), (Vec2.mk 1.6 0.8)]
  match Delaunay.triangulate points with
  | none => ensure false "expected triangulation"
  | some tri =>
    validateTriangulation points tri
    validateDelaunay points tri

test "convex hexagon triangulates to n-2 triangles" := do
  let points := regularPolygon 6 1.0
  match Delaunay.triangulate points with
  | none => ensure false "expected triangulation"
  | some tri =>
    ensure (tri.triangles.size == 12) "expected 4 triangles"
    ensure (tri.hull.size == 6) "expected 6 hull points"
    ensure (isPermutation tri.hull 6) "hull should contain all vertices"
    validateTriangulation points tri
    validateDelaunay points tri

test "duplicate points do not break triangulation" := do
  let points := #[
    (Vec2.mk 0.0 0.0),
    (Vec2.mk 1.0 0.0),
    (Vec2.mk 0.0 1.0),
    (Vec2.mk 0.0 1.0),
    (Vec2.mk 1.0 0.0)
  ]
  match Delaunay.triangulate points with
  | none => ensure false "expected triangulation"
  | some tri =>
    ensure (tri.triangles.size >= 3) "expected at least one triangle"
    validateTriangulation points tri
    validateDelaunay points tri

test "collinear points return none" := do
  let points := #[(Vec2.mk 0.0 0.0), (Vec2.mk 1.0 0.0), (Vec2.mk 2.0 0.0), (Vec2.mk 3.0 0.0)]
  match Delaunay.triangulate points with
  | none => pure ()
  | some _ => ensure false "expected none for collinear input"

test "less than three points return none" := do
  let points0 : Array Vec2 := #[]
  let points1 := #[(Vec2.mk 0.0 0.0)]
  let points2 := #[(Vec2.mk 0.0 0.0), (Vec2.mk 1.0 0.0)]
  ensure (Delaunay.triangulate points0).isNone "expected none for empty input"
  ensure (Delaunay.triangulate points1).isNone "expected none for single point"
  ensure (Delaunay.triangulate points2).isNone "expected none for two points"

test "triangle produces single triangle and hull" := do
  let points := #[(Vec2.mk 0.0 0.0), (Vec2.mk 2.0 0.0), (Vec2.mk 0.5 1.0)]
  match Delaunay.triangulate points with
  | none => ensure false "expected triangulation"
  | some tri =>
    ensure (tri.triangles.size == 3) "expected 1 triangle"
    ensure (tri.hull.size == 3) "expected 3 hull points"
    ensure (isPermutation tri.hull 3) "hull should contain all vertices"
    ensure (isPermutation tri.triangles 3) "triangle should contain all vertices"
    validateTriangulation points tri
    validateDelaunay points tri

test "triangle with interior point forms a fan" := do
  let points := #[(Vec2.mk 0.0 0.0), (Vec2.mk 2.0 0.0), (Vec2.mk 0.0 2.0), (Vec2.mk 0.6 0.6)]
  match Delaunay.triangulate points with
  | none => ensure false "expected triangulation"
  | some tri =>
    ensure (tri.triangles.size == 9) "expected 3 triangles"
    ensure (tri.hull.size == 3) "expected triangular hull"
    for t in [:tri.triangles.size / 3] do
      let base := t * 3
      let i0 := tri.triangles[base]!
      let i1 := tri.triangles[base + 1]!
      let i2 := tri.triangles[base + 2]!
      ensure (i0 == 3 || i1 == 3 || i2 == 3) "each triangle should include interior point"
    validateTriangulation points tri
    validateDelaunay points tri

test "convex polygon has hull-sized boundary" := do
  let points := regularPolygon 5 1.5
  match Delaunay.triangulate points with
  | none => ensure false "expected triangulation"
  | some tri =>
    ensure (tri.hull.size == 5) "expected 5 hull points"
    ensure (isPermutation tri.hull 5) "hull should contain all vertices"
    ensure (countHullEdges tri.halfedges == tri.hull.size) "boundary edges should match hull size"
    validateTriangulation points tri
    validateDelaunay points tri



end LinalgTests.DelaunayTests
