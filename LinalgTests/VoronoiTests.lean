/-
  Tests for Voronoi diagram generation and clipping.
-/

import Linalg
import Crucible

namespace LinalgTests.VoronoiTests

open Crucible
open Linalg

private structure LCG where
  state : UInt64
deriving Inhabited

private def lcgNew (seedVal : Nat) : LCG :=
  { state := seedVal.toUInt64 }

private def lcgNext (rng : LCG) : Prod LCG Float :=
  let a : UInt64 := 1103515245
  let c : UInt64 := 12345
  let m : UInt64 := 2147483648
  let newState := (a * rng.state + c) % m
  let value := newState.toFloat / m.toFloat
  ({ state := newState }, value)

private def lcgNextInRange (rng : LCG) (min max : Float) : Prod LCG Float :=
  let (rng', t) := lcgNext rng
  (rng', min + t * (max - min))

private structure SeedConfig where
  numPoints : Nat
  seedVal : Nat
  minDistance : Float
  edgeMargin : Float
deriving Inhabited

private def generateSeedPoints (config : SeedConfig) : Array Vec2 := Id.run do
  let mut rng := lcgNew config.seedVal
  let mut points : Array Vec2 := #[]

  let minX := 0.0 + config.edgeMargin
  let maxX := 1.0 - config.edgeMargin
  let minY := 0.0 + config.edgeMargin
  let maxY := 1.0 - config.edgeMargin

  let minDistSq := config.minDistance * config.minDistance
  let mut attempts := 0
  let maxAttempts := config.numPoints * 100

  while points.size < config.numPoints && attempts < maxAttempts do
    let (rng', x) := lcgNextInRange rng minX maxX
    rng := rng'
    let (rng', y) := lcgNextInRange rng minY maxY
    rng := rng'

    let candidate := Vec2.mk x y

    let mut tooClose := false
    for p in points do
      if candidate.distanceSquared p < minDistSq then
        tooClose := true
        break

    if !tooClose then
      points := points.push candidate

    attempts := attempts + 1

  return points

private def sumArea (polys : Array Polygon2D) : Float :=
  polys.foldl (fun acc poly => acc + poly.area) 0.0

testSuite "Voronoi"

test "unit square corners produce quadrant cells" := do
  let points := #[
    Vec2.mk 0.0 0.0,
    Vec2.mk 1.0 0.0,
    Vec2.mk 1.0 1.0,
    Vec2.mk 0.0 1.0
  ]
  let expectedCentroids := #[
    Vec2.mk 0.25 0.25,
    Vec2.mk 0.75 0.25,
    Vec2.mk 0.75 0.75,
    Vec2.mk 0.25 0.75
  ]
  let bounds := AABB2D.fromMinMax Vec2.zero Vec2.one
  match Voronoi.generate points bounds with
  | none => ensure false "expected voronoi polygons"
  | some polys =>
    ensure (polys.size == 4) "expected 4 polygons"
    for i in [:polys.size] do
      let poly := polys[i]!
      ensure (floatNear poly.area 0.25 1e-6) "quadrant area should be 0.25"
      let c := poly.centroid
      let expected := expectedCentroids[i]!
      ensure (floatNear c.x expected.x 1e-6) "centroid.x should match quadrant"
      ensure (floatNear c.y expected.y 1e-6) "centroid.y should match quadrant"
      ensure (poly.containsPointInclusive points[i]! 1e-6) "site should be inside cell"

test "lloyd relaxation moves corners to quadrant centroids" := do
  let points := #[
    Vec2.mk 0.0 0.0,
    Vec2.mk 1.0 0.0,
    Vec2.mk 1.0 1.0,
    Vec2.mk 0.0 1.0
  ]
  let expected := #[
    Vec2.mk 0.25 0.25,
    Vec2.mk 0.75 0.25,
    Vec2.mk 0.75 0.75,
    Vec2.mk 0.25 0.75
  ]
  let bounds := AABB2D.fromMinMax Vec2.zero Vec2.one
  match Voronoi.lloydRelaxation points bounds 1 with
  | none => ensure false "expected relaxed points"
  | some relaxed =>
    ensure (relaxed.size == points.size) "expected same number of points"
    for i in [:relaxed.size] do
      let p := relaxed[i]!
      let e := expected[i]!
      ensure (floatNear p.x e.x 1e-6) "relaxed.x should match quadrant centroid"
      ensure (floatNear p.y e.y 1e-6) "relaxed.y should match quadrant centroid"
      ensure (bounds.containsPoint p) "relaxed point should be inside bounds"

test "lloyd relaxation with zero iterations returns input" := do
  let points := #[
    Vec2.mk 0.2 0.1,
    Vec2.mk 0.9 0.3,
    Vec2.mk 0.4 0.8,
    Vec2.mk 0.7 0.7
  ]
  let bounds := AABB2D.fromMinMax Vec2.zero Vec2.one
  match Voronoi.lloydRelaxation points bounds 0 with
  | none => ensure false "expected relaxed points"
  | some relaxed =>
    ensure (relaxed.size == points.size) "expected same number of points"
    for i in [:points.size] do
      let a := points[i]!
      let b := relaxed[i]!
      ensure (floatNear a.x b.x 1e-12) "point.x should match input"
      ensure (floatNear a.y b.y 1e-12) "point.y should match input"

test "voronoi cells are convex for deterministic generator" := do
  let points := generateSeedPoints {
    numPoints := 24
    seedVal := 42
    minDistance := 0.12
    edgeMargin := 0.03
  }
  let bounds := AABB2D.fromMinMax Vec2.zero Vec2.one
  match Voronoi.generate points bounds with
  | none => ensure false "expected voronoi polygons"
  | some polys =>
    ensure (polys.size == points.size) "expected one polygon per site"
    for poly in polys do
      ensure (poly.vertices.size >= 3) "polygon should have >= 3 vertices"
      ensure poly.isConvex "polygon should be convex"

test "sites are inside cells and cover bounds" := do
  let points := generateSeedPoints {
    numPoints := 24
    seedVal := 7
    minDistance := 0.12
    edgeMargin := 0.03
  }
  let bounds := AABB2D.fromMinMax Vec2.zero Vec2.one
  match Voronoi.generate points bounds with
  | none => ensure false "expected voronoi polygons"
  | some polys =>
    ensure (polys.size == points.size) "expected one polygon per site"
    let area := sumArea polys
    ensure (Float.abs (area - 1.0) <= 1e-3) "cells should cover bounds area"
    for i in [:points.size] do
      let poly := polys[i]!
      ensure (poly.containsPointInclusive points[i]! 1e-6) "site should be inside cell"

test "collinear points return none" := do
  let points := #[
    Vec2.mk 0.0 0.0,
    Vec2.mk 0.5 0.0,
    Vec2.mk 1.0 0.0
  ]
  let bounds := AABB2D.fromMinMax Vec2.zero Vec2.one
  match Voronoi.generate points bounds with
  | none => pure ()
  | some _ => ensure false "expected none for collinear points"



end LinalgTests.VoronoiTests
