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

testSuite "Voronoi"

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

#generate_tests

end LinalgTests.VoronoiTests
