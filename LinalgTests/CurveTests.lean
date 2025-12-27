/-
  Tests for curves and splines.
-/

import Linalg
import Crucible

namespace LinalgTests.CurveTests

open Crucible
open Linalg

testSuite "Bezier2 (Quadratic)"

test "bezier2 at t=0 returns start point" := do
  let b := Bezier2.mk (Vec2.mk 0.0 0.0) (Vec2.mk 1.0 2.0) (Vec2.mk 2.0 0.0)
  let p := b.evalVec2 0.0
  ensure (floatNear p.x 0.0 0.0001) "x should be 0"
  ensure (floatNear p.y 0.0 0.0001) "y should be 0"

test "bezier2 at t=1 returns end point" := do
  let b := Bezier2.mk (Vec2.mk 0.0 0.0) (Vec2.mk 1.0 2.0) (Vec2.mk 2.0 0.0)
  let p := b.evalVec2 1.0
  ensure (floatNear p.x 2.0 0.0001) "x should be 2"
  ensure (floatNear p.y 0.0 0.0001) "y should be 0"

test "bezier2 at t=0.5 is influenced by control point" := do
  let b := Bezier2.mk (Vec2.mk 0.0 0.0) (Vec2.mk 1.0 2.0) (Vec2.mk 2.0 0.0)
  let p := b.evalVec2 0.5
  ensure (floatNear p.x 1.0 0.0001) "x should be 1"
  ensure (p.y > 0.5) "y should be pulled up by control point"

test "bezier2 tangent at t=0 points toward control" := do
  let b := Bezier2.mk Vec2.zero (Vec2.mk 1.0 1.0) (Vec2.mk 2.0 0.0)
  let d := b.derivativeVec2 0.0
  -- Derivative at t=0 is 2*(p1 - p0) = 2*(1,1) = (2,2)
  ensure (floatNear d.x 2.0 0.0001) "dx should be 2"
  ensure (floatNear d.y 2.0 0.0001) "dy should be 2"

test "bezier2 split produces valid subcurves" := do
  let b := Bezier2.mk Vec2.zero (Vec2.mk 1.0 2.0) (Vec2.mk 2.0 0.0)
  let (left, right) := b.splitVec2 0.5
  -- Left curve should start at p0 and end at midpoint
  let leftStart := left.evalVec2 0.0
  let leftEnd := left.evalVec2 1.0
  let mid := b.evalVec2 0.5
  ensure (floatNear leftStart.x 0.0 0.0001) "left start x"
  ensure (floatNear leftEnd.x mid.x 0.0001) "left end should be at midpoint"
  -- Right curve should start at midpoint and end at p2
  let rightStart := right.evalVec2 0.0
  let rightEnd := right.evalVec2 1.0
  ensure (floatNear rightStart.x mid.x 0.0001) "right start should be at midpoint"
  ensure (floatNear rightEnd.x 2.0 0.0001) "right end x"

testSuite "Bezier3 (Cubic)"

test "bezier3 at t=0 returns start point" := do
  let b := Bezier3.mk Vec2.zero (Vec2.mk 0.0 1.0) (Vec2.mk 1.0 1.0) (Vec2.mk 1.0 0.0)
  let p := b.evalVec2 0.0
  ensure (floatNear p.x 0.0 0.0001) "x should be 0"
  ensure (floatNear p.y 0.0 0.0001) "y should be 0"

test "bezier3 at t=1 returns end point" := do
  let b := Bezier3.mk Vec2.zero (Vec2.mk 0.0 1.0) (Vec2.mk 1.0 1.0) (Vec2.mk 1.0 0.0)
  let p := b.evalVec2 1.0
  ensure (floatNear p.x 1.0 0.0001) "x should be 1"
  ensure (floatNear p.y 0.0 0.0001) "y should be 0"

test "bezier3 tangent at t=0" := do
  let b := Bezier3.mk Vec2.zero (Vec2.mk 1.0 0.0) (Vec2.mk 2.0 0.0) (Vec2.mk 3.0 0.0)
  let d := b.derivativeVec2 0.0
  -- Derivative at t=0 is 3*(p1 - p0) = 3*(1,0) = (3,0)
  ensure (floatNear d.x 3.0 0.0001) "dx should be 3"
  ensure (floatNear d.y 0.0 0.0001) "dy should be 0"

test "bezier3 curvature of straight line is zero" := do
  let b := Bezier3.mk Vec2.zero (Vec2.mk 1.0 0.0) (Vec2.mk 2.0 0.0) (Vec2.mk 3.0 0.0)
  let k := b.curvatureVec2 0.5
  ensure (floatNear k 0.0 0.0001) "curvature of line should be 0"

test "bezier3 split preserves curve" := do
  let b := Bezier3.mk Vec2.zero (Vec2.mk 0.0 1.0) (Vec2.mk 1.0 1.0) (Vec2.mk 1.0 0.0)
  let (left, right) := b.splitVec2 0.5
  let mid := b.evalVec2 0.5
  let leftEnd := left.evalVec2 1.0
  let rightStart := right.evalVec2 0.0
  ensure (floatNear leftEnd.x mid.x 0.0001) "left end matches mid"
  ensure (floatNear rightStart.x mid.x 0.0001) "right start matches mid"

testSuite "Bezier3 Vec3"

test "bezier3 vec3 evaluation" := do
  let b := Bezier3.mk Vec3.zero (Vec3.mk 0.0 1.0 0.0) (Vec3.mk 1.0 1.0 0.0) (Vec3.mk 1.0 0.0 0.0)
  let start := b.evalVec3 0.0
  let endPt := b.evalVec3 1.0
  ensure (floatNear start.x 0.0 0.0001) "start x"
  ensure (floatNear endPt.x 1.0 0.0001) "end x"

testSuite "CatmullRom"

test "catmullrom passes through control points" := do
  let s : CatmullRom Vec2 := { p0 := Vec2.mk (-1.0) 0.0, p1 := Vec2.mk 0.0 0.0, p2 := Vec2.mk 1.0 0.0, p3 := Vec2.mk 2.0 0.0 }
  let start := s.evalVec2 0.0
  let endPt := s.evalVec2 1.0
  ensure (floatNear start.x 0.0 0.0001) "should pass through p1"
  ensure (floatNear endPt.x 1.0 0.0001) "should pass through p2"

test "catmullrom is smooth at endpoints" := do
  let s : CatmullRom Vec2 := { p0 := Vec2.mk (-1.0) 0.0, p1 := Vec2.mk 0.0 0.0, p2 := Vec2.mk 1.0 0.0, p3 := Vec2.mk 2.0 0.0 }
  let t0 := s.tangentVec2 0.0
  let t1 := s.tangentVec2 1.0
  -- For collinear points, tangent should be along x-axis
  ensure (floatNear t0.x 1.0 0.01) "tangent at 0 should point along x"
  ensure (floatNear t1.x 1.0 0.01) "tangent at 1 should point along x"

testSuite "SplinePath2"

test "spline path with 3 points" := do
  let path := SplinePath2.mk #[Vec2.zero, Vec2.mk 1.0 1.0, Vec2.mk 2.0 0.0] false
  let start := path.eval 0.0
  let endPt := path.eval 1.0
  ensure (floatNear start.x 0.0 0.01) "should start near first point"
  ensure (floatNear endPt.x 2.0 0.01) "should end near last point"

test "closed spline path" := do
  let path := SplinePath2.mk #[Vec2.zero, Vec2.mk 1.0 0.0, Vec2.mk 1.0 1.0, Vec2.mk 0.0 1.0] true
  let start := path.eval 0.0
  let end_ := path.eval 1.0
  -- For a closed path, start and end should be close (wraps around)
  ensure (start.distance end_ < 0.1) "closed path should wrap"

testSuite "Arc Length"

test "arc length table for straight line" := do
  -- A straight line from (0,0) to (1,0)
  let eval := fun t => Vec2.mk t 0.0
  let table := ArcLengthTable.build eval 100
  ensure (floatNear table.totalLength 1.0 0.01) "length should be 1"

test "arc length uToT gives uniform spacing" := do
  let eval := fun t => Vec2.mk t 0.0
  let table := ArcLengthTable.build eval 100
  -- For a straight line, u should equal t
  let t := table.uToT 0.5
  ensure (floatNear t 0.5 0.01) "uToT(0.5) should be ~0.5 for straight line"

test "arc length of curved path is longer than chord" := do
  -- Curved path that goes up and back down
  let b := Bezier3.mk Vec2.zero (Vec2.mk 0.0 1.0) (Vec2.mk 1.0 1.0) (Vec2.mk 1.0 0.0)
  let length := b.arcLengthVec2 0.001
  let chord := Vec2.zero.distance (Vec2.mk 1.0 0.0)
  ensure (length > chord) "arc length should exceed chord length"

#generate_tests

end LinalgTests.CurveTests
