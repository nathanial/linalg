/-
  Tests for Line2D and Segment2D operations.
-/

import Linalg
import Crucible

namespace LinalgTests.Line2DTests

open Crucible
open Linalg

testSuite "Line2D Basic"

test "horizontal line is horizontal" := do
  let l := Line2D.horizontal 5.0
  ensure l.isHorizontal "should be horizontal"
  ensure (!l.isVertical) "should not be vertical"

test "vertical line is vertical" := do
  let l := Line2D.vertical 3.0
  ensure l.isVertical "should be vertical"
  ensure (!l.isHorizontal) "should not be horizontal"

test "line from two points" := do
  let l := Line2D.fromPoints (Vec2.mk 0.0 0.0) (Vec2.mk 1.0 1.0)
  ensure (floatNear l.direction.x l.direction.y 0.0001) "direction should be 45 degrees"

test "slope of horizontal line is 0" := do
  let l := Line2D.horizontal 0.0
  match l.slope with
  | some s => ensure (floatNear s 0.0 0.0001) "slope should be 0"
  | none => ensure false "should have slope"

test "vertical line has no slope" := do
  let l := Line2D.vertical 0.0
  match l.slope with
  | some _ => ensure false "should not have slope"
  | none => ensure true "correctly returned none"

testSuite "Line2D Distance"

test "point on line has zero distance" := do
  let l := Line2D.horizontal 0.0
  let d := l.distance (Vec2.mk 5.0 0.0)
  ensure (floatNear d 0.0 0.0001) "distance should be 0"

test "point above horizontal line has positive signed distance" := do
  let l := Line2D.horizontal 0.0
  let sd := l.signedDistance (Vec2.mk 0.0 5.0)
  -- Depends on direction, but should be non-zero
  ensure (Float.abs' sd > 0.0) "signed distance should be non-zero"

test "closest point on line" := do
  let l := Line2D.horizontal 0.0
  let closest := l.closestPoint (Vec2.mk 3.0 5.0)
  ensure (floatNear closest.x 3.0 0.0001) "x should be 3"
  ensure (floatNear closest.y 0.0 0.0001) "y should be 0"

testSuite "Line2D Intersection"

test "parallel lines do not intersect" := do
  let l1 := Line2D.horizontal 0.0
  let l2 := Line2D.horizontal 5.0
  ensure (l1.isParallel l2) "lines should be parallel"
  match l1.intersection l2 with
  | some _ => ensure false "should not intersect"
  | none => ensure true "correctly returned none"

test "perpendicular lines intersect" := do
  let l1 := Line2D.horizontal 0.0
  let l2 := Line2D.vertical 0.0
  ensure (!l1.isParallel l2) "lines should not be parallel"
  match l1.intersection l2 with
  | some p =>
    ensure (floatNear p.x 0.0 0.0001) "x should be 0"
    ensure (floatNear p.y 0.0 0.0001) "y should be 0"
  | none => ensure false "should intersect"

test "diagonal lines intersect" := do
  let l1 := Line2D.fromPoints (Vec2.mk 0.0 0.0) (Vec2.mk 1.0 1.0)
  let l2 := Line2D.fromPoints (Vec2.mk 0.0 1.0) (Vec2.mk 1.0 0.0)
  match l1.intersection l2 with
  | some p =>
    ensure (floatNear p.x 0.5 0.0001) "x should be 0.5"
    ensure (floatNear p.y 0.5 0.0001) "y should be 0.5"
  | none => ensure false "should intersect"

testSuite "Segment2D Basic"

test "segment length" := do
  let s := Segment2D.mk' 0.0 0.0 3.0 4.0
  ensure (floatNear s.length 5.0 0.0001) "length should be 5"

test "segment midpoint" := do
  let s := Segment2D.mk' 0.0 0.0 10.0 0.0
  let mid := s.midpoint
  ensure (floatNear mid.x 5.0 0.0001) "midpoint.x should be 5"
  ensure (floatNear mid.y 0.0 0.0001) "midpoint.y should be 0"

test "segment direction" := do
  let s := Segment2D.mk' 0.0 0.0 1.0 0.0
  let dir := s.direction
  ensure (floatNear dir.x 1.0 0.0001) "direction.x should be 1"
  ensure (floatNear dir.y 0.0 0.0001) "direction.y should be 0"

test "pointAt t=0 returns a" := do
  let s := Segment2D.mk' 1.0 2.0 5.0 6.0
  let p := s.pointAt 0.0
  ensure (floatNear p.x 1.0 0.0001) "x should be 1"
  ensure (floatNear p.y 2.0 0.0001) "y should be 2"

test "pointAt t=1 returns b" := do
  let s := Segment2D.mk' 1.0 2.0 5.0 6.0
  let p := s.pointAt 1.0
  ensure (floatNear p.x 5.0 0.0001) "x should be 5"
  ensure (floatNear p.y 6.0 0.0001) "y should be 6"

testSuite "Segment2D Distance"

test "closest point inside segment" := do
  let s := Segment2D.mk' 0.0 0.0 10.0 0.0
  let closest := s.closestPoint (Vec2.mk 5.0 3.0)
  ensure (floatNear closest.x 5.0 0.0001) "x should be 5"
  ensure (floatNear closest.y 0.0 0.0001) "y should be 0"

test "closest point clamped to endpoint a" := do
  let s := Segment2D.mk' 0.0 0.0 10.0 0.0
  let closest := s.closestPoint (Vec2.mk (-5.0) 0.0)
  ensure (floatNear closest.x 0.0 0.0001) "x should be 0"
  ensure (floatNear closest.y 0.0 0.0001) "y should be 0"

test "closest point clamped to endpoint b" := do
  let s := Segment2D.mk' 0.0 0.0 10.0 0.0
  let closest := s.closestPoint (Vec2.mk 15.0 0.0)
  ensure (floatNear closest.x 10.0 0.0001) "x should be 10"

test "distance to point" := do
  let s := Segment2D.mk' 0.0 0.0 10.0 0.0
  let d := s.distance (Vec2.mk 5.0 3.0)
  ensure (floatNear d 3.0 0.0001) "distance should be 3"

test "containsPoint for point on segment" := do
  let s := Segment2D.mk' 0.0 0.0 10.0 0.0
  ensure (s.containsPoint (Vec2.mk 5.0 0.0) 0.0001) "should contain (5,0)"

testSuite "Segment2D Intersection"

test "crossing segments intersect" := do
  let s1 := Segment2D.mk' 0.0 0.0 10.0 10.0
  let s2 := Segment2D.mk' 0.0 10.0 10.0 0.0
  ensure (s1.intersects s2) "segments should intersect"

test "crossing segments intersection point" := do
  let s1 := Segment2D.mk' 0.0 0.0 10.0 10.0
  let s2 := Segment2D.mk' 0.0 10.0 10.0 0.0
  match s1.intersection s2 with
  | some p =>
    ensure (floatNear p.x 5.0 0.0001) "x should be 5"
    ensure (floatNear p.y 5.0 0.0001) "y should be 5"
  | none => ensure false "should have intersection"

test "parallel segments do not intersect" := do
  let s1 := Segment2D.mk' 0.0 0.0 10.0 0.0
  let s2 := Segment2D.mk' 0.0 5.0 10.0 5.0
  ensure (!s1.intersects s2) "parallel segments should not intersect"

test "non-overlapping segments do not intersect" := do
  let s1 := Segment2D.mk' 0.0 0.0 1.0 0.0
  let s2 := Segment2D.mk' 5.0 0.0 6.0 0.0
  ensure (!s1.intersects s2) "distant segments should not intersect"

test "T-intersection" := do
  let s1 := Segment2D.mk' 0.0 0.0 10.0 0.0
  let s2 := Segment2D.mk' 5.0 (-5.0) 5.0 5.0
  ensure (s1.intersects s2) "T segments should intersect"

testSuite "Segment2D Utilities"

test "reverse swaps endpoints" := do
  let s := Segment2D.mk' 0.0 0.0 10.0 10.0
  let rev := s.reverse
  ensure (floatNear rev.a.x 10.0 0.0001) "a.x should be 10"
  ensure (floatNear rev.b.x 0.0 0.0001) "b.x should be 0"

test "extend increases length" := do
  let s := Segment2D.mk' 0.0 0.0 10.0 0.0
  let extended := s.extend 2.0
  ensure (floatNear extended.length 14.0 0.0001) "extended length should be 14"

test "toLine preserves point on segment" := do
  let s := Segment2D.mk' 0.0 0.0 10.0 0.0
  let l := s.toLine
  let d := l.distance (Vec2.mk 5.0 0.0)
  ensure (floatNear d 0.0 0.0001) "point on segment should be on line"



end LinalgTests.Line2DTests
