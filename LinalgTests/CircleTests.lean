/-
  Tests for Circle operations.
-/

import Linalg
import Crucible

namespace LinalgTests.CircleTests

open Crucible
open Linalg

testSuite "Circle Basic"

test "unit circle has radius 1" := do
  let c := Circle.unit
  ensure (floatNear c.radius 1.0 0.0001) "radius should be 1"
  ensure (floatNear c.center.x 0.0 0.0001) "center.x should be 0"
  ensure (floatNear c.center.y 0.0 0.0001) "center.y should be 0"

test "diameter is twice radius" := do
  let c := Circle.atOrigin 5.0
  ensure (floatNear c.diameter 10.0 0.0001) "diameter should be 10"

test "circumference formula" := do
  let c := Circle.atOrigin 1.0
  ensure (floatNear c.circumference (2.0 * Float.pi) 0.0001) "circumference should be 2*pi"

test "area formula" := do
  let c := Circle.atOrigin 2.0
  ensure (floatNear c.area (4.0 * Float.pi) 0.0001) "area should be 4*pi"

testSuite "Circle Containment"

test "circle contains its center" := do
  let c := Circle.mk' 5.0 5.0 3.0
  ensure (c.containsPoint (Vec2.mk 5.0 5.0)) "should contain center"

test "circle contains point inside" := do
  let c := Circle.atOrigin 5.0
  ensure (c.containsPoint (Vec2.mk 1.0 1.0)) "should contain (1,1)"

test "circle does not contain point outside" := do
  let c := Circle.atOrigin 1.0
  ensure (!c.containsPoint (Vec2.mk 2.0 0.0)) "should not contain (2,0)"

test "containsPointInclusive includes boundary" := do
  let c := Circle.atOrigin 1.0
  ensure (c.containsPointInclusive (Vec2.mk 1.0 0.0)) "should include boundary point"

test "signed distance inside is negative" := do
  let c := Circle.atOrigin 2.0
  let d := c.signedDistance (Vec2.mk 0.5 0.0)
  ensure (d < 0.0) "signed distance should be negative inside"

test "signed distance outside is positive" := do
  let c := Circle.atOrigin 1.0
  let d := c.signedDistance (Vec2.mk 3.0 0.0)
  ensure (d > 0.0) "signed distance should be positive outside"

testSuite "Circle Geometric Queries"

test "closest point on boundary" := do
  let c := Circle.atOrigin 2.0
  let closest := c.closestPoint (Vec2.mk 5.0 0.0)
  ensure (floatNear closest.x 2.0 0.0001) "x should be 2"
  ensure (floatNear closest.y 0.0 0.0001) "y should be 0"

test "point at angle 0 is on positive x" := do
  let c := Circle.atOrigin 3.0
  let p := c.pointAtAngle 0.0
  ensure (floatNear p.x 3.0 0.0001) "x should be 3"
  ensure (floatNear p.y 0.0 0.0001) "y should be 0"

test "point at angle pi/2 is on positive y" := do
  let c := Circle.atOrigin 2.0
  let p := c.pointAtAngle (Float.pi / 2.0)
  ensure (floatNear p.x 0.0 0.0001) "x should be 0"
  ensure (floatNear p.y 2.0 0.0001) "y should be 2"

testSuite "Circle-Circle Intersection"

test "overlapping circles intersect" := do
  let c1 := Circle.atOrigin 2.0
  let c2 := Circle.mk' 3.0 0.0 2.0
  ensure (c1.intersects c2) "circles should intersect"

test "distant circles do not intersect" := do
  let c1 := Circle.atOrigin 1.0
  let c2 := Circle.mk' 10.0 0.0 1.0
  ensure (!c1.intersects c2) "circles should not intersect"

test "touching circles intersect" := do
  let c1 := Circle.atOrigin 1.0
  let c2 := Circle.mk' 2.0 0.0 1.0
  ensure (c1.intersects c2) "touching circles should intersect"

test "tangent circles have one intersection point" := do
  let c1 := Circle.atOrigin 1.0
  let c2 := Circle.mk' 2.0 0.0 1.0
  let points := c1.intersectionPoints c2
  ensure (points.size == 1) "should have 1 intersection point"

test "overlapping circles have two intersection points" := do
  let c1 := Circle.atOrigin 2.0
  let c2 := Circle.mk' 2.0 0.0 2.0
  let points := c1.intersectionPoints c2
  ensure (points.size == 2) "should have 2 intersection points"

test "separate circles have no intersection points" := do
  let c1 := Circle.atOrigin 1.0
  let c2 := Circle.mk' 10.0 0.0 1.0
  let points := c1.intersectionPoints c2
  ensure (points.size == 0) "should have 0 intersection points"

test "outer circle contains inner" := do
  let outer := Circle.atOrigin 5.0
  let inner := Circle.mk' 1.0 0.0 1.0
  ensure (outer.containsCircle inner) "outer should contain inner"

testSuite "Circle From Three Points"

test "fromThreePoints creates valid circle" := do
  let p1 := Vec2.mk 0.0 1.0
  let p2 := Vec2.mk 1.0 0.0
  let p3 := Vec2.mk (-1.0) 0.0
  match Circle.fromThreePoints p1 p2 p3 with
  | some c =>
    ensure (floatNear c.center.y 0.0 0.0001) "center.y should be 0"
    ensure (floatNear c.radius 1.0 0.0001) "radius should be 1"
  | none => ensure false "should create valid circle"

test "collinear points return none" := do
  let p1 := Vec2.mk 0.0 0.0
  let p2 := Vec2.mk 1.0 0.0
  let p3 := Vec2.mk 2.0 0.0
  match Circle.fromThreePoints p1 p2 p3 with
  | some _ => ensure false "should return none for collinear points"
  | none => ensure true "correctly returned none"

testSuite "Circle Transformations"

test "translate moves center" := do
  let c := Circle.atOrigin 1.0
  let translated := c.translate (Vec2.mk 5.0 3.0)
  ensure (floatNear translated.center.x 5.0 0.0001) "x should be 5"
  ensure (floatNear translated.center.y 3.0 0.0001) "y should be 3"
  ensure (floatNear translated.radius 1.0 0.0001) "radius unchanged"

test "scale changes radius" := do
  let c := Circle.atOrigin 2.0
  let scaled := c.scale 3.0
  ensure (floatNear scaled.radius 6.0 0.0001) "radius should be 6"

test "bounding box is correct" := do
  let c := Circle.mk' 5.0 5.0 2.0
  let (minV, maxV) := c.boundingBox
  ensure (floatNear minV.x 3.0 0.0001) "min.x should be 3"
  ensure (floatNear minV.y 3.0 0.0001) "min.y should be 3"
  ensure (floatNear maxV.x 7.0 0.0001) "max.x should be 7"
  ensure (floatNear maxV.y 7.0 0.0001) "max.y should be 7"



end LinalgTests.CircleTests
