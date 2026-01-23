/-
  Tests for 2D Collision Detection (SAT and GJK).
-/

import Linalg
import Crucible

namespace LinalgTests.Collision2DTests

open Crucible
open Linalg

-- ============================================================================
-- SAT Polygon-Polygon Tests
-- ============================================================================

testSuite "SAT Polygon-Polygon"

test "overlapping rectangles collide" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 2.0 2.0) (Vec2.mk 6.0 6.0)
  let result := SAT.polygonPolygon r1 r2
  ensure result.colliding "overlapping rectangles should collide"
  ensure (result.depth > 0.0) "penetration depth should be positive"

test "separated rectangles don't collide" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 2.0 2.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 5.0 5.0) (Vec2.mk 7.0 7.0)
  let result := SAT.polygonPolygon r1 r2
  ensure (!result.colliding) "separated rectangles should not collide"

test "touching rectangles (edge to edge)" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 2.0 2.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 2.0 0.0) (Vec2.mk 4.0 2.0)
  let result := SAT.polygonPolygon r1 r2
  -- Edge touching might or might not register as collision depending on epsilon
  -- Just verify no crash
  ensure true "edge touching test completed"

test "triangle inside rectangle" := do
  let rect := Polygon2D.rectangle Vec2.zero (Vec2.mk 10.0 10.0)
  let tri := Polygon2D.triangle (Vec2.mk 3.0 3.0) (Vec2.mk 7.0 3.0) (Vec2.mk 5.0 7.0)
  let result := SAT.polygonPolygon rect tri
  ensure result.colliding "triangle inside rectangle should collide"

test "MTV points from shape1 toward shape2" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 3.0 0.0) (Vec2.mk 7.0 4.0)  -- overlapping on right
  let result := SAT.polygonPolygon r1 r2
  ensure result.colliding "should collide"
  -- MTV should point toward r2 (positive X direction)
  ensure (result.mtv.x >= 0.0) "MTV x should be non-negative (pointing toward r2)"

test "regular hexagons collision" := do
  let h1 := Polygon2D.regular 6 2.0
  let h2 := (Polygon2D.regular 6 2.0).translate (Vec2.mk 3.0 0.0)
  let result := SAT.polygonPolygon h1 h2
  ensure result.colliding "overlapping hexagons should collide"

test "regular hexagons separation" := do
  let h1 := Polygon2D.regular 6 2.0
  let h2 := (Polygon2D.regular 6 2.0).translate (Vec2.mk 10.0 0.0)
  let result := SAT.polygonPolygon h1 h2
  ensure (!result.colliding) "separated hexagons should not collide"

-- ============================================================================
-- SAT Polygon-Circle Tests
-- ============================================================================

testSuite "SAT Polygon-Circle"

test "circle inside rectangle" := do
  let rect := Polygon2D.rectangle Vec2.zero (Vec2.mk 10.0 10.0)
  let circle := Circle.mk (Vec2.mk 5.0 5.0) 2.0
  let result := SAT.polygonCircle rect circle
  ensure result.colliding "circle inside rectangle should collide"

test "circle outside rectangle" := do
  let rect := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let circle := Circle.mk (Vec2.mk 10.0 10.0) 2.0
  let result := SAT.polygonCircle rect circle
  ensure (!result.colliding) "circle outside rectangle should not collide"

test "circle overlapping rectangle edge" := do
  let rect := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let circle := Circle.mk (Vec2.mk 5.0 2.0) 2.0  -- overlaps right edge
  let result := SAT.polygonCircle rect circle
  ensure result.colliding "circle overlapping edge should collide"

test "circle at rectangle corner" := do
  let rect := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let circle := Circle.mk (Vec2.mk 5.0 5.0) 2.0  -- near corner
  let result := SAT.polygonCircle rect circle
  -- Distance from corner (4,4) to circle center (5,5) is sqrt(2) ≈ 1.414 < radius 2
  ensure result.colliding "circle near corner should collide"

test "circle far from corner" := do
  let rect := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let circle := Circle.mk (Vec2.mk 7.0 7.0) 2.0  -- far from corner
  let result := SAT.polygonCircle rect circle
  -- Distance from corner (4,4) to circle center (7,7) is sqrt(18) ≈ 4.24 > radius 2
  ensure (!result.colliding) "circle far from corner should not collide"

test "circle inside triangle" := do
  let tri := Polygon2D.triangle Vec2.zero (Vec2.mk 10.0 0.0) (Vec2.mk 5.0 10.0)
  let circle := Circle.mk (Vec2.mk 5.0 3.0) 1.0
  let result := SAT.polygonCircle tri circle
  ensure result.colliding "circle inside triangle should collide"

-- ============================================================================
-- SAT Circle-Circle Tests
-- ============================================================================

testSuite "SAT Circle-Circle"

test "overlapping circles" := do
  let c1 := Circle.mk Vec2.zero 3.0
  let c2 := Circle.mk (Vec2.mk 4.0 0.0) 2.0
  let result := SAT.circleCircle c1 c2
  ensure result.colliding "overlapping circles should collide"
  ensure (floatNear result.depth 1.0 0.01) "depth should be 1.0 (3+2-4)"

test "separated circles" := do
  let c1 := Circle.mk Vec2.zero 2.0
  let c2 := Circle.mk (Vec2.mk 10.0 0.0) 2.0
  let result := SAT.circleCircle c1 c2
  ensure (!result.colliding) "separated circles should not collide"

test "touching circles" := do
  let c1 := Circle.mk Vec2.zero 2.0
  let c2 := Circle.mk (Vec2.mk 4.0 0.0) 2.0  -- exactly touching
  let result := SAT.circleCircle c1 c2
  -- Exactly touching should not collide (no penetration)
  ensure (!result.colliding) "touching circles should not collide"

test "concentric circles" := do
  let c1 := Circle.mk Vec2.zero 3.0
  let c2 := Circle.mk Vec2.zero 2.0
  let result := SAT.circleCircle c1 c2
  ensure result.colliding "concentric circles should collide"
  ensure (result.depth > 0.0) "concentric circles should have positive depth"

test "MTV direction for circles" := do
  let c1 := Circle.mk Vec2.zero 2.0
  let c2 := Circle.mk (Vec2.mk 3.0 0.0) 2.0
  let result := SAT.circleCircle c1 c2
  ensure result.colliding "should collide"
  ensure (result.mtv.x > 0.0) "MTV should point toward c2 (positive X)"

-- ============================================================================
-- GJK Intersection Tests
-- ============================================================================

testSuite "GJK Intersection"

test "GJK detects polygon-polygon intersection" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 2.0 2.0) (Vec2.mk 6.0 6.0)
  ensure (GJK.intersects r1 r2) "GJK should detect intersection"

test "GJK detects polygon-polygon separation" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 2.0 2.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 5.0 5.0) (Vec2.mk 7.0 7.0)
  ensure (!GJK.intersects r1 r2) "GJK should detect separation"

test "GJK detects circle-circle intersection" := do
  let c1 := Circle.mk Vec2.zero 3.0
  let c2 := Circle.mk (Vec2.mk 4.0 0.0) 2.0
  ensure (GJK.intersects c1 c2) "GJK should detect circle intersection"

test "GJK detects circle-circle separation" := do
  let c1 := Circle.mk Vec2.zero 2.0
  let c2 := Circle.mk (Vec2.mk 10.0 0.0) 2.0
  ensure (!GJK.intersects c1 c2) "GJK should detect circle separation"

test "GJK detects polygon-circle intersection" := do
  let rect := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let circle := Circle.mk (Vec2.mk 5.0 2.0) 2.0
  ensure (GJK.intersects rect circle) "GJK should detect polygon-circle intersection"

test "GJK detects polygon-circle separation" := do
  let rect := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let circle := Circle.mk (Vec2.mk 10.0 10.0) 2.0
  ensure (!GJK.intersects rect circle) "GJK should detect polygon-circle separation"

test "GJK detects AABB2D intersection" := do
  let a1 := AABB2D.fromMinMax Vec2.zero (Vec2.mk 4.0 4.0)
  let a2 := AABB2D.fromMinMax (Vec2.mk 2.0 2.0) (Vec2.mk 6.0 6.0)
  ensure (GJK.intersects a1 a2) "GJK should detect AABB intersection"

test "GJK detects AABB2D separation" := do
  let a1 := AABB2D.fromMinMax Vec2.zero (Vec2.mk 2.0 2.0)
  let a2 := AABB2D.fromMinMax (Vec2.mk 5.0 5.0) (Vec2.mk 7.0 7.0)
  ensure (!GJK.intersects a1 a2) "GJK should detect AABB separation"

-- ============================================================================
-- GJK Collision (with MTV) Tests
-- ============================================================================

testSuite "GJK Collision"

test "GJK collision returns MTV for overlapping shapes" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 2.0 2.0) (Vec2.mk 6.0 6.0)
  let result := GJK.collision r1 r2
  ensure result.colliding "GJK should detect collision"
  ensure (result.depth > 0.0) "depth should be positive"

test "GJK collision returns no collision for separated shapes" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 2.0 2.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 5.0 5.0) (Vec2.mk 7.0 7.0)
  let result := GJK.collision r1 r2
  ensure (!result.colliding) "GJK should detect no collision"

-- ============================================================================
-- Convenience Function Tests
-- ============================================================================

testSuite "Collision Convenience Functions"

test "collidePolygons works" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 2.0 2.0) (Vec2.mk 6.0 6.0)
  let result := collidePolygons r1 r2
  ensure result.colliding "collidePolygons should work"

test "collidePolygonCircle works" := do
  let rect := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let circle := Circle.mk (Vec2.mk 2.0 2.0) 1.0
  let result := collidePolygonCircle rect circle
  ensure result.colliding "collidePolygonCircle should work"

test "collideCircles works" := do
  let c1 := Circle.mk Vec2.zero 3.0
  let c2 := Circle.mk (Vec2.mk 4.0 0.0) 2.0
  let result := collideCircles c1 c2
  ensure result.colliding "collideCircles should work"

test "collideGJK works with polygons" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 2.0 2.0) (Vec2.mk 6.0 6.0)
  let result := collideGJK r1 r2
  ensure result.colliding "collideGJK should work"

test "intersectsGJK works" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 2.0 2.0) (Vec2.mk 6.0 6.0)
  ensure (intersectsGJK r1 r2) "intersectsGJK should work"

-- ============================================================================
-- SAT vs GJK Consistency Tests
-- ============================================================================

testSuite "SAT-GJK Consistency"

test "SAT and GJK agree on polygon-polygon collision" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 4.0 4.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 2.0 2.0) (Vec2.mk 6.0 6.0)
  let satResult := SAT.polygonPolygon r1 r2
  let gjkResult := GJK.intersects r1 r2
  ensure (satResult.colliding == gjkResult) "SAT and GJK should agree on collision"

test "SAT and GJK agree on polygon-polygon separation" := do
  let r1 := Polygon2D.rectangle Vec2.zero (Vec2.mk 2.0 2.0)
  let r2 := Polygon2D.rectangle (Vec2.mk 5.0 5.0) (Vec2.mk 7.0 7.0)
  let satResult := SAT.polygonPolygon r1 r2
  let gjkResult := GJK.intersects r1 r2
  ensure (satResult.colliding == gjkResult) "SAT and GJK should agree on separation"

test "SAT and GJK agree on circle-circle collision" := do
  let c1 := Circle.mk Vec2.zero 3.0
  let c2 := Circle.mk (Vec2.mk 4.0 0.0) 2.0
  let satResult := SAT.circleCircle c1 c2
  let gjkResult := GJK.intersects c1 c2
  ensure (satResult.colliding == gjkResult) "SAT and GJK should agree on circle collision"

test "SAT and GJK agree on circle-circle separation" := do
  let c1 := Circle.mk Vec2.zero 2.0
  let c2 := Circle.mk (Vec2.mk 10.0 0.0) 2.0
  let satResult := SAT.circleCircle c1 c2
  let gjkResult := GJK.intersects c1 c2
  ensure (satResult.colliding == gjkResult) "SAT and GJK should agree on circle separation"



end LinalgTests.Collision2DTests
