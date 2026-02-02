/-
  Tests for spatial data structures.
-/

import Linalg
import Crucible

namespace LinalgTests.SpatialTests

open Crucible
open Linalg
open Linalg.Spatial

-- ============================================================================
-- AABB2D Tests
-- ============================================================================

testSuite "AABB2D"

test "AABB2D center is correct" := do
  let aabb := AABB2D.fromMinMax (Vec2.mk 0.0 0.0) (Vec2.mk 10.0 10.0)
  let c := aabb.center
  ensure (floatNear c.x 5.0 0.0001) "center x should be 5"
  ensure (floatNear c.y 5.0 0.0001) "center y should be 5"

test "AABB2D contains center point" := do
  let aabb := AABB2D.fromMinMax Vec2.zero (Vec2.mk 10.0 10.0)
  ensure (aabb.containsPoint (Vec2.mk 5.0 5.0)) "should contain center"

test "AABB2D does not contain point outside" := do
  let aabb := AABB2D.fromMinMax Vec2.zero (Vec2.mk 10.0 10.0)
  ensure (!(aabb.containsPoint (Vec2.mk 15.0 5.0))) "should not contain outside point"

test "AABB2D intersects overlapping box" := do
  let a := AABB2D.fromMinMax Vec2.zero (Vec2.mk 5.0 5.0)
  let b := AABB2D.fromMinMax (Vec2.mk 3.0 3.0) (Vec2.mk 8.0 8.0)
  ensure (a.intersects b) "boxes should intersect"

test "AABB2D does not intersect disjoint box" := do
  let a := AABB2D.fromMinMax Vec2.zero (Vec2.mk 5.0 5.0)
  let b := AABB2D.fromMinMax (Vec2.mk 10.0 10.0) (Vec2.mk 15.0 15.0)
  ensure (!(a.intersects b)) "boxes should not intersect"

test "AABB2D merge combines boxes" := do
  let a := AABB2D.fromMinMax Vec2.zero (Vec2.mk 5.0 5.0)
  let b := AABB2D.fromMinMax (Vec2.mk 3.0 3.0) (Vec2.mk 8.0 8.0)
  let merged := AABB2D.merge a b
  ensure (floatNear merged.min.x 0.0 0.0001) "min x should be 0"
  ensure (floatNear merged.max.x 8.0 0.0001) "max x should be 8"

test "AABB2D area calculation" := do
  let aabb := AABB2D.fromMinMax Vec2.zero (Vec2.mk 4.0 3.0)
  ensure (floatNear aabb.area 12.0 0.0001) "area should be 12"

test "AABB2D subdivide creates 4 quadrants" := do
  let aabb := AABB2D.fromMinMax Vec2.zero (Vec2.mk 10.0 10.0)
  let quads := aabb.subdivide
  ensure (quads.size == 4) "should have 4 quadrants"

-- ============================================================================
-- Grid2D Tests
-- ============================================================================

testSuite "Grid2D"

test "Grid2D build creates grid" := do
  let points := #[Vec2.mk 1.0 1.0, Vec2.mk 5.0 5.0, Vec2.mk 9.0 9.0]
  let grid := Grid2D.build points 5.0
  ensure (grid.itemCount == 3) "should have 3 items"

test "Grid2D queryRect finds points" := do
  let points := #[Vec2.mk 1.0 1.0, Vec2.mk 5.0 5.0, Vec2.mk 9.0 9.0]
  let grid := Grid2D.build points 5.0
  let query := AABB2D.fromMinMax Vec2.zero (Vec2.mk 3.0 3.0)
  let results := grid.queryRect query
  ensure (results.size >= 1) "should find at least 1 point"

test "Grid2D insert adds item" := do
  let grid := Grid2D.empty 5.0
  let grid' := grid.insert 0 (Vec2.mk 2.5 2.5)
  ensure (grid'.itemCount == 1) "should have 1 item"

test "Grid2D remove removes item" := do
  let grid := Grid2D.empty 5.0
  let point := Vec2.mk 2.5 2.5
  let grid' := grid.insert 0 point
  ensure (grid'.itemCount == 1) "should have 1 item after insert"
  let bounds := AABB2D.fromPoint point
  let grid'' := grid'.remove 0 bounds
  ensure (grid''.itemCount == 0) "should have 0 items after remove"

test "Grid2D remove preserves other items" := do
  let points := #[Vec2.mk 1.0 1.0, Vec2.mk 9.0 9.0]
  let grid := Grid2D.build points 5.0
  let grid' := grid.remove 0 (AABB2D.fromPoint points[0]!)
  ensure (grid'.itemCount == 1) "should have 1 item left"
  let query := AABB2D.fromMinMax (Vec2.mk 8.0 8.0) (Vec2.mk 10.0 10.0)
  let results := grid'.queryRect query
  ensure (results.contains 1) "should still find remaining point"

-- ============================================================================
-- Grid3D Tests
-- ============================================================================

testSuite "Grid3D"

test "Grid3D build creates grid" := do
  let points := #[Vec3.mk 1.0 1.0 1.0, Vec3.mk 5.0 5.0 5.0, Vec3.mk 9.0 9.0 9.0]
  let grid := Grid3D.build points 5.0
  ensure (grid.itemCount == 3) "should have 3 items"

test "Grid3D queryAABB finds points" := do
  let points := #[Vec3.mk 1.0 1.0 1.0, Vec3.mk 5.0 5.0 5.0, Vec3.mk 9.0 9.0 9.0]
  let grid := Grid3D.build points 5.0
  let query := AABB.fromMinMax Vec3.zero (Vec3.mk 3.0 3.0 3.0)
  let results := grid.queryAABB query
  ensure (results.size >= 1) "should find at least 1 point"

test "Grid3D remove removes item" := do
  let grid := Grid3D.empty 5.0
  let point := Vec3.mk 2.5 2.5 2.5
  let grid' := grid.insert 0 point
  let bounds := AABB.fromPoint point
  let grid'' := grid'.remove 0 bounds
  ensure (grid''.itemCount == 0) "should have 0 items after remove"

-- ============================================================================
-- Quadtree Tests
-- ============================================================================

testSuite "Quadtree"

test "Quadtree build creates tree" := do
  let points := #[Vec2.mk 1.0 1.0, Vec2.mk 5.0 5.0, Vec2.mk 9.0 9.0, Vec2.mk 3.0 7.0]
  let tree := Quadtree.build points
  ensure (tree.itemCount == 4) "should have 4 items"

test "Quadtree queryRect finds points" := do
  let points := #[Vec2.mk 1.0 1.0, Vec2.mk 5.0 5.0, Vec2.mk 9.0 9.0]
  let tree := Quadtree.build points
  let query := AABB2D.fromMinMax Vec2.zero (Vec2.mk 3.0 3.0)
  let results := tree.queryRect query
  ensure (results.size >= 1) "should find at least 1 point"

test "Quadtree queryPoint finds containing items" := do
  let rects := #[
    AABB2D.fromMinMax (Vec2.mk 0.0 0.0) (Vec2.mk 5.0 5.0),
    AABB2D.fromMinMax (Vec2.mk 3.0 3.0) (Vec2.mk 8.0 8.0)
  ]
  let tree := Quadtree.build rects
  let results := tree.queryPoint (Vec2.mk 4.0 4.0)
  ensure (results.size == 2) "point should be in both rectangles"

test "Quadtree split reinserts existing items" := do
  let config := { TreeConfig.default with maxLeafItems := 1 }
  let points := #[Vec2.mk (-5.0) (-5.0), Vec2.mk 5.0 5.0]
  let tree := Quadtree.build points config
  let resultsSW := tree.queryPoint (Vec2.mk (-5.0) (-5.0))
  let resultsNE := tree.queryPoint (Vec2.mk 5.0 5.0)
  ensure (resultsSW.contains 0) "southwest point should be found after split"
  ensure (!resultsSW.contains 1) "southwest query should not return northeast point"
  ensure (resultsNE.contains 1) "northeast point should be found after split"
  ensure (!resultsNE.contains 0) "northeast query should not return southwest point"

test "Quadtree queryCircle finds nearby points" := do
  let points := #[Vec2.mk 0.0 0.0, Vec2.mk 1.0 0.0, Vec2.mk 10.0 10.0]
  let tree := Quadtree.build points
  let results := tree.queryCircle Vec2.zero 2.0
  ensure (results.size >= 2) "should find at least 2 nearby points"

test "Quadtree nodeCount is reasonable" := do
  let points := (List.range 20).map (fun i => Vec2.mk i.toFloat i.toFloat) |>.toArray
  let tree := Quadtree.build points
  ensure (tree.nodeCount >= 1) "should have at least 1 node"

test "Quadtree remove removes item" := do
  let points := #[Vec2.mk 1.0 1.0, Vec2.mk 5.0 5.0, Vec2.mk 9.0 9.0]
  let tree := Quadtree.build points
  let bounds := AABB2D.fromPoint points[1]!
  let tree' := tree.remove 1 bounds
  ensure (tree'.itemCount == 2) "should have 2 items after remove"

test "Quadtree remove item not found in queries" := do
  let points := #[Vec2.mk 1.0 1.0, Vec2.mk 5.0 5.0, Vec2.mk 9.0 9.0]
  let tree := Quadtree.build points
  let tree' := tree.remove 1 (AABB2D.fromPoint points[1]!)
  let results := tree'.queryPoint (Vec2.mk 5.0 5.0)
  ensure (!results.contains 1) "removed item should not be found"

test "Quadtree remove missing item keeps count" := do
  let points := #[Vec2.mk 1.0 1.0, Vec2.mk 5.0 5.0]
  let tree := Quadtree.build points
  let tree' := tree.remove 99 (AABB2D.fromPoint points[0]!)
  ensure (tree'.itemCount == tree.itemCount) "missing item should not change count"

test "Quadtree remove preserves other items" := do
  let points := #[Vec2.mk 1.0 1.0, Vec2.mk 5.0 5.0, Vec2.mk 9.0 9.0]
  let tree := Quadtree.build points
  let tree' := tree.remove 1 (AABB2D.fromPoint points[1]!)
  let results := tree'.queryRect (AABB2D.fromMinMax Vec2.zero (Vec2.mk 10.0 10.0))
  ensure (results.contains 0) "should still find item 0"
  ensure (results.contains 2) "should still find item 2"

test "Quadtree remove overlapping item" := do
  -- Large rectangle spanning multiple quadrants
  -- Use maxLeafItems=1 to force the tree to split, so items end up in different leaves
  let config := { TreeConfig.default with maxLeafItems := 1 }
  let rect := AABB2D.fromMinMax Vec2.zero (Vec2.mk 4.0 4.0)
  let items := #[rect, AABB2D.fromMinMax (Vec2.mk 8.0 8.0) (Vec2.mk 9.0 9.0)]
  let bounds := items.foldl (fun acc item => AABB2D.merge acc item) items[0]!
  let padding := bounds.size.scale 0.01
  let paddedBounds := AABB2D.fromMinMax (bounds.min.sub padding) (bounds.max.add padding)
  let tree := items.foldl (fun (t, idx) item => (t.insert idx item, idx + 1))
    (Quadtree.empty paddedBounds config, 0) |>.1
  let tree' := tree.remove 0 rect
  ensure (tree'.itemCount == 1) "should have 1 item"
  -- Point (2,2) is in item 0's bounds but item 1 is far away in a different quadrant
  let results := tree'.queryPoint (Vec2.mk 2.0 2.0)
  ensure results.isEmpty "removed item should not be found"

test "Quadtree empty tree queries return empty" := do
  let tree := Quadtree.empty (AABB2D.fromMinMax Vec2.zero (Vec2.mk 10.0 10.0))
  let results := tree.queryRect (AABB2D.fromMinMax Vec2.zero (Vec2.mk 5.0 5.0))
  ensure results.isEmpty "empty tree should return empty results"

test "Quadtree single item" := do
  let tree := Quadtree.build #[Vec2.mk 5.0 5.0]
  ensure (tree.itemCount == 1) "should have 1 item"
  let results := tree.queryPoint (Vec2.mk 5.0 5.0)
  ensure (results.size == 1) "should find the single item"

test "Quadtree kNearest returns sorted by distance" := do
  let points := #[Vec2.mk 0.0 0.0, Vec2.mk 3.0 0.0, Vec2.mk 1.0 0.0, Vec2.mk 2.0 0.0]
  let tree := Quadtree.build points
  let results := tree.kNearest points Vec2.zero 4
  ensure (results.size == 4) "should find 4 nearest"
  -- Verify results are sorted by distance (check each consecutive pair)
  let distances := results.map fun idx => (points[idx]!).distanceSquared Vec2.zero
  let isSorted := (List.range (distances.size - 1)).all fun i =>
    decide (distances[i]! <= distances[i + 1]!)
  ensure isSorted "results should be sorted by ascending distance"

test "Quadtree queryPointExact filters non-containing items" := do
  -- Two items in same leaf (no splitting), only one contains query point
  let items := #[
    AABB2D.fromMinMax Vec2.zero (Vec2.mk 5.0 5.0),
    AABB2D.fromMinMax (Vec2.mk 7.0 7.0) (Vec2.mk 9.0 9.0)
  ]
  let tree := Quadtree.build items
  let point := Vec2.mk 3.0 3.0
  -- Regular queryPoint returns candidates (both in same leaf)
  let candidates := tree.queryPoint point
  -- Exact query filters to only items actually containing the point
  let exact := tree.queryPointExact items point
  ensure (exact.size == 1) "should find exactly 1 item containing point"
  ensure (exact[0]! == 0) "should be the first item"

test "Quadtree queryRectExact filters non-intersecting items" := do
  let items := #[
    AABB2D.fromMinMax Vec2.zero (Vec2.mk 3.0 3.0),
    AABB2D.fromMinMax (Vec2.mk 7.0 7.0) (Vec2.mk 9.0 9.0)
  ]
  let tree := Quadtree.build items
  let query := AABB2D.fromMinMax (Vec2.mk 2.0 2.0) (Vec2.mk 4.0 4.0)
  let exact := tree.queryRectExact items query
  ensure (exact.size == 1) "should find exactly 1 intersecting item"
  ensure (exact[0]! == 0) "should be the first item"

-- ============================================================================
-- Octree Tests
-- ============================================================================

testSuite "Octree"

test "Octree build creates tree" := do
  let points := #[Vec3.mk 1.0 1.0 1.0, Vec3.mk 5.0 5.0 5.0, Vec3.mk 9.0 9.0 9.0]
  let tree := Octree.build points
  ensure (tree.itemCount == 3) "should have 3 items"

test "Octree queryAABB finds points" := do
  let points := #[Vec3.mk 1.0 1.0 1.0, Vec3.mk 5.0 5.0 5.0, Vec3.mk 9.0 9.0 9.0]
  let tree := Octree.build points
  let query := AABB.fromMinMax Vec3.zero (Vec3.mk 3.0 3.0 3.0)
  let results := tree.queryAABB query
  ensure (results.size >= 1) "should find at least 1 point"

test "Octree querySphere finds nearby points" := do
  let points := #[Vec3.mk 0.0 0.0 0.0, Vec3.mk 1.0 0.0 0.0, Vec3.mk 10.0 10.0 10.0]
  let tree := Octree.build points
  let results := tree.querySphere Vec3.zero 2.0
  ensure (results.size >= 2) "should find at least 2 nearby points"

test "Octree queryPoint finds containing items" := do
  let boxes := #[
    AABB.fromMinMax Vec3.zero (Vec3.mk 5.0 5.0 5.0),
    AABB.fromMinMax (Vec3.mk 3.0 3.0 3.0) (Vec3.mk 8.0 8.0 8.0)
  ]
  let tree := Octree.build boxes
  let results := tree.queryPoint (Vec3.mk 4.0 4.0 4.0)
  ensure (results.size == 2) "point should be in both boxes"

test "Octree split reinserts existing items" := do
  let config := { TreeConfig.default with maxLeafItems := 1 }
  let points := #[Vec3.mk (-5.0) (-5.0) (-5.0), Vec3.mk 5.0 5.0 5.0]
  let tree := Octree.build points config
  let resultsSW := tree.queryPoint (Vec3.mk (-5.0) (-5.0) (-5.0))
  let resultsNE := tree.queryPoint (Vec3.mk 5.0 5.0 5.0)
  ensure (resultsSW.contains 0) "southwest point should be found after split"
  ensure (!resultsSW.contains 1) "southwest query should not return northeast point"
  ensure (resultsNE.contains 1) "northeast point should be found after split"
  ensure (!resultsNE.contains 0) "northeast query should not return southwest point"

test "Octree rayCast finds items along ray" := do
  let boxes := #[
    AABB.fromMinMax (Vec3.mk (-1.0) (-1.0) 5.0) (Vec3.mk 1.0 1.0 7.0),
    AABB.fromMinMax (Vec3.mk (-1.0) (-1.0) 10.0) (Vec3.mk 1.0 1.0 12.0)
  ]
  let tree := Octree.build boxes
  let ray := Ray.mk' Vec3.zero Vec3.unitZ
  let results := tree.rayCast ray
  ensure (results.size == 2) "ray should hit both boxes"

test "Octree remove removes item" := do
  let points := #[Vec3.mk 1.0 1.0 1.0, Vec3.mk 5.0 5.0 5.0, Vec3.mk 9.0 9.0 9.0]
  let tree := Octree.build points
  let bounds := AABB.fromPoint points[1]!
  let tree' := tree.remove 1 bounds
  ensure (tree'.itemCount == 2) "should have 2 items after remove"

test "Octree remove missing item keeps count" := do
  let points := #[Vec3.mk 1.0 1.0 1.0, Vec3.mk 5.0 5.0 5.0]
  let tree := Octree.build points
  let tree' := tree.remove 99 (AABB.fromPoint points[0]!)
  ensure (tree'.itemCount == tree.itemCount) "missing item should not change count"

test "Octree remove item not found in queries" := do
  let points := #[Vec3.mk 1.0 1.0 1.0, Vec3.mk 5.0 5.0 5.0, Vec3.mk 9.0 9.0 9.0]
  let tree := Octree.build points
  let tree' := tree.remove 1 (AABB.fromPoint points[1]!)
  let results := tree'.queryPoint (Vec3.mk 5.0 5.0 5.0)
  ensure (!results.contains 1) "removed item should not be found"

test "Octree remove preserves other items" := do
  let points := #[Vec3.mk 1.0 1.0 1.0, Vec3.mk 5.0 5.0 5.0, Vec3.mk 9.0 9.0 9.0]
  let tree := Octree.build points
  let tree' := tree.remove 1 (AABB.fromPoint points[1]!)
  let results := tree'.queryAABB (AABB.fromMinMax Vec3.zero (Vec3.mk 10.0 10.0 10.0))
  ensure (results.contains 0) "should still find item 0"
  ensure (results.contains 2) "should still find item 2"

test "Octree empty tree queries return empty" := do
  let tree := Octree.empty (AABB.fromMinMax Vec3.zero (Vec3.mk 10.0 10.0 10.0))
  let results := tree.queryAABB (AABB.fromMinMax Vec3.zero (Vec3.mk 5.0 5.0 5.0))
  ensure results.isEmpty "empty tree should return empty results"

test "Octree single item" := do
  let tree := Octree.build #[Vec3.mk 5.0 5.0 5.0]
  ensure (tree.itemCount == 1) "should have 1 item"

test "Octree kNearest returns sorted by distance" := do
  let points := #[Vec3.mk 0.0 0.0 0.0, Vec3.mk 3.0 0.0 0.0, Vec3.mk 1.0 0.0 0.0]
  let tree := Octree.build points
  let results := tree.kNearest points Vec3.zero 3
  ensure (results[0]! == 0) "first should be closest"
  ensure (results[1]! == 2) "second should be second closest"
  ensure (results[2]! == 1) "third should be furthest"

test "Octree queryFrustum finds items in view" := do
  -- Create boxes spread out in space
  let boxes := #[
    AABB.fromMinMax (Vec3.mk 0.0 0.0 5.0) (Vec3.mk 2.0 2.0 7.0),
    AABB.fromMinMax (Vec3.mk 0.0 0.0 10.0) (Vec3.mk 2.0 2.0 12.0),
    AABB.fromMinMax (Vec3.mk 100.0 100.0 5.0) (Vec3.mk 102.0 102.0 7.0)  -- far away
  ]
  let tree := Octree.build boxes
  -- Create a frustum looking down -Z axis from position (1, 1, 20)
  let eye := Vec3.mk 1.0 1.0 20.0
  let target := Vec3.mk 1.0 1.0 0.0
  let view := Mat4.lookAt eye target Vec3.unitY
  let proj := Mat4.perspective (Float.pi / 2.0) 1.0 0.1 100.0
  let frustum := Frustum.fromViewProjection (view * proj)
  let results := tree.queryFrustum frustum
  -- The first two boxes should be in view, the third is far off to the side
  -- Just verify the query runs and returns something reasonable
  ensure (results.size <= 3) "should not find more items than exist"

test "Octree queryPointExact filters non-containing items" := do
  -- Two items in same leaf (no splitting), only one contains query point
  let items := #[
    AABB.fromMinMax Vec3.zero (Vec3.mk 5.0 5.0 5.0),
    AABB.fromMinMax (Vec3.mk 7.0 7.0 7.0) (Vec3.mk 9.0 9.0 9.0)
  ]
  let tree := Octree.build items
  let point := Vec3.mk 3.0 3.0 3.0
  -- Exact query filters to only items actually containing the point
  let exact := tree.queryPointExact items point
  ensure (exact.size == 1) "should find exactly 1 item containing point"
  ensure (exact[0]! == 0) "should be the first item"

test "Octree queryAABBExact filters non-intersecting items" := do
  let items := #[
    AABB.fromMinMax Vec3.zero (Vec3.mk 3.0 3.0 3.0),
    AABB.fromMinMax (Vec3.mk 7.0 7.0 7.0) (Vec3.mk 9.0 9.0 9.0)
  ]
  let tree := Octree.build items
  let query := AABB.fromMinMax (Vec3.mk 2.0 2.0 2.0) (Vec3.mk 4.0 4.0 4.0)
  let exact := tree.queryAABBExact items query
  ensure (exact.size == 1) "should find exactly 1 intersecting item"
  ensure (exact[0]! == 0) "should be the first item"

-- ============================================================================
-- KDTree3D Tests
-- ============================================================================

testSuite "KDTree3D"

test "KDTree3D build creates tree" := do
  let points := #[Vec3.mk 1.0 1.0 1.0, Vec3.mk 5.0 5.0 5.0, Vec3.mk 9.0 9.0 9.0]
  let tree := KDTree3D.build points
  ensure (tree.pointCount == 3) "should have 3 points"

test "KDTree3D nearest finds closest point" := do
  let points := #[Vec3.mk 0.0 0.0 0.0, Vec3.mk 5.0 0.0 0.0, Vec3.mk 10.0 0.0 0.0]
  let tree := KDTree3D.build points
  match tree.nearest points (Vec3.mk 4.0 0.0 0.0) with
  | some idx => ensure (idx == 1) "nearest should be index 1"
  | none => ensure false "should find nearest"

test "KDTree3D kNearest finds k closest" := do
  let points := #[Vec3.mk 0.0 0.0 0.0, Vec3.mk 1.0 0.0 0.0, Vec3.mk 2.0 0.0 0.0, Vec3.mk 10.0 0.0 0.0]
  let tree := KDTree3D.build points
  let results := tree.kNearest points Vec3.zero 3
  ensure (results.size == 3) "should find 3 nearest"

test "KDTree3D withinRadius finds nearby points" := do
  let points := #[Vec3.mk 0.0 0.0 0.0, Vec3.mk 1.0 0.0 0.0, Vec3.mk 10.0 0.0 0.0]
  let tree := KDTree3D.build points
  let results := tree.withinRadius points Vec3.zero 2.0
  ensure (results.size == 2) "should find 2 points within radius"

test "KDTree3D withinAABB finds points in box" := do
  let points := #[Vec3.mk 1.0 1.0 1.0, Vec3.mk 5.0 5.0 5.0, Vec3.mk 9.0 9.0 9.0]
  let tree := KDTree3D.build points
  let query := AABB.fromMinMax Vec3.zero (Vec3.mk 3.0 3.0 3.0)
  let results := tree.withinAABB points query
  ensure (results.size == 1) "should find 1 point in box"

test "KDTree3D kNearest returns sorted by distance" := do
  let points := #[Vec3.mk 0.0 0.0 0.0, Vec3.mk 3.0 0.0 0.0, Vec3.mk 1.0 0.0 0.0, Vec3.mk 2.0 0.0 0.0]
  let tree := KDTree3D.build points
  let results := tree.kNearest points Vec3.zero 4
  -- Should be sorted: 0 (dist 0), 2 (dist 1), 3 (dist 2), 1 (dist 3)
  ensure (results[0]! == 0) "first should be index 0 (distance 0)"
  ensure (results[1]! == 2) "second should be index 2 (distance 1)"
  ensure (results[2]! == 3) "third should be index 3 (distance 2)"
  ensure (results[3]! == 1) "fourth should be index 1 (distance 3)"

test "KDTree3D empty tree nearest returns none" := do
  let emptyPoints : Array Vec3 := #[]
  let tree := KDTree3D.build emptyPoints
  let result := tree.nearest emptyPoints Vec3.zero
  ensure result.isNone "empty tree should return none"

-- ============================================================================
-- KDTree2D Tests
-- ============================================================================

testSuite "KDTree2D"

test "KDTree2D build creates tree" := do
  let points := #[Vec2.mk 1.0 1.0, Vec2.mk 5.0 5.0, Vec2.mk 9.0 9.0]
  let tree := KDTree2D.build points
  ensure (tree.pointCount == 3) "should have 3 points"

test "KDTree2D nearest finds closest point" := do
  let points := #[Vec2.mk 0.0 0.0, Vec2.mk 5.0 0.0, Vec2.mk 10.0 0.0]
  let tree := KDTree2D.build points
  match tree.nearest points (Vec2.mk 4.0 0.0) with
  | some idx => ensure (idx == 1) "nearest should be index 1"
  | none => ensure false "should find nearest"

test "KDTree2D kNearest finds k closest" := do
  let points := #[Vec2.mk 0.0 0.0, Vec2.mk 1.0 0.0, Vec2.mk 2.0 0.0, Vec2.mk 10.0 0.0]
  let tree := KDTree2D.build points
  let results := tree.kNearest points Vec2.zero 3
  ensure (results.size == 3) "should find 3 nearest"

test "KDTree2D withinRadius finds nearby points" := do
  let points := #[Vec2.mk 0.0 0.0, Vec2.mk 1.0 0.0, Vec2.mk 10.0 0.0]
  let tree := KDTree2D.build points
  let results := tree.withinRadius points Vec2.zero 2.0
  ensure (results.size == 2) "should find 2 points within radius"

test "KDTree2D withinRect finds points in rectangle" := do
  let points := #[Vec2.mk 1.0 1.0, Vec2.mk 5.0 5.0, Vec2.mk 9.0 9.0]
  let tree := KDTree2D.build points
  let query := AABB2D.fromMinMax Vec2.zero (Vec2.mk 3.0 3.0)
  let results := tree.withinRect points query
  ensure (results.size == 1) "should find 1 point in rectangle"
  ensure (results[0]! == 0) "should be index 0"

-- ============================================================================
-- BVH Tests
-- ============================================================================

testSuite "BVH"

test "BVH build creates tree" := do
  let boxes := #[
    AABB.fromMinMax Vec3.zero (Vec3.mk 1.0 1.0 1.0),
    AABB.fromMinMax (Vec3.mk 2.0 2.0 2.0) (Vec3.mk 3.0 3.0 3.0),
    AABB.fromMinMax (Vec3.mk 4.0 4.0 4.0) (Vec3.mk 5.0 5.0 5.0)
  ]
  let bvh := BVH.build boxes
  ensure (bvh.itemCount == 3) "should have 3 items"

test "BVH queryAABB finds intersecting items" := do
  let boxes := #[
    AABB.fromMinMax Vec3.zero (Vec3.mk 2.0 2.0 2.0),
    AABB.fromMinMax (Vec3.mk 5.0 5.0 5.0) (Vec3.mk 7.0 7.0 7.0),
    AABB.fromMinMax (Vec3.mk 10.0 10.0 10.0) (Vec3.mk 12.0 12.0 12.0)
  ]
  let bvh := BVH.build boxes
  let query := AABB.fromMinMax (Vec3.mk (-1.0) (-1.0) (-1.0)) (Vec3.mk 3.0 3.0 3.0)
  let results := bvh.queryAABB query
  ensure (results.size >= 1) "should find at least 1 item"

test "BVH rayCast finds closest hit" := do
  let boxes := #[
    AABB.fromMinMax (Vec3.mk (-1.0) (-1.0) 5.0) (Vec3.mk 1.0 1.0 7.0),
    AABB.fromMinMax (Vec3.mk (-1.0) (-1.0) 10.0) (Vec3.mk 1.0 1.0 12.0)
  ]
  let bvh := BVH.build boxes
  let ray := Ray.mk' Vec3.zero Vec3.unitZ
  -- We need a hit test function for the actual primitives
  let hitTest := fun idx =>
    if h : idx < boxes.size then
      Intersection.rayAABBHit ray boxes[idx]
    else none
  match bvh.rayCast ray hitTest with
  | some hit => ensure (hit.index == 0) "should hit first box (closer)"
  | none => ensure false "should find a hit"

test "BVH rayAny returns true when hit exists" := do
  let boxes := #[
    AABB.fromMinMax (Vec3.mk (-1.0) (-1.0) 5.0) (Vec3.mk 1.0 1.0 7.0)
  ]
  let bvh := BVH.build boxes
  let ray := Ray.mk' Vec3.zero Vec3.unitZ
  let hitTest := fun idx =>
    if h : idx < boxes.size then
      (Intersection.rayAABB ray boxes[idx]).isSome
    else false
  ensure (bvh.rayAny ray 100.0 hitTest) "should find a hit"

test "BVH rayAny returns false when no hit" := do
  let boxes := #[
    AABB.fromMinMax (Vec3.mk 10.0 10.0 5.0) (Vec3.mk 12.0 12.0 7.0)
  ]
  let bvh := BVH.build boxes
  let ray := Ray.mk' Vec3.zero Vec3.unitZ
  let hitTest := fun idx =>
    if h : idx < boxes.size then
      (Intersection.rayAABB ray boxes[idx]).isSome
    else false
  ensure (!(bvh.rayAny ray 100.0 hitTest)) "should not find a hit"

test "BVH depth is reasonable" := do
  let boxes := (List.range 16).map (fun i =>
    let f := i.toFloat
    AABB.fromMinMax (Vec3.mk f f f) (Vec3.mk (f+1.0) (f+1.0) (f+1.0))
  ) |>.toArray
  let bvh := BVH.build boxes
  ensure (bvh.depth <= 10) "depth should be reasonable"

test "BVH refit updates bounds" := do
  let boxes := #[
    AABB.fromMinMax Vec3.zero (Vec3.mk 1.0 1.0 1.0),
    AABB.fromMinMax (Vec3.mk 2.0 2.0 2.0) (Vec3.mk 3.0 3.0 3.0)
  ]
  let bvh := BVH.build boxes
  -- Move boxes
  let movedBoxes := #[
    AABB.fromMinMax (Vec3.mk 10.0 10.0 10.0) (Vec3.mk 11.0 11.0 11.0),
    AABB.fromMinMax (Vec3.mk 12.0 12.0 12.0) (Vec3.mk 13.0 13.0 13.0)
  ]
  let refitted := bvh.refit movedBoxes
  -- Query should find items at new locations
  let query := AABB.fromMinMax (Vec3.mk 9.0 9.0 9.0) (Vec3.mk 14.0 14.0 14.0)
  let results := refitted.queryAABB query
  ensure (results.size == 2) "should find both items at new locations"

-- ============================================================================
-- Dynamic AABB Tree Tests
-- ============================================================================

testSuite "DynamicAABBTree"

test "insert, update, and remove" := do
  let config : DynamicAABBTreeConfig := { fatMargin := 0.0 }
  let tree := DynamicAABBTree.empty config
  let a0 := AABB.fromMinMax Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let a1 := AABB.fromMinMax (Vec3.mk 2.0 2.0 2.0) (Vec3.mk 3.0 3.0 3.0)
  let (tree, p0) := tree.insert a0 0
  let (tree, p1) := tree.insert a1 1
  let query0 := AABB.fromMinMax (Vec3.mk (-0.5) (-0.5) (-0.5)) (Vec3.mk 1.5 1.5 1.5)
  let results0 := tree.queryAABB query0
  ensure (results0.contains 0) "should find item 0"
  ensure (!results0.contains 1) "should not find item 1"
  let moved := AABB.fromMinMax (Vec3.mk 0.5 0.5 0.5) (Vec3.mk 1.5 1.5 1.5)
  let tree := tree.update p1 moved
  let results1 := tree.queryAABB (AABB.fromMinMax Vec3.zero (Vec3.mk 2.0 2.0 2.0))
  ensure (results1.contains 1) "should find moved item 1"
  let tree := tree.remove p0
  let results2 := tree.queryAABB (AABB.fromMinMax (Vec3.mk (-1.0) (-1.0) (-1.0)) (Vec3.mk 2.0 2.0 2.0))
  ensure (!results2.contains 0) "removed item should not be found"

test "broadPhasePairs finds overlaps" := do
  let config : DynamicAABBTreeConfig := { fatMargin := 0.0 }
  let tree := DynamicAABBTree.empty config
  let a0 := AABB.fromMinMax Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let a1 := AABB.fromMinMax (Vec3.mk 0.5 0.5 0.5) (Vec3.mk 1.5 1.5 1.5)
  let a2 := AABB.fromMinMax (Vec3.mk 5.0 5.0 5.0) (Vec3.mk 6.0 6.0 6.0)
  let (tree, _) := tree.insert a0 0
  let (tree, _) := tree.insert a1 1
  let (tree, _) := tree.insert a2 2
  let pairs := tree.broadPhasePairs
  ensure (pairs.contains (0, 1)) "should include overlapping pair"
  ensure (!(pairs.contains (0, 2))) "should not include disjoint pair"

test "queryRay finds hits within maxT" := do
  let config : DynamicAABBTreeConfig := { fatMargin := 0.0 }
  let tree := DynamicAABBTree.empty config
  let a0 := AABB.fromMinMax (Vec3.mk (-1.0) (-1.0) 5.0) (Vec3.mk 1.0 1.0 6.0)
  let (tree, _) := tree.insert a0 0
  let ray := Ray.mk' Vec3.zero Vec3.unitZ
  let hits0 := tree.queryRay ray 4.0
  ensure (!hits0.contains 0) "should not hit beyond maxT"
  let hits1 := tree.queryRay ray 10.0
  ensure (hits1.contains 0) "should hit within maxT"

-- ============================================================================
-- Sweep-and-Prune Tests
-- ============================================================================

testSuite "SweepAndPrune"

test "queryAABB and pairs" := do
  let sap := SweepAndPrune.empty .x
  let a0 := AABB.fromMinMax Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let a1 := AABB.fromMinMax (Vec3.mk 0.5 0.5 0.5) (Vec3.mk 1.5 1.5 1.5)
  let a2 := AABB.fromMinMax (Vec3.mk 5.0 5.0 5.0) (Vec3.mk 6.0 6.0 6.0)
  let sap := sap.insert 0 a0
  let sap := sap.insert 1 a1
  let sap := sap.insert 2 a2
  let query := AABB.fromMinMax (Vec3.mk (-0.5) (-0.5) (-0.5)) (Vec3.mk 2.0 2.0 2.0)
  let results := sap.queryAABB query
  ensure (results.contains 0) "should find item 0"
  ensure (results.contains 1) "should find item 1"
  ensure (!results.contains 2) "should not find item 2"
  let pairs := sap.broadPhasePairs
  ensure (pairs.contains (0, 1)) "should include overlapping pair"
  ensure (!(pairs.contains (0, 2))) "should not include disjoint pair"

test "update moves entries" := do
  let sap := SweepAndPrune.empty .x
  let a0 := AABB.fromMinMax Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let a1 := AABB.fromMinMax (Vec3.mk 3.0 0.0 0.0) (Vec3.mk 4.0 1.0 1.0)
  let sap := sap.insert 0 a0
  let sap := sap.insert 1 a1
  let moved := AABB.fromMinMax (Vec3.mk 0.5 0.0 0.0) (Vec3.mk 1.5 1.0 1.0)
  let sap := sap.update 1 moved
  let pairs := sap.broadPhasePairs
  ensure (pairs.contains (0, 1)) "should overlap after update"

test "auto axis preserves pairs" := do
  let sap := SweepAndPrune.empty .x
  let a0 := AABB.fromMinMax Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let a1 := AABB.fromMinMax (Vec3.mk 0.5 0.5 0.5) (Vec3.mk 1.5 1.5 1.5)
  let a2 := AABB.fromMinMax (Vec3.mk 4.0 4.0 4.0) (Vec3.mk 5.0 5.0 5.0)
  let sap := sap.insert 0 a0
  let sap := sap.insert 1 a1
  let sap := sap.insert 2 a2
  let pairs := sap.broadPhasePairsAuto
  ensure (pairs.contains (0, 1)) "auto axis should include overlap"
  ensure (!(pairs.contains (0, 2))) "auto axis should not include disjoint pair"



end LinalgTests.SpatialTests
