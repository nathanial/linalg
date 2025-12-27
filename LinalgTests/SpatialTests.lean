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

test "Quadtree queryCircle finds nearby points" := do
  let points := #[Vec2.mk 0.0 0.0, Vec2.mk 1.0 0.0, Vec2.mk 10.0 10.0]
  let tree := Quadtree.build points
  let results := tree.queryCircle Vec2.zero 2.0
  ensure (results.size >= 2) "should find at least 2 nearby points"

test "Quadtree nodeCount is reasonable" := do
  let points := (List.range 20).map (fun i => Vec2.mk i.toFloat i.toFloat) |>.toArray
  let tree := Quadtree.build points
  ensure (tree.nodeCount >= 1) "should have at least 1 node"

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

test "Octree rayCast finds items along ray" := do
  let boxes := #[
    AABB.fromMinMax (Vec3.mk (-1.0) (-1.0) 5.0) (Vec3.mk 1.0 1.0 7.0),
    AABB.fromMinMax (Vec3.mk (-1.0) (-1.0) 10.0) (Vec3.mk 1.0 1.0 12.0)
  ]
  let tree := Octree.build boxes
  let ray := Ray.mk' Vec3.zero Vec3.unitZ
  let results := tree.rayCast ray
  ensure (results.size == 2) "ray should hit both boxes"

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

#generate_tests

end LinalgTests.SpatialTests
