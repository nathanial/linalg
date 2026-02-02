/-
  Tests for geometry primitives and intersection tests.
-/

import Linalg
import Crucible

namespace LinalgTests.GeometryTests

open Crucible
open Linalg

testSuite "Ray"

test "ray pointAt works correctly" := do
  let ray := Ray.mk' Vec3.zero Vec3.unitZ
  let p := ray.pointAt 5.0
  ensure (floatNear p.z 5.0 0.0001) "z should be 5"

testSuite "AABB"

test "AABB center is correct" := do
  let aabb := AABB.fromMinMax (Vec3.mk 0.0 0.0 0.0) (Vec3.mk 10.0 10.0 10.0)
  let c := aabb.center
  ensure (floatNear c.x 5.0 0.0001) "center x should be 5"
  ensure (floatNear c.y 5.0 0.0001) "center y should be 5"
  ensure (floatNear c.z 5.0 0.0001) "center z should be 5"

test "AABB contains center point" := do
  let aabb := AABB.fromMinMax Vec3.zero (Vec3.mk 10.0 10.0 10.0)
  ensure (aabb.containsPoint (Vec3.mk 5.0 5.0 5.0)) "should contain center"

test "AABB does not contain point outside" := do
  let aabb := AABB.fromMinMax Vec3.zero (Vec3.mk 10.0 10.0 10.0)
  ensure (!(aabb.containsPoint (Vec3.mk 15.0 5.0 5.0))) "should not contain outside point"

testSuite "Sphere"

test "sphere contains its center" := do
  let s := Sphere.mk' (Vec3.mk 5.0 5.0 5.0) 3.0
  ensure (s.containsPoint s.center) "should contain center"

test "sphere does not contain distant point" := do
  let s := Sphere.mk' Vec3.zero 1.0
  ensure (!(s.containsPoint (Vec3.mk 10.0 0.0 0.0))) "should not contain distant point"

testSuite "Plane"

test "point on XY plane has zero distance" := do
  let p := Plane.xy
  let point := Vec3.mk 5.0 5.0 0.0
  ensure (floatNear (p.distanceToPoint point) 0.0 0.0001) "distance should be 0"

test "point above XY plane has positive signed distance" := do
  let p := Plane.xy
  let point := Vec3.mk 0.0 0.0 5.0
  ensure (floatNear (p.signedDistance point) 5.0 0.0001) "signed distance should be 5"

testSuite "Intersection"

test "ray hits sphere from outside" := do
  let ray := Ray.mk' (Vec3.mk 0.0 0.0 (-10.0)) Vec3.unitZ
  let sphere := Sphere.mk' Vec3.zero 1.0
  match Intersection.raySphere ray sphere with
  | some hit => ensure (floatNear hit.t 9.0 0.0001) "t should be 9"
  | none => ensure false "expected ray to hit sphere"

test "ray misses sphere" := do
  let ray := Ray.mk' (Vec3.mk 10.0 0.0 (-10.0)) Vec3.unitZ
  let sphere := Sphere.mk' Vec3.zero 1.0
  match Intersection.raySphere ray sphere with
  | some _ => ensure false "expected ray to miss sphere"
  | none => pure ()

test "ray hits AABB" := do
  let ray := Ray.mk' (Vec3.mk 0.0 0.0 (-10.0)) Vec3.unitZ
  let aabb := AABB.fromMinMax (Vec3.mk (-1.0) (-1.0) (-1.0)) (Vec3.mk 1.0 1.0 1.0)
  match Intersection.rayAABB ray aabb with
  | some (tMin, _) => ensure (floatNear tMin 9.0 0.0001) "tMin should be 9"
  | none => ensure false "expected ray to hit AABB"

test "ray hits plane" := do
  let ray := Ray.mk' (Vec3.mk 0.0 0.0 (-5.0)) Vec3.unitZ
  let plane := Plane.xy
  match Intersection.rayPlane ray plane with
  | some hit => ensure (floatNear hit.t 5.0 0.0001) "t should be 5"
  | none => ensure false "expected ray to hit plane"

test "sphere-sphere intersection" := do
  let a := Sphere.mk' Vec3.zero 2.0
  let b := Sphere.mk' (Vec3.mk 3.0 0.0 0.0) 2.0
  ensure (Intersection.sphereSphere a b) "spheres should intersect"

test "sphere-sphere no intersection" := do
  let a := Sphere.mk' Vec3.zero 1.0
  let b := Sphere.mk' (Vec3.mk 10.0 0.0 0.0) 1.0
  ensure (!(Intersection.sphereSphere a b)) "spheres should not intersect"

test "AABB-AABB intersection" := do
  let a := AABB.fromMinMax Vec3.zero (Vec3.mk 5.0 5.0 5.0)
  let b := AABB.fromMinMax (Vec3.mk 3.0 3.0 3.0) (Vec3.mk 8.0 8.0 8.0)
  ensure (Intersection.aabbAABB a b) "AABBs should intersect"

test "AABB-AABB no intersection" := do
  let a := AABB.fromMinMax Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let b := AABB.fromMinMax (Vec3.mk 5.0 5.0 5.0) (Vec3.mk 6.0 6.0 6.0)
  ensure (!(Intersection.aabbAABB a b)) "AABBs should not intersect"

testSuite "Mesh"

test "mesh rayCast hits triangle" := do
  let vertices := #[
    Vec3.mk 0.0 0.0 0.0,
    Vec3.mk 1.0 0.0 0.0,
    Vec3.mk 0.0 1.0 0.0
  ]
  let indices := #[0, 1, 2]
  let mesh : Mesh := { vertices, indices }
  let ray := Ray.mk' (Vec3.mk 0.25 0.25 (-1.0)) Vec3.unitZ
  match Mesh.rayCast mesh ray with
  | some hit => ensure (hit.t > 0.0) "expected hit with positive t"
  | none => ensure false "expected ray to hit mesh"

test "mesh BVH rayCast hits triangle" := do
  let vertices := #[
    Vec3.mk 0.0 0.0 0.0,
    Vec3.mk 1.0 0.0 0.0,
    Vec3.mk 0.0 1.0 0.0
  ]
  let indices := #[0, 1, 2]
  let mesh : Mesh := { vertices, indices }
  let ray := Ray.mk' (Vec3.mk 0.25 0.25 (-1.0)) Vec3.unitZ
  match Linalg.Spatial.MeshBVH.build mesh with
  | none => ensure false "expected BVH to build"
  | some bvh =>
    match Linalg.Spatial.MeshBVH.rayCast mesh bvh ray with
    | some hit => ensure (hit.t > 0.0) "expected BVH hit with positive t"
    | none => ensure false "expected BVH ray to hit mesh"

testSuite "ConvexHull3D"

test "quickhull builds tetrahedron" := do
  let points := #[
    Vec3.mk 0.0 0.0 0.0,
    Vec3.mk 1.0 0.0 0.0,
    Vec3.mk 0.0 1.0 0.0,
    Vec3.mk 0.0 0.0 1.0
  ]
  let hull := ConvexHull3D.quickHull points
  ensure (hull.faceCount == 4) "tetrahedron should have 4 faces"

testSuite "Convex Decomposition"

test "decomposition splits separated triangles" := do
  let vertices := #[
    Vec3.mk 0.0 0.0 0.0,
    Vec3.mk 1.0 0.0 0.0,
    Vec3.mk 0.0 1.0 0.0,
    Vec3.mk 10.0 0.0 0.0,
    Vec3.mk 11.0 0.0 0.0,
    Vec3.mk 10.0 1.0 0.0
  ]
  let indices := #[0, 1, 2, 3, 4, 5]
  let mesh : Mesh := { vertices, indices }
  let config : ConvexDecompositionConfig := { maxTrianglesPerPart := 1, maxDepth := 4 }
  let parts := ConvexDecomposition.decompose mesh config
  ensure (parts.size == 2) "should split into 2 parts"
  if parts.size == 2 then
    let c0 := parts[0]!.bounds.center
    let c1 := parts[1]!.bounds.center
    ensure (c0.distance c1 > 5.0) "parts should be far apart"
  else
    pure ()

test "decomposition of empty mesh is empty" := do
  let mesh : Mesh := { vertices := #[], indices := #[] }
  let parts := ConvexDecomposition.decompose mesh
  ensure (parts.isEmpty) "empty mesh should produce no parts"

testSuite "Collision3D"

test "GJK3D AABB-AABB intersects" := do
  let a := AABB.fromMinMax Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let b := AABB.fromMinMax (Vec3.mk 0.5 0.5 0.5) (Vec3.mk 1.5 1.5 1.5)
  ensure (intersectsGJK3D a b) "AABBs should intersect"

test "GJK3D AABB-AABB no intersection" := do
  let a := AABB.fromMinMax Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let b := AABB.fromMinMax (Vec3.mk 3.0 3.0 3.0) (Vec3.mk 4.0 4.0 4.0)
  ensure (!(intersectsGJK3D a b)) "AABBs should not intersect"



end LinalgTests.GeometryTests
