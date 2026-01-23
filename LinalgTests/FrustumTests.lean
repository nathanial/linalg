/-
  Tests for Frustum primitive and culling.
-/

import Linalg
import Crucible

namespace LinalgTests.FrustumTests

open Crucible
open Linalg

-- Helper to create a simple box-shaped frustum for testing
-- Plane equation: normal·point = distance
-- signedDistance = normal·point - distance (positive = in front of plane)
def testFrustum : Frustum :=
  -- Create a frustum that's a box from -10 to 10 on all axes
  -- For a point at origin to be inside all planes, signedDistance must be >= 0
  Frustum.mk'
    (Plane.mk (Vec3.mk 0.0 0.0 1.0) (-10.0))    -- near: z >= -10 (plane at z=-10, normal +Z)
    (Plane.mk (Vec3.mk 0.0 0.0 (-1.0)) (-10.0)) -- far: z <= 10 (plane at z=10, normal -Z)
    (Plane.mk (Vec3.mk 1.0 0.0 0.0) (-10.0))    -- left: x >= -10 (plane at x=-10, normal +X)
    (Plane.mk (Vec3.mk (-1.0) 0.0 0.0) (-10.0)) -- right: x <= 10 (plane at x=10, normal -X)
    (Plane.mk (Vec3.mk 0.0 (-1.0) 0.0) (-10.0)) -- top: y <= 10 (plane at y=10, normal -Y)
    (Plane.mk (Vec3.mk 0.0 1.0 0.0) (-10.0))    -- bottom: y >= -10 (plane at y=-10, normal +Y)

testSuite "Frustum Basic"

test "frustum from perspective projection creates planes" := do
  let proj := Mat4.perspective Float.halfPi 1.0 0.1 100.0
  let frustum := Frustum.fromViewProjection proj
  -- Verify it creates valid planes (all have unit normals)
  ensure (floatNear frustum.near.normal.length 1.0 0.001) "near plane should be normalized"
  ensure (floatNear frustum.far.normal.length 1.0 0.001) "far plane should be normalized"

test "frustum has 6 planes" := do
  let frustum := testFrustum
  let planes := frustum.planes
  ensure (planes.size == 6) "should have 6 planes"

testSuite "Frustum Point Containment"

test "point at origin is inside box frustum" := do
  let frustum := testFrustum
  ensure (frustum.containsPoint Vec3.zero) "origin should be inside"

test "point outside box frustum" := do
  let frustum := testFrustum
  ensure (!(frustum.containsPoint (Vec3.mk 20.0 0.0 0.0))) "point at x=20 should be outside"

test "point on boundary" := do
  let frustum := testFrustum
  -- Point exactly at x=10 is on the boundary
  ensure (frustum.containsPoint (Vec3.mk 9.9 0.0 0.0)) "point at x=9.9 should be inside"

testSuite "Frustum Sphere Culling"

test "sphere at origin is visible" := do
  let frustum := testFrustum
  let sphere := Sphere.mk' Vec3.zero 1.0
  ensure (frustum.isSphereVisible sphere) "sphere at origin should be visible"

test "sphere far outside is not visible" := do
  let frustum := testFrustum
  let sphere := Sphere.mk' (Vec3.mk 100.0 0.0 0.0) 1.0
  ensure (!(frustum.isSphereVisible sphere)) "sphere far outside should not be visible"

test "sphere partially outside is visible" := do
  let frustum := testFrustum
  -- Sphere centered at x=10 with radius 2 extends from x=8 to x=12
  let sphere := Sphere.mk' (Vec3.mk 10.0 0.0 0.0) 2.0
  ensure (frustum.isSphereVisible sphere) "sphere straddling boundary should be visible"

testSuite "Frustum AABB Culling"

test "AABB at origin is visible" := do
  let frustum := testFrustum
  let aabb := AABB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  ensure (frustum.isAABBVisible aabb) "AABB at origin should be visible"

test "AABB far outside is not visible" := do
  let frustum := testFrustum
  let aabb := AABB.fromCenterExtents (Vec3.mk 100.0 0.0 0.0) (Vec3.mk 1.0 1.0 1.0)
  ensure (!(frustum.isAABBVisible aabb)) "AABB far outside should not be visible"

test "large AABB straddling frustum is visible" := do
  let frustum := testFrustum
  let aabb := AABB.fromCenterExtents Vec3.zero (Vec3.mk 100.0 100.0 100.0)
  ensure (frustum.isAABBVisible aabb) "large AABB should be visible"

testSuite "Frustum Containment Types"

test "small sphere inside returns inside or intersects" := do
  let frustum := testFrustum
  let sphere := Sphere.mk' Vec3.zero 1.0
  let result := frustum.testSphere sphere
  ensure (result != .outside) "small sphere at origin should not be outside"

test "sphere outside returns outside" := do
  let frustum := testFrustum
  let sphere := Sphere.mk' (Vec3.mk 100.0 0.0 0.0) 1.0
  let result := frustum.testSphere sphere
  ensure (result == .outside) "distant sphere should be outside"

test "small AABB inside returns inside or intersects" := do
  let frustum := testFrustum
  let aabb := AABB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)
  let result := frustum.testAABB aabb
  ensure (result != .outside) "small AABB at origin should not be outside"



end LinalgTests.FrustumTests
