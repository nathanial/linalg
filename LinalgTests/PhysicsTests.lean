/-
  Tests for Physics helpers.
-/

import Linalg
import Crucible

namespace LinalgTests.PhysicsTests

open Crucible
open Linalg

-- ============================================================================
-- Particle Tests
-- ============================================================================

testSuite "Particle"

test "create particle with mass" := do
  let p := Particle.create (Vec3.mk 1.0 2.0 3.0) Vec3.zero Vec3.zero 2.0
  ensure (floatNear p.mass 2.0 0.0001) "mass should be 2"
  ensure (floatNear p.inverseMass 0.5 0.0001) "inverse mass should be 0.5"
  ensure (floatNear p.position.x 1.0 0.0001) "position x"

test "static particle has zero inverse mass" := do
  let p := Particle.static (Vec3.mk 5.0 0.0 0.0)
  ensure (floatNear p.inverseMass 0.0 0.0001) "inverse mass should be 0"
  ensure p.isStatic "should be static"

test "apply force to particle" := do
  let p := Particle.create Vec3.zero Vec3.zero Vec3.zero 2.0
  let p' := p.applyForce (Vec3.mk 10.0 0.0 0.0)
  -- a = F/m = 10/2 = 5
  ensure (floatNear p'.acceleration.x 5.0 0.0001) "acceleration x should be 5"

test "apply force to static particle does nothing" := do
  let p := Particle.static Vec3.zero
  let p' := p.applyForce (Vec3.mk 100.0 0.0 0.0)
  ensure (floatNear p'.acceleration.x 0.0 0.0001) "static particle should not accelerate"

-- ============================================================================
-- Rigid Body Tests
-- ============================================================================

testSuite "Rigid Body"

test "create rigid body" := do
  let inertia := InertiaTensor.solidSphere 1.0 1.0
  let body := RigidBody.create Vec3.zero 1.0 inertia
  ensure (floatNear body.mass 1.0 0.0001) "mass"
  ensure (floatNear body.inverseMass 1.0 0.0001) "inverse mass"

test "static rigid body" := do
  let body := RigidBody.static Vec3.zero
  ensure body.isStatic "should be static"
  ensure (floatNear body.inverseMass 0.0 0.0001) "inverse mass should be 0"

-- ============================================================================
-- Inertia Tensor Tests
-- ============================================================================

testSuite "Inertia Tensors"

test "solid sphere inertia is 2/5 * m * r^2" := do
  let mass := 10.0
  let radius := 2.0
  let expected := 0.4 * mass * radius * radius  -- 16.0
  let tensor := InertiaTensor.solidSphere mass radius
  ensure (floatNear (tensor.get 0 0) expected 0.0001) "Ixx"
  ensure (floatNear (tensor.get 1 1) expected 0.0001) "Iyy"
  ensure (floatNear (tensor.get 2 2) expected 0.0001) "Izz"

test "hollow sphere inertia is 2/3 * m * r^2" := do
  let mass := 10.0
  let radius := 2.0
  let expected := (2.0 / 3.0) * mass * radius * radius  -- 26.67
  let tensor := InertiaTensor.hollowSphere mass radius
  ensure (floatNear (tensor.get 0 0) expected 0.0001) "Ixx"
  ensure (floatNear (tensor.get 1 1) expected 0.0001) "Iyy"

test "solid box with equal sides has equal inertia on all axes" := do
  let mass := 12.0
  let halfExtent := 1.0
  let tensor := InertiaTensor.solidBox mass (Vec3.mk halfExtent halfExtent halfExtent)
  -- For a cube: I = 1/12 * m * (2h)^2 * 2 = 1/6 * m * (2h)^2
  -- With half-extent 1, full side is 2, so I = 1/12 * 12 * 8 = 8
  let expected := (1.0 / 12.0) * mass * 8.0
  ensure (floatNear (tensor.get 0 0) expected 0.0001) "Ixx"
  ensure (floatNear (tensor.get 1 1) expected 0.0001) "Iyy"
  ensure (floatNear (tensor.get 2 2) expected 0.0001) "Izz"

test "thin rod inertia is 1/12 * m * L^2" := do
  let mass := 6.0
  let length := 3.0
  let expected := (1.0 / 12.0) * mass * length * length  -- 4.5
  let tensor := InertiaTensor.thinRod mass length
  ensure (floatNear (tensor.get 0 0) expected 0.0001) "Ixx (perpendicular)"
  ensure (floatNear (tensor.get 1 1) 0.0 0.0001) "Iyy (along rod)"
  ensure (floatNear (tensor.get 2 2) expected 0.0001) "Izz (perpendicular)"

test "cylinder has different inertia along axis vs perpendicular" := do
  let mass := 10.0
  let radius := 1.0
  let height := 4.0
  let tensor := InertiaTensor.solidCylinder mass radius height
  let iy := 0.5 * mass * radius * radius  -- 5.0
  let ixz := (1.0 / 12.0) * mass * (3.0 * radius * radius + height * height)  -- 14.17
  ensure (floatNear (tensor.get 1 1) iy 0.0001) "Iy (axis)"
  ensure (floatNear (tensor.get 0 0) ixz 0.01) "Ix (perpendicular)"

test "parallel axis increases inertia" := do
  let original := InertiaTensor.solidSphere 1.0 1.0
  let shifted := InertiaTensor.parallelAxis original 1.0 (Vec3.mk 2.0 0.0 0.0)
  -- Shifting in x adds m*d^2 = 4 to Iyy and Izz
  ensure (shifted.get 1 1 > original.get 1 1) "Iyy should increase"
  ensure (shifted.get 2 2 > original.get 2 2) "Izz should increase"

-- ============================================================================
-- Integration Tests
-- ============================================================================

testSuite "Euler Integration"

test "euler step with zero dt is identity" := do
  let p := Particle.create (Vec3.mk 1.0 2.0 3.0) (Vec3.mk 10.0 0.0 0.0) Vec3.zero 1.0
  let p' := Integration.eulerStep p 0.0
  ensure (floatNear p'.position.x 1.0 0.0001) "position unchanged"
  ensure (floatNear p'.velocity.x 10.0 0.0001) "velocity unchanged"

test "euler step updates position" := do
  let p := Particle.create Vec3.zero (Vec3.mk 10.0 0.0 0.0) Vec3.zero 1.0
  let p' := Integration.eulerStep p 0.1
  -- x = x0 + v*dt = 0 + 10*0.1 = 1.0
  ensure (floatNear p'.position.x 1.0 0.0001) "position should be 1.0"

test "euler step updates velocity from acceleration" := do
  let p := Particle.create Vec3.zero Vec3.zero (Vec3.mk 10.0 0.0 0.0) 1.0
  let p' := Integration.eulerStep p 0.1
  -- v = v0 + a*dt = 0 + 10*0.1 = 1.0
  ensure (floatNear p'.velocity.x 1.0 0.0001) "velocity should be 1.0"

test "semi-implicit euler is more accurate for gravity" := do
  -- Free fall: after 1 second with g=10, should fall approximately 5 meters
  let gravity := Vec3.mk 0.0 (-10.0) 0.0
  let p := Particle.create Vec3.zero Vec3.zero gravity 1.0

  -- 10 steps of 0.1s each
  let mut particle := p
  for _ in [:10] do
    particle := Integration.semiImplicitEulerStep particle 0.1

  -- Analytical: y = -0.5 * g * t^2 = -0.5 * 10 * 1 = -5.0
  -- Semi-implicit Euler gives -5.5 (first-order method has O(dt) error)
  ensure (floatNear particle.position.y (-5.5) 0.1) "should fall approximately 5.5 meters (semi-implicit)"

testSuite "Verlet Integration"

test "verlet step with constant acceleration" := do
  -- Starting at origin, no initial velocity, acceleration = (2, 0, 0)
  let pos := Vec3.zero
  let prevPos := Vec3.zero  -- No initial velocity
  let accel := Vec3.mk 2.0 0.0 0.0
  let dt := 0.1

  let newPos := Integration.verletStep pos prevPos accel dt
  -- x = 2*0 - 0 + 2*0.01 = 0.02
  ensure (floatNear newPos.x 0.02 0.0001) "position after one step"

test "velocity verlet matches analytical for constant acceleration" := do
  let pos := Vec3.zero
  let vel := Vec3.zero
  let accel := Vec3.mk 10.0 0.0 0.0
  let getAccel := fun _ => accel

  let (newPos, newVel) := Integration.velocityVerletStep pos vel accel getAccel 0.1
  -- x = v*t + 0.5*a*t^2 = 0 + 0.5*10*0.01 = 0.05
  -- v = v + a*t = 0 + 10*0.1 = 1.0
  ensure (floatNear newPos.x 0.05 0.0001) "position"
  ensure (floatNear newVel.x 1.0 0.0001) "velocity"

testSuite "RK4 Integration"

test "rk4 with constant acceleration" := do
  let pos := Vec3.zero
  let vel := Vec3.zero
  let getAccel := fun _ _ => Vec3.mk 10.0 0.0 0.0

  let (newPos, newVel) := Integration.rk4Step pos vel getAccel 0.1
  -- Should be very close to analytical solution
  ensure (floatNear newPos.x 0.05 0.001) "position"
  ensure (floatNear newVel.x 1.0 0.001) "velocity"

testSuite "Rigid Body Integration"

test "integrate rigid body updates position" := do
  let inertia := InertiaTensor.solidSphere 1.0 1.0
  let mut body := RigidBody.create Vec3.zero 1.0 inertia
  body := { body with velocity := Vec3.mk 10.0 0.0 0.0 }

  let body' := Integration.integrateRigidBody body 0.1
  ensure (floatNear body'.position.x 1.0 0.0001) "position x"

test "integrate rigid body with angular velocity rotates" := do
  let inertia := InertiaTensor.solidSphere 1.0 1.0
  let mut body := RigidBody.create Vec3.zero 1.0 inertia
  body := { body with angularVelocity := Vec3.mk 0.0 1.0 0.0 }

  let body' := Integration.integrateRigidBody body 0.1
  -- Orientation should have changed
  ensure (body'.orientation.w != 1.0 || body'.orientation.y != 0.0) "orientation should change"

-- ============================================================================
-- Collision Response Tests
-- ============================================================================

testSuite "Particle Collision"

test "equal mass head-on collision exchanges velocities" := do
  let a := Particle.create (Vec3.mk (-1.0) 0.0 0.0) (Vec3.mk 5.0 0.0 0.0) Vec3.zero 1.0
  let b := Particle.create (Vec3.mk 1.0 0.0 0.0) (Vec3.mk (-5.0) 0.0 0.0) Vec3.zero 1.0
  let contact : Contact := { point := Vec3.zero, normal := Vec3.unitX, penetration := 0.0 }

  let (a', b') := CollisionResponse.resolveParticleCollision a b contact 1.0
  -- Perfect elastic collision between equal masses exchanges velocities
  ensure (floatNear a'.velocity.x (-5.0) 0.01) "a should have b's velocity"
  ensure (floatNear b'.velocity.x 5.0 0.01) "b should have a's velocity"

test "perfectly inelastic collision reduces relative velocity to zero" := do
  let a := Particle.create (Vec3.mk (-1.0) 0.0 0.0) (Vec3.mk 10.0 0.0 0.0) Vec3.zero 1.0
  let b := Particle.create (Vec3.mk 1.0 0.0 0.0) Vec3.zero Vec3.zero 1.0
  let contact : Contact := { point := Vec3.zero, normal := Vec3.unitX, penetration := 0.0 }

  let (a', b') := CollisionResponse.resolveParticleCollision a b contact 0.0
  -- After inelastic collision, both move together at 5 m/s (momentum conservation)
  ensure (floatNear a'.velocity.x 5.0 0.01) "a velocity"
  ensure (floatNear b'.velocity.x 5.0 0.01) "b velocity"

test "collision with static particle reflects moving particle" := do
  let a := Particle.create (Vec3.mk (-1.0) 0.0 0.0) (Vec3.mk 10.0 0.0 0.0) Vec3.zero 1.0
  let b := Particle.static (Vec3.mk 1.0 0.0 0.0)
  let contact : Contact := { point := Vec3.zero, normal := Vec3.unitX, penetration := 0.0 }

  let (a', b') := CollisionResponse.resolveParticleCollision a b contact 1.0
  -- Perfect elastic collision with immovable object reverses velocity
  ensure (floatNear a'.velocity.x (-10.0) 0.01) "a should bounce back"
  ensure (floatNear b'.velocity.x 0.0 0.0001) "b should stay still"

test "impulse is along contact normal" := do
  let a := Particle.create Vec3.zero (Vec3.mk 10.0 5.0 0.0) Vec3.zero 1.0
  let b := Particle.create (Vec3.mk 2.0 0.0 0.0) Vec3.zero Vec3.zero 1.0
  let contact : Contact := { point := Vec3.mk 1.0 0.0 0.0, normal := Vec3.unitX, penetration := 0.0 }

  let impulse := CollisionResponse.particleImpulse a b contact 1.0
  -- Impulse should only be in X direction (along normal)
  ensure (floatNear impulse.y 0.0 0.0001) "impulse y should be 0"
  ensure (floatNear impulse.z 0.0 0.0001) "impulse z should be 0"

testSuite "Positional Correction"

test "correction separates overlapping particles" := do
  let a := Particle.create (Vec3.mk (-0.05) 0.0 0.0) Vec3.zero Vec3.zero 1.0
  let b := Particle.create (Vec3.mk 0.05 0.0 0.0) Vec3.zero Vec3.zero 1.0
  let contact : Contact := { point := Vec3.zero, normal := Vec3.unitX, penetration := 0.1 }

  let (a', b') := CollisionResponse.positionalCorrection a b contact 1.0 0.0
  -- They should be pushed apart
  ensure (a'.position.x < a.position.x) "a should move left"
  ensure (b'.position.x > b.position.x) "b should move right"

test "correction with slop prevents jitter" := do
  let a := Particle.create (Vec3.mk (-0.05) 0.0 0.0) Vec3.zero Vec3.zero 1.0
  let b := Particle.create (Vec3.mk 0.05 0.0 0.0) Vec3.zero Vec3.zero 1.0
  let contact : Contact := { point := Vec3.zero, normal := Vec3.unitX, penetration := 0.005 }

  -- With slop of 0.01, penetration of 0.005 should not cause correction
  let (a', b') := CollisionResponse.positionalCorrection a b contact 0.8 0.01
  ensure (floatNear a'.position.x a.position.x 0.0001) "a should not move"
  ensure (floatNear b'.position.x b.position.x 0.0001) "b should not move"

-- ============================================================================
-- Swept Collision Tests
-- ============================================================================

testSuite "Swept Sphere Collision"

test "sphere moving toward sphere finds hit" := do
  let start := Vec3.mk (-5.0) 0.0 0.0
  let end_ := Vec3.mk 5.0 0.0 0.0
  let target := Sphere.mk' Vec3.zero 1.0

  match SweptCollision.sphereVsSphere start end_ 1.0 target with
  | some hit =>
    -- Should hit when centers are 2 units apart (at x = -2)
    -- t = (0 - (-5) - 2) / 10 = 3/10 = 0.3
    ensure (floatNear hit.t 0.3 0.01) "hit time"
  | none => ensure false "should find hit"

test "sphere moving away finds no hit" := do
  let start := Vec3.mk 5.0 0.0 0.0
  let end_ := Vec3.mk 10.0 0.0 0.0
  let target := Sphere.mk' Vec3.zero 1.0

  match SweptCollision.sphereVsSphere start end_ 1.0 target with
  | some _ => ensure false "should not hit"
  | none => pure ()

test "parallel movement misses" := do
  let start := Vec3.mk (-5.0) 5.0 0.0  -- Above the sphere
  let end_ := Vec3.mk 5.0 5.0 0.0
  let target := Sphere.mk' Vec3.zero 1.0

  match SweptCollision.sphereVsSphere start end_ 1.0 target with
  | some _ => ensure false "should miss"
  | none => pure ()

test "already overlapping returns t=0" := do
  let start := Vec3.mk 0.5 0.0 0.0  -- Inside target
  let end_ := Vec3.mk 5.0 0.0 0.0
  let target := Sphere.mk' Vec3.zero 1.0

  match SweptCollision.sphereVsSphere start end_ 1.0 target with
  | some hit =>
    ensure (floatNear hit.t 0.0 0.0001) "should return t=0 for overlap"
  | none => ensure false "should detect overlap"

testSuite "Swept Plane Collision"

test "sphere sweeping into plane finds hit" := do
  let start := Vec3.mk 0.0 5.0 0.0
  let end_ := Vec3.mk 0.0 (-5.0) 0.0
  let plane := Plane.xz  -- y = 0 plane

  match SweptCollision.sphereVsPlane start end_ 1.0 plane with
  | some hit =>
    -- Should hit when sphere center is at y=1 (radius above plane)
    -- t = (5 - 1) / 10 = 0.4
    ensure (floatNear hit.t 0.4 0.01) "hit time"
  | none => ensure false "should find hit"

testSuite "Swept AABB Collision"

test "sphere sweeping into AABB face" := do
  let start := Vec3.mk (-5.0) 0.0 0.0
  let end_ := Vec3.mk 0.0 0.0 0.0
  let aabb := AABB.fromCenterExtents Vec3.zero (Vec3.mk 1.0 1.0 1.0)

  match SweptCollision.sphereVsAABB start end_ 0.5 aabb with
  | some hit =>
    ensure (hit.t > 0.0) "should hit after start"
    ensure (hit.t < 1.0) "should hit before end"
  | none => ensure false "should find hit"

-- ============================================================================
-- Continuous Collision Tests
-- ============================================================================

testSuite "Time of Impact"

test "sphere-sphere TOI matches analytical" := do
  let sphereA := Sphere.mk' (Vec3.mk (-5.0) 0.0 0.0) 1.0
  let velocityA := Vec3.mk 10.0 0.0 0.0
  let sphereB := Sphere.mk' (Vec3.mk 5.0 0.0 0.0) 1.0
  let velocityB := Vec3.zero

  match ContinuousCollision.sphereVsSphereTOI sphereA velocityA sphereB velocityB 1.0 with
  | some toi =>
    -- Distance = 10, combined radius = 2, relative velocity = 10
    -- TOI = (10 - 2) / 10 = 0.8
    ensure (floatNear toi 0.8 0.01) "TOI"
  | none => ensure false "should find TOI"

test "separating objects return none" := do
  let sphereA := Sphere.mk' (Vec3.mk (-5.0) 0.0 0.0) 1.0
  let velocityA := Vec3.mk (-10.0) 0.0 0.0  -- Moving away
  let sphereB := Sphere.mk' (Vec3.mk 5.0 0.0 0.0) 1.0
  let velocityB := Vec3.zero

  match ContinuousCollision.sphereVsSphereTOI sphereA velocityA sphereB velocityB 1.0 with
  | some _ => ensure false "should not find TOI for separating objects"
  | none => pure ()

test "stationary objects that don't overlap return none" := do
  let sphereA := Sphere.mk' (Vec3.mk (-5.0) 0.0 0.0) 1.0
  let sphereB := Sphere.mk' (Vec3.mk 5.0 0.0 0.0) 1.0

  match ContinuousCollision.sphereVsSphereTOI sphereA Vec3.zero sphereB Vec3.zero 1.0 with
  | some _ => ensure false "should not find TOI for stationary non-overlapping objects"
  | none => pure ()

test "already overlapping returns 0" := do
  let sphereA := Sphere.mk' Vec3.zero 1.0
  let sphereB := Sphere.mk' (Vec3.mk 0.5 0.0 0.0) 1.0

  match ContinuousCollision.sphereVsSphereTOI sphereA Vec3.zero sphereB Vec3.zero 1.0 with
  | some toi =>
    ensure (floatNear toi 0.0 0.0001) "overlapping should return 0"
  | none => ensure false "should detect overlap"

test "refineTOI improves accuracy" := do
  -- Simple distance function that crosses zero at t=0.5
  let distFunc := fun t => 0.5 - t
  let refined := ContinuousCollision.refineTOI distFunc 0.0 1.0 10
  ensure (floatNear refined 0.5 0.001) "should refine to 0.5"



end LinalgTests.PhysicsTests
