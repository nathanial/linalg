/-
  Physics helpers for game simulation.

  Includes velocity integration (Euler, Verlet, RK4), impulse-based
  collision response, inertia tensor calculations, and swept/continuous
  collision detection.
-/

import Linalg.Core
import Linalg.Vec3
import Linalg.Mat3
import Linalg.Quat
import Linalg.Geometry.Sphere
import Linalg.Geometry.AABB
import Linalg.Geometry.Plane

namespace Linalg

-- ============================================================================
-- Core Structures
-- ============================================================================

/-- State of a point mass particle for physics simulation. -/
structure Particle where
  position : Vec3
  velocity : Vec3
  acceleration : Vec3
  mass : Float
  inverseMass : Float
  deriving Repr, Inhabited

namespace Particle

/-- Create a particle with the given mass. Computes inverseMass automatically. -/
def create (position : Vec3) (velocity : Vec3 := Vec3.zero)
    (acceleration : Vec3 := Vec3.zero) (mass : Float := 1.0) : Particle :=
  let invMass := if mass > 0.0 then 1.0 / mass else 0.0
  { position, velocity, acceleration, mass, inverseMass := invMass }

/-- Create a static (immovable) particle with infinite mass. -/
def static (position : Vec3) : Particle :=
  { position, velocity := Vec3.zero, acceleration := Vec3.zero,
    mass := 0.0, inverseMass := 0.0 }

/-- Check if particle is static (infinite mass). -/
@[inline]
def isStatic (p : Particle) : Bool := p.inverseMass == 0.0

/-- Apply a force to the particle (F = ma, so a = F/m). -/
@[inline]
def applyForce (p : Particle) (force : Vec3) : Particle :=
  if p.isStatic then p
  else { p with acceleration := p.acceleration.add (force.scale p.inverseMass) }

/-- Clear accumulated acceleration. -/
@[inline]
def clearAcceleration (p : Particle) : Particle :=
  { p with acceleration := Vec3.zero }

end Particle

/-- State of a rigid body for physics simulation. -/
structure RigidBody where
  -- Linear properties
  position : Vec3
  velocity : Vec3
  acceleration : Vec3
  mass : Float
  inverseMass : Float
  -- Angular properties
  orientation : Quat
  angularVelocity : Vec3
  angularAcceleration : Vec3
  -- Inertia properties
  inertiaTensor : Mat3
  inverseInertiaTensor : Mat3
  deriving Repr, Inhabited

namespace RigidBody

/-- Create a rigid body. -/
def create (position : Vec3) (mass : Float) (inertiaTensor : Mat3)
    (orientation : Quat := Quat.identity) : RigidBody :=
  let invMass := if mass > 0.0 then 1.0 / mass else 0.0
  let invInertia := if mass > 0.0 then inertiaTensor.inverse.getD Mat3.zero else Mat3.zero
  { position
    velocity := Vec3.zero
    acceleration := Vec3.zero
    mass
    inverseMass := invMass
    orientation
    angularVelocity := Vec3.zero
    angularAcceleration := Vec3.zero
    inertiaTensor
    inverseInertiaTensor := invInertia }

/-- Create a static (immovable) rigid body. -/
def static (position : Vec3) (orientation : Quat := Quat.identity) : RigidBody :=
  { position
    velocity := Vec3.zero
    acceleration := Vec3.zero
    mass := 0.0
    inverseMass := 0.0
    orientation
    angularVelocity := Vec3.zero
    angularAcceleration := Vec3.zero
    inertiaTensor := Mat3.zero
    inverseInertiaTensor := Mat3.zero }

/-- Check if rigid body is static (infinite mass). -/
@[inline]
def isStatic (b : RigidBody) : Bool := b.inverseMass == 0.0

/-- Apply a force at the center of mass. -/
@[inline]
def applyForce (b : RigidBody) (force : Vec3) : RigidBody :=
  if b.isStatic then b
  else { b with acceleration := b.acceleration.add (force.scale b.inverseMass) }

/-- Apply a force at a world point, generating both linear and angular acceleration. -/
def applyForceAtPoint (b : RigidBody) (force : Vec3) (point : Vec3) : RigidBody :=
  if b.isStatic then b
  else
    let newAccel := b.acceleration.add (force.scale b.inverseMass)
    let r := point.sub b.position
    let torque := r.cross force
    let angAccel := b.inverseInertiaTensor.transformVec3 torque
    { b with
      acceleration := newAccel
      angularAcceleration := b.angularAcceleration.add angAccel }

/-- Apply a torque. -/
def applyTorque (b : RigidBody) (torque : Vec3) : RigidBody :=
  if b.isStatic then b
  else
    let angAccel := b.inverseInertiaTensor.transformVec3 torque
    { b with angularAcceleration := b.angularAcceleration.add angAccel }

/-- Clear accumulated accelerations. -/
def clearAccelerations (b : RigidBody) : RigidBody :=
  { b with acceleration := Vec3.zero, angularAcceleration := Vec3.zero }

end RigidBody

/-- Contact information from collision detection. -/
structure Contact where
  point : Vec3
  normal : Vec3   -- Points from A to B
  penetration : Float
  deriving Repr, Inhabited, BEq

/-- Material properties for collision response. -/
structure Material where
  restitution : Float := 0.5       -- Bounciness (0 = inelastic, 1 = elastic)
  staticFriction : Float := 0.6
  dynamicFriction : Float := 0.4
  deriving Repr, Inhabited

namespace Material

/-- Perfectly inelastic material (no bounce). -/
def inelastic : Material := { restitution := 0.0, staticFriction := 0.6, dynamicFriction := 0.4 }

/-- Bouncy material. -/
def bouncy : Material := { restitution := 0.9, staticFriction := 0.3, dynamicFriction := 0.2 }

/-- Icy/slippery material. -/
def slippery : Material := { restitution := 0.1, staticFriction := 0.05, dynamicFriction := 0.02 }

/-- Combine two materials (average). -/
def combine (a b : Material) : Material :=
  { restitution := (a.restitution + b.restitution) * 0.5
    staticFriction := Float.sqrt (a.staticFriction * b.staticFriction)
    dynamicFriction := Float.sqrt (a.dynamicFriction * b.dynamicFriction) }

end Material

-- ============================================================================
-- Inertia Tensors
-- ============================================================================

namespace InertiaTensor

/-- Create a diagonal Mat3. -/
private def diagonalMat3 (ix iy iz : Float) : Mat3 :=
  { data := #[ix, 0, 0, 0, iy, 0, 0, 0, iz] }

/-- Inertia tensor for a solid sphere. I = 2/5 * m * r^2 -/
def solidSphere (mass radius : Float) : Mat3 :=
  let i := 0.4 * mass * radius * radius
  diagonalMat3 i i i

/-- Inertia tensor for a hollow sphere (thin shell). I = 2/3 * m * r^2 -/
def hollowSphere (mass radius : Float) : Mat3 :=
  let i := (2.0 / 3.0) * mass * radius * radius
  diagonalMat3 i i i

/-- Inertia tensor for a solid box (cuboid).
    halfExtents are half-widths in each axis direction. -/
def solidBox (mass : Float) (halfExtents : Vec3) : Mat3 :=
  let x2 := 4.0 * halfExtents.x * halfExtents.x
  let y2 := 4.0 * halfExtents.y * halfExtents.y
  let z2 := 4.0 * halfExtents.z * halfExtents.z
  let ix := (1.0 / 12.0) * mass * (y2 + z2)
  let iy := (1.0 / 12.0) * mass * (x2 + z2)
  let iz := (1.0 / 12.0) * mass * (x2 + y2)
  diagonalMat3 ix iy iz

/-- Inertia tensor for a solid cylinder (axis along Y).
    height is total height, radius is cylinder radius. -/
def solidCylinder (mass radius height : Float) : Mat3 :=
  let r2 := radius * radius
  let h2 := height * height
  let iy := 0.5 * mass * r2                           -- Axis of rotation
  let ixz := (1.0 / 12.0) * mass * (3.0 * r2 + h2)   -- Perpendicular axes
  diagonalMat3 ixz iy ixz

/-- Inertia tensor for a solid cone (axis along Y, apex at origin).
    height is total height, radius is base radius. -/
def solidCone (mass radius height : Float) : Mat3 :=
  let r2 := radius * radius
  let h2 := height * height
  let iy := 0.3 * mass * r2                                    -- Axis
  let ixz := (3.0 / 80.0) * mass * (4.0 * r2 + h2)             -- Perpendicular
  diagonalMat3 ixz iy ixz

/-- Inertia tensor for a solid capsule (axis along Y).
    height is the cylinder part height (not including hemispherical caps).
    radius is the capsule radius. -/
def solidCapsule (mass radius height : Float) : Mat3 :=
  -- Split mass between cylinder and hemispheres
  let cylVolume := Float.pi * radius * radius * height
  let sphereVolume := (4.0 / 3.0) * Float.pi * radius * radius * radius
  let totalVolume := cylVolume + sphereVolume
  let cylMass := mass * cylVolume / totalVolume
  let sphereMass := mass * sphereVolume / totalVolume

  -- Cylinder contribution
  let cylR2 := radius * radius
  let cylH2 := height * height
  let cylIy := 0.5 * cylMass * cylR2
  let cylIxz := (1.0 / 12.0) * cylMass * (3.0 * cylR2 + cylH2)

  -- Hemisphere contributions (use parallel axis theorem)
  -- Each hemisphere is at distance (height/2 + 3r/8) from center
  let sphereILocal := 0.4 * sphereMass * cylR2  -- Sphere inertia
  let offsetY := height / 2.0 + 3.0 * radius / 8.0
  let sphereIy := sphereILocal  -- No offset along y axis
  let sphereIxz := sphereILocal + sphereMass * offsetY * offsetY  -- Parallel axis

  let iy := cylIy + sphereIy
  let ixz := cylIxz + sphereIxz
  diagonalMat3 ixz iy ixz

/-- Inertia tensor for a thin rod (axis along Y, centered at origin).
    I = 1/12 * m * L^2 for perpendicular axes, 0 for axis along rod. -/
def thinRod (mass length : Float) : Mat3 :=
  let i := (1.0 / 12.0) * mass * length * length
  diagonalMat3 i 0.0 i

/-- Transform inertia tensor from local to world space.
    worldTensor = R * localTensor * R^T -/
def toWorldSpace (localTensor : Mat3) (orientation : Quat) : Mat3 :=
  let rotMat := orientation.toMat3
  let rotMatT := rotMat.transpose
  rotMat.multiply localTensor |>.multiply rotMatT

/-- Parallel axis theorem: compute inertia tensor about a new axis.
    newTensor = tensor + m * (offset^T * offset * I - offset * offset^T)
    where offset is the displacement from original to new axis. -/
def parallelAxis (tensor : Mat3) (mass : Float) (offset : Vec3) : Mat3 :=
  -- Diagonal terms increase by m * (sum of other offsets squared)
  let offX2 := offset.x * offset.x
  let offY2 := offset.y * offset.y
  let offZ2 := offset.z * offset.z
  let addIx := mass * (offY2 + offZ2)
  let addIy := mass * (offX2 + offZ2)
  let addIz := mass * (offX2 + offY2)
  -- Off-diagonal terms
  let addXY := -mass * offset.x * offset.y
  let addXZ := -mass * offset.x * offset.z
  let addYZ := -mass * offset.y * offset.z

  -- Column-major order: [col0, col1, col2]
  let adjustment : Mat3 := { data := #[
    addIx, addXY, addXZ,  -- column 0
    addXY, addIy, addYZ,  -- column 1
    addXZ, addYZ, addIz   -- column 2
  ]}

  tensor.add adjustment

/-- Combine inertia tensors (simple addition, assumes same origin). -/
def combine (tensors : Array Mat3) : Mat3 :=
  tensors.foldl Mat3.add Mat3.zero

end InertiaTensor

-- ============================================================================
-- Velocity Integration
-- ============================================================================

namespace Integration

/-- Explicit (forward) Euler integration step for a particle.
    Simple but can be unstable for stiff systems.
    x(t+dt) = x(t) + v(t) * dt
    v(t+dt) = v(t) + a(t) * dt -/
def eulerStep (p : Particle) (dt : Float) : Particle :=
  { p with
    position := p.position.add (p.velocity.scale dt)
    velocity := p.velocity.add (p.acceleration.scale dt) }

/-- Semi-implicit (symplectic) Euler integration.
    Updates velocity first, then uses new velocity for position.
    More stable than explicit Euler, conserves energy better.
    v(t+dt) = v(t) + a(t) * dt
    x(t+dt) = x(t) + v(t+dt) * dt -/
def semiImplicitEulerStep (p : Particle) (dt : Float) : Particle :=
  let newVel := p.velocity.add (p.acceleration.scale dt)
  { p with
    velocity := newVel
    position := p.position.add (newVel.scale dt) }

/-- Stormer-Verlet integration step.
    Position-based, doesn't track velocity directly.
    x(t+dt) = 2*x(t) - x(t-dt) + a(t) * dt^2
    Returns new position. -/
def verletStep (position previousPosition acceleration : Vec3) (dt : Float) : Vec3 :=
  let dt2 := dt * dt
  position.scale 2.0 |>.sub previousPosition |>.add (acceleration.scale dt2)

/-- Velocity Verlet integration step.
    Provides both position and velocity, second-order accurate.
    x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
    v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
    getAcceleration takes position and returns acceleration. -/
def velocityVerletStep (position velocity acceleration : Vec3)
    (getAcceleration : Vec3 → Vec3) (dt : Float) : Vec3 × Vec3 :=
  let dt2 := dt * dt
  -- New position
  let newPos := position.add (velocity.scale dt) |>.add (acceleration.scale (0.5 * dt2))
  -- New acceleration at new position
  let newAccel := getAcceleration newPos
  -- New velocity using average of old and new acceleration
  let newVel := velocity.add ((acceleration.add newAccel).scale (0.5 * dt))
  (newPos, newVel)

/-- RK4 (Runge-Kutta 4th order) integration step.
    High accuracy but computationally expensive (4 force evaluations).
    getAcceleration takes (position, velocity) and returns acceleration. -/
def rk4Step (position velocity : Vec3)
    (getAcceleration : Vec3 → Vec3 → Vec3) (dt : Float) : Vec3 × Vec3 :=
  -- k1
  let a1 := getAcceleration position velocity
  let k1x := velocity
  let k1v := a1

  -- k2
  let p2 := position.add (k1x.scale (dt * 0.5))
  let v2 := velocity.add (k1v.scale (dt * 0.5))
  let a2 := getAcceleration p2 v2
  let k2x := v2
  let k2v := a2

  -- k3
  let p3 := position.add (k2x.scale (dt * 0.5))
  let v3 := velocity.add (k2v.scale (dt * 0.5))
  let a3 := getAcceleration p3 v3
  let k3x := v3
  let k3v := a3

  -- k4
  let p4 := position.add (k3x.scale dt)
  let v4 := velocity.add (k3v.scale dt)
  let a4 := getAcceleration p4 v4
  let k4x := v4
  let k4v := a4

  -- Weighted average
  let dx := (k1x.add (k2x.scale 2.0) |>.add (k3x.scale 2.0) |>.add k4x).scale (dt / 6.0)
  let dv := (k1v.add (k2v.scale 2.0) |>.add (k3v.scale 2.0) |>.add k4v).scale (dt / 6.0)

  (position.add dx, velocity.add dv)

/-- Integrate a rigid body using semi-implicit Euler.
    Handles both linear and angular motion. -/
def integrateRigidBody (body : RigidBody) (dt : Float) : RigidBody :=
  if body.isStatic then body
  else
    -- Linear integration
    let newVel := body.velocity.add (body.acceleration.scale dt)
    let newPos := body.position.add (newVel.scale dt)

    -- Angular integration
    let newAngVel := body.angularVelocity.add (body.angularAcceleration.scale dt)

    -- Update orientation: q' = q + 0.5 * dt * omega * q
    -- omega as pure quaternion (x, y, z, 0)
    let omega := Quat.mk newAngVel.x newAngVel.y newAngVel.z 0.0
    -- spin = omega * q (using HMul instance)
    let spinFull := omega * body.orientation
    -- Scale by 0.5 * dt and add to current orientation
    let s := 0.5 * dt
    let newOrientation := Quat.mk
      (body.orientation.x + spinFull.x * s)
      (body.orientation.y + spinFull.y * s)
      (body.orientation.z + spinFull.z * s)
      (body.orientation.w + spinFull.w * s)
      |>.normalize

    { body with
      position := newPos
      velocity := newVel
      orientation := newOrientation
      angularVelocity := newAngVel }

end Integration

-- ============================================================================
-- Collision Response
-- ============================================================================

namespace CollisionResponse

/-- Compute relative velocity at contact point between two particles. -/
@[inline]
def relativeVelocity (a b : Particle) : Vec3 :=
  b.velocity.sub a.velocity

/-- Compute impulse magnitude for collision between two particles.
    Returns impulse to apply to particle A (negate for B). -/
def particleImpulse (a b : Particle) (contact : Contact) (restitution : Float) : Vec3 :=
  let relVel := relativeVelocity a b
  let velAlongNormal := relVel.dot contact.normal

  -- Don't resolve if velocities are separating
  if velAlongNormal > 0.0 then Vec3.zero
  else
    -- Impulse magnitude: j = -(1 + e) * v_rel . n / (1/m_a + 1/m_b)
    let e := restitution
    let j := -(1.0 + e) * velAlongNormal / (a.inverseMass + b.inverseMass)
    contact.normal.scale j

/-- Apply an impulse to a particle. -/
@[inline]
def applyImpulse (p : Particle) (impulse : Vec3) : Particle :=
  if p.isStatic then p
  else { p with velocity := p.velocity.add (impulse.scale p.inverseMass) }

/-- Resolve collision between two particles.
    Returns updated particles. -/
def resolveParticleCollision (a b : Particle) (contact : Contact)
    (restitution : Float) : Particle × Particle :=
  let impulse := particleImpulse a b contact restitution
  let newA := applyImpulse a impulse.neg
  let newB := applyImpulse b impulse
  (newA, newB)

/-- Positional correction to prevent sinking (Baumgarte stabilization).
    percent: how much of the penetration to correct (0.2 - 0.8 typical)
    slop: penetration tolerance to prevent jitter -/
def positionalCorrection (a b : Particle) (contact : Contact)
    (percent : Float := 0.8) (slop : Float := 0.01) : Particle × Particle :=
  let penetration := Float.max (contact.penetration - slop) 0.0
  let totalInvMass := a.inverseMass + b.inverseMass
  if totalInvMass == 0.0 then (a, b)
  else
    let correction := contact.normal.scale (penetration / totalInvMass * percent)
    let newA := { a with position := a.position.sub (correction.scale a.inverseMass) }
    let newB := { b with position := b.position.add (correction.scale b.inverseMass) }
    (newA, newB)

/-- Compute impulse for rigid body collision.
    Returns (linearImpulse, angularImpulseA, angularImpulseB). -/
def rigidBodyImpulse (a b : RigidBody) (contact : Contact)
    (material : Material) : Vec3 × Vec3 × Vec3 :=
  -- Vectors from center of mass to contact point
  let ra := contact.point.sub a.position
  let rb := contact.point.sub b.position

  -- Relative velocity at contact point
  let relVel := b.velocity.add (b.angularVelocity.cross rb)
    |>.sub (a.velocity.add (a.angularVelocity.cross ra))

  let velAlongNormal := relVel.dot contact.normal

  -- Don't resolve if separating
  if velAlongNormal > 0.0 then (Vec3.zero, Vec3.zero, Vec3.zero)
  else
    -- Cross products for angular contribution
    let raCrossN := ra.cross contact.normal
    let rbCrossN := rb.cross contact.normal

    -- Effective mass denominator
    let angularA := a.inverseInertiaTensor.transformVec3 raCrossN |>.cross ra
    let angularB := b.inverseInertiaTensor.transformVec3 rbCrossN |>.cross rb
    let angularSum := (angularA.add angularB).dot contact.normal

    let denom := a.inverseMass + b.inverseMass + angularSum

    -- Impulse magnitude
    let j := -(1.0 + material.restitution) * velAlongNormal / denom

    let linearImpulse := contact.normal.scale j
    let angularImpulseA := ra.cross linearImpulse
    let angularImpulseB := rb.cross linearImpulse

    (linearImpulse, angularImpulseA, angularImpulseB)

/-- Apply linear and angular impulse to a rigid body. -/
def applyRigidBodyImpulse (body : RigidBody)
    (linearImpulse angularImpulse : Vec3) : RigidBody :=
  if body.isStatic then body
  else
    let newVel := body.velocity.add (linearImpulse.scale body.inverseMass)
    let newAngVel := body.angularVelocity.add (body.inverseInertiaTensor.transformVec3 angularImpulse)
    { body with velocity := newVel, angularVelocity := newAngVel }

/-- Resolve collision between two rigid bodies. -/
def resolveRigidBodyCollision (a b : RigidBody) (contact : Contact)
    (material : Material) : RigidBody × RigidBody :=
  let (impulse, angA, angB) := rigidBodyImpulse a b contact material
  let newA := applyRigidBodyImpulse a impulse.neg angA.neg
  let newB := applyRigidBodyImpulse b impulse angB
  (newA, newB)

/-- Compute friction impulse (tangent to contact normal). -/
def frictionImpulse (relativeVelocity contactNormal : Vec3)
    (normalImpulseMag : Float) (friction : Float) : Vec3 :=
  -- Get tangent velocity
  let velAlongNormal := relativeVelocity.dot contactNormal
  let tangentVel := relativeVelocity.sub (contactNormal.scale velAlongNormal)
  let tangentSpeed := tangentVel.length

  if tangentSpeed < 0.0001 then Vec3.zero
  else
    let tangent := tangentVel.scale (1.0 / tangentSpeed)
    -- Coulomb friction: |f_t| <= mu * |f_n|
    let maxFriction := friction * Float.abs normalImpulseMag
    let frictionMag := Float.min tangentSpeed maxFriction
    tangent.scale (-frictionMag)

end CollisionResponse

-- ============================================================================
-- Swept Collision Detection
-- ============================================================================

/-- Result of swept collision detection. -/
structure SweptHit where
  t : Float        -- Time of impact (0..1)
  point : Vec3     -- Contact point at time of impact
  normal : Vec3    -- Contact normal at time of impact
  deriving Repr, Inhabited

namespace SweptCollision

/-- Swept sphere vs static sphere collision.
    start/end_ are the center positions of the moving sphere.
    Returns hit info if collision occurs during the sweep. -/
def sphereVsSphere (start end_ : Vec3) (radius : Float)
    (target : Sphere) : Option SweptHit :=
  -- This is a ray-sphere intersection where the ray has thickness
  let movement := end_.sub start
  let toSphere := start.sub target.center
  let combinedRadius := radius + target.radius

  -- Quadratic coefficients: at^2 + bt + c = 0
  let a := movement.dot movement
  let b := 2.0 * toSphere.dot movement
  let c := toSphere.dot toSphere - combinedRadius * combinedRadius

  -- Already overlapping?
  if c < 0.0 then
    some { t := 0.0, point := start, normal := toSphere.normalize }
  else if a < 0.0001 then
    none  -- Not moving
  else
    let discriminant := b * b - 4.0 * a * c
    if discriminant < 0.0 then none
    else
      let t := (-b - Float.sqrt discriminant) / (2.0 * a)
      if t < 0.0 || t > 1.0 then none
      else
        let hitPos := start.add (movement.scale t)
        let normal := hitPos.sub target.center |>.normalize
        let contactPoint := target.center.add (normal.scale target.radius)
        some { t, point := contactPoint, normal }

/-- Swept sphere vs static plane collision. -/
def sphereVsPlane (start end_ : Vec3) (radius : Float)
    (plane : Plane) : Option SweptHit :=
  let movement := end_.sub start
  let startDist := plane.signedDistance start

  -- Check if starting inside (touching or penetrating)
  if Float.abs startDist <= radius then
    let contactPoint := start.sub (plane.normal.scale startDist)
    some { t := 0.0, point := contactPoint, normal := plane.normal }
  -- Moving parallel
  else
    let velocityAlongNormal := movement.dot plane.normal
    if Float.abs velocityAlongNormal < 0.0001 then none
    -- Moving away from plane
    else if (startDist > 0.0 && velocityAlongNormal > 0.0) ||
            (startDist < 0.0 && velocityAlongNormal < 0.0) then none
    else
      -- Find intersection time: startDist + t * velocityAlongNormal = ±radius
      let targetDist := if startDist > 0.0 then radius else -radius
      let t := (targetDist - startDist) / velocityAlongNormal
      if t < 0.0 || t > 1.0 then none
      else
        let hitPos := start.add (movement.scale t)
        let contactNormal := if startDist > 0.0 then plane.normal else plane.normal.neg
        let contactPoint := hitPos.sub (contactNormal.scale radius)
        some { t, point := contactPoint, normal := contactNormal }

/-- Swept sphere vs static AABB collision using Minkowski sum approach. -/
def sphereVsAABB (start end_ : Vec3) (radius : Float)
    (aabb : AABB) : Option SweptHit :=
  -- Expand the AABB by the sphere radius (Minkowski sum)
  let expandedMin := aabb.min.sub (Vec3.mk radius radius radius)
  let expandedMax := aabb.max.add (Vec3.mk radius radius radius)

  let movement := end_.sub start

  -- Ray vs expanded AABB (slab method)
  let invDx := if Float.abs movement.x > 0.0001 then 1.0 / movement.x else Float.infinity
  let invDy := if Float.abs movement.y > 0.0001 then 1.0 / movement.y else Float.infinity
  let invDz := if Float.abs movement.z > 0.0001 then 1.0 / movement.z else Float.infinity

  let tx1 := (expandedMin.x - start.x) * invDx
  let tx2 := (expandedMax.x - start.x) * invDx
  let ty1 := (expandedMin.y - start.y) * invDy
  let ty2 := (expandedMax.y - start.y) * invDy
  let tz1 := (expandedMin.z - start.z) * invDz
  let tz2 := (expandedMax.z - start.z) * invDz

  let tMinX := Float.min tx1 tx2
  let tMaxX := Float.max tx1 tx2
  let tMinY := Float.min ty1 ty2
  let tMaxY := Float.max ty1 ty2
  let tMinZ := Float.min tz1 tz2
  let tMaxZ := Float.max tz1 tz2

  let tMin := Float.max tMinX (Float.max tMinY tMinZ)
  let tMax := Float.min tMaxX (Float.min tMaxY tMaxZ)

  if tMax < 0.0 || tMin > tMax || tMin > 1.0 then none
  else
    let t := if tMin < 0.0 then 0.0 else tMin
    let hitPos := start.add (movement.scale t)

    -- Find closest point on original AABB to determine normal
    let closest := Vec3.mk
      (Float.clamp hitPos.x aabb.min.x aabb.max.x)
      (Float.clamp hitPos.y aabb.min.y aabb.max.y)
      (Float.clamp hitPos.z aabb.min.z aabb.max.z)

    let toHit := hitPos.sub closest
    let dist := toHit.length
    let normal := if dist > 0.0001 then toHit.scale (1.0 / dist)
      else
        -- Determine face from entry direction
        let ax := Float.abs (hitPos.x - aabb.center.x) / (aabb.extents.x + 0.0001)
        let ay := Float.abs (hitPos.y - aabb.center.y) / (aabb.extents.y + 0.0001)
        let az := Float.abs (hitPos.z - aabb.center.z) / (aabb.extents.z + 0.0001)
        if ax > ay && ax > az then
          if hitPos.x > aabb.center.x then Vec3.unitX else Vec3.left
        else if ay > az then
          if hitPos.y > aabb.center.y then Vec3.unitY else Vec3.down
        else
          if hitPos.z > aabb.center.z then Vec3.unitZ else Vec3.forward

    some { t, point := closest, normal }

/-- Swept AABB vs static AABB collision (Minkowski difference). -/
def aabbVsAABB (start end_ halfExtents : Vec3)
    (staticAABB : AABB) : Option SweptHit :=
  -- Minkowski difference: expand static AABB by moving AABB's half-extents
  let expandedMin := staticAABB.min.sub halfExtents
  let expandedMax := staticAABB.max.add halfExtents

  let movement := end_.sub start

  -- Ray from start through expanded AABB
  let invDx := if Float.abs movement.x > 0.0001 then 1.0 / movement.x else Float.infinity
  let invDy := if Float.abs movement.y > 0.0001 then 1.0 / movement.y else Float.infinity
  let invDz := if Float.abs movement.z > 0.0001 then 1.0 / movement.z else Float.infinity

  let tx1 := (expandedMin.x - start.x) * invDx
  let tx2 := (expandedMax.x - start.x) * invDx
  let ty1 := (expandedMin.y - start.y) * invDy
  let ty2 := (expandedMax.y - start.y) * invDy
  let tz1 := (expandedMin.z - start.z) * invDz
  let tz2 := (expandedMax.z - start.z) * invDz

  let tMinX := Float.min tx1 tx2
  let tMaxX := Float.max tx1 tx2
  let tMinY := Float.min ty1 ty2
  let tMaxY := Float.max ty1 ty2
  let tMinZ := Float.min tz1 tz2
  let tMaxZ := Float.max tz1 tz2

  let tMin := Float.max tMinX (Float.max tMinY tMinZ)
  let tMax := Float.min tMaxX (Float.min tMaxY tMaxZ)

  if tMax < 0.0 || tMin > tMax || tMin > 1.0 then none
  else
    let t := if tMin < 0.0 then 0.0 else tMin

    -- Determine which face was hit (the one with the latest tMin)
    let normal :=
      if tMinX >= tMinY && tMinX >= tMinZ then
        if movement.x > 0.0 then Vec3.left else Vec3.unitX
      else if tMinY >= tMinZ then
        if movement.y > 0.0 then Vec3.down else Vec3.unitY
      else
        if movement.z > 0.0 then Vec3.forward else Vec3.unitZ

    let hitPos := start.add (movement.scale t)
    some { t, point := hitPos, normal }

end SweptCollision

-- ============================================================================
-- Continuous Collision Detection (Time of Impact)
-- ============================================================================

namespace ContinuousCollision

/-- Parameters for conservative advancement algorithm. -/
structure CAParams where
  maxIterations : Nat := 32
  tolerance : Float := 0.0001
  deriving Repr, Inhabited

/-- Compute time of impact between two moving spheres.
    Uses analytical solution for sphere-sphere. -/
def sphereVsSphereTOI (sphereA : Sphere) (velocityA : Vec3)
    (sphereB : Sphere) (velocityB : Vec3)
    (maxTime : Float := 1.0) : Option Float :=
  -- Relative motion: treat A as moving, B as stationary
  let relVel := velocityA.sub velocityB
  let relPos := sphereA.center.sub sphereB.center
  let combinedRadius := sphereA.radius + sphereB.radius

  -- Quadratic: |relPos + t * relVel|^2 = combinedRadius^2
  let a := relVel.dot relVel
  let b := 2.0 * relPos.dot relVel
  let c := relPos.dot relPos - combinedRadius * combinedRadius

  -- Already overlapping
  if c < 0.0 then some 0.0
  -- Moving apart or parallel
  else if a < 0.0001 then none
  else
    let discriminant := b * b - 4.0 * a * c
    if discriminant < 0.0 then none
    else
      let t := (-b - Float.sqrt discriminant) / (2.0 * a)
      if t < 0.0 || t > maxTime then none
      else some t

/-- Compute time of impact between a moving sphere and static AABB.
    Uses ray-casting against expanded AABB. -/
def sphereVsAABBTOI (sphere : Sphere) (velocity : Vec3)
    (aabb : AABB) (maxTime : Float := 1.0) : Option Float :=
  let startPos := sphere.center
  let endPos := startPos.add (velocity.scale maxTime)

  match SweptCollision.sphereVsAABB startPos endPos sphere.radius aabb with
  | none => none
  | some hit => some (hit.t * maxTime)

/-- Conservative advancement for general shapes.
    distFunc: returns signed distance between shapes at time t
    velBound: upper bound on the rate of distance change -/
def conservativeAdvancement (distFunc : Float → Float)
    (velBound : Float) (params : CAParams := {})
    (maxTime : Float := 1.0) : Option Float :=
  if velBound <= 0.0 then none
  else Id.run do
    let mut t := 0.0

    for _ in [:params.maxIterations] do
      let dist := distFunc t
      if dist < params.tolerance then
        return some t
      if dist < 0.0 then
        return some t  -- Penetrating
      -- Advance conservatively
      let dt := dist / velBound
      t := t + dt
      if t > maxTime then
        return none

    return none

/-- Binary search refinement for time of impact.
    Assumes distFunc(t0) > 0 and distFunc(t1) <= 0 (or close to 0). -/
def refineTOI (distFunc : Float → Float) (t0 t1 : Float)
    (iterations : Nat := 8) : Float := Id.run do
  let mut lo := t0
  let mut hi := t1

  for _ in [:iterations] do
    let mid := (lo + hi) * 0.5
    if distFunc mid > 0.0 then
      lo := mid
    else
      hi := mid

  return (lo + hi) * 0.5

/-- Compute time of impact between two moving AABBs. -/
def aabbVsAABBTOI (movingAABB : AABB) (velocity : Vec3)
    (staticAABB : AABB) (maxTime : Float := 1.0) : Option Float :=
  let startPos := movingAABB.center
  let endPos := startPos.add (velocity.scale maxTime)
  let halfExtents := movingAABB.extents

  match SweptCollision.aabbVsAABB startPos endPos halfExtents staticAABB with
  | none => none
  | some hit => some (hit.t * maxTime)

end ContinuousCollision

end Linalg
