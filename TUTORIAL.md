# Linalg Tutorial

A hands-on guide to game math with Linalg, covering practical scenarios from basic vectors to physics simulation.

## Table of Contents

1. [Getting Started](#getting-started)
2. [Working with Vectors](#working-with-vectors)
3. [Transforms and Matrices](#transforms-and-matrices)
4. [Rotations with Quaternions](#rotations-with-quaternions)
5. [Collision Detection](#collision-detection)
6. [Animation and Easing](#animation-and-easing)
7. [Curves and Paths](#curves-and-paths)
8. [Procedural Generation with Noise](#procedural-generation-with-noise)
9. [Spatial Data Structures](#spatial-data-structures)
10. [Physics Simulation](#physics-simulation)

---

## Getting Started

Add Linalg to your `lakefile.lean`:

```lean
require linalg from git "https://github.com/nathanial/linalg" @ "v0.0.1"
```

Then import and open the namespace:

```lean
import Linalg

open Linalg
```

All examples assume these imports are in scope.

---

## Working with Vectors

### Creating Vectors

```lean
-- Explicit construction
let position := Vec3.mk 10.0 5.0 (-3.0)

-- Named constants
let origin := Vec3.zero          -- (0, 0, 0)
let direction := Vec3.forward    -- (0, 0, -1) right-handed
let up := Vec3.up                -- (0, 1, 0)

-- 2D vectors
let screenPos := Vec2.mk 400.0 300.0
```

### Basic Operations

```lean
-- Movement: position + velocity * dt
let velocity := Vec3.mk 5.0 0.0 2.0
let dt := 0.016  -- 60 FPS
let newPosition := position.add (velocity.scale dt)

-- Distance between two points
let target := Vec3.mk 20.0 5.0 10.0
let dist := position.distance target

-- Direction to target (normalized)
let toTarget := target.sub position |>.normalize
```

### Dot Product: Angles and Projections

The dot product tells you how aligned two vectors are:

```lean
-- Check if enemy is in front of player
let playerForward := Vec3.forward
let toEnemy := enemyPos.sub playerPos |>.normalize
let dot := playerForward.dot toEnemy

if dot > 0.0 then
  IO.println "Enemy is in front"
else
  IO.println "Enemy is behind"

-- Angle between vectors (in radians)
let angle := Float.acos (playerForward.dot toEnemy)

-- Project velocity onto a surface normal (for sliding collision)
let normal := Vec3.up
let projectedVelocity := velocity.project normal
let slideVelocity := velocity.sub projectedVelocity
```

### Cross Product: Perpendicular Vectors

The cross product gives a vector perpendicular to two input vectors:

```lean
-- Build a coordinate frame from a forward direction
let forward := Vec3.mk 1.0 0.0 1.0 |>.normalize
let worldUp := Vec3.up
let right := forward.cross worldUp |>.normalize
let up := right.cross forward

-- Calculate surface normal from triangle vertices
let v0 := Vec3.mk 0.0 0.0 0.0
let v1 := Vec3.mk 1.0 0.0 0.0
let v2 := Vec3.mk 0.0 1.0 0.0
let edge1 := v1.sub v0
let edge2 := v2.sub v0
let normal := edge1.cross edge2 |>.normalize
```

### Reflection

```lean
-- Reflect a velocity off a wall
let wallNormal := Vec3.mk (-1.0) 0.0 0.0 |>.normalize
let incomingVelocity := Vec3.mk 5.0 0.0 2.0
let reflectedVelocity := incomingVelocity.reflect wallNormal
```

### Linear Interpolation (Lerp)

Smoothly blend between two values:

```lean
-- Interpolate position for smooth movement
let startPos := Vec3.mk 0.0 0.0 0.0
let endPos := Vec3.mk 10.0 5.0 0.0
let t := 0.5  -- 0 = start, 1 = end
let midpoint := startPos.lerp endPos t

-- Interpolate colors
let red := Vec3.mk 1.0 0.0 0.0
let blue := Vec3.mk 0.0 0.0 1.0
let purple := red.lerp blue 0.5
```

---

## Transforms and Matrices

### The Transform Type

For most game objects, use the `Transform` type which combines position, rotation, and scale:

```lean
-- Create a transform
let transform := Transform.create
  (position := Vec3.mk 10.0 0.0 5.0)
  (rotation := Quat.fromAxisAngle Vec3.up Float.halfPi)
  (scale := Vec3.one)

-- Get local axes
let forward := transform.forward
let right := transform.right
let up := transform.up

-- Transform a point from local to world space
let localPoint := Vec3.mk 1.0 0.0 0.0
let worldPoint := transform.transformPoint localPoint

-- Transform a direction (ignores position)
let localDir := Vec3.forward
let worldDir := transform.transformDirection localDir
```

### Transform Hierarchy

Parent-child relationships for scene graphs:

```lean
-- Car body at world position
let carTransform := Transform.create
  (position := Vec3.mk 100.0 0.0 50.0)
  (rotation := Quat.fromAxisAngle Vec3.up 0.5)

-- Wheel positioned relative to car
let wheelLocalTransform := Transform.create
  (position := Vec3.mk 1.5 (-0.5) 2.0)  -- Front-right wheel offset
  (rotation := Quat.identity)

-- Compute wheel's world transform
let wheelWorldTransform := carTransform.compose wheelLocalTransform

-- Get wheel position in world space
let wheelWorldPos := wheelWorldTransform.position
```

### Matrix Operations

For low-level control or GPU upload, use matrices directly:

```lean
-- Build a model matrix manually
let translation := Mat4.translation 10.0 5.0 3.0
let rotation := Mat4.rotationY Float.halfPi
let scale := Mat4.scaling 2.0 2.0 2.0

-- Order matters! Scale -> Rotate -> Translate
let modelMatrix := translation * rotation * scale

-- Transform a point
let localPos := Vec3.mk 1.0 0.0 0.0
let worldPos := modelMatrix.transformPoint localPos

-- Transform a direction (no translation)
let localDir := Vec3.forward
let worldDir := modelMatrix.transformDirection localDir
```

### View and Projection Matrices

For camera setup:

```lean
-- Camera parameters
let cameraPos := Vec3.mk 0.0 10.0 20.0
let target := Vec3.zero
let up := Vec3.up

-- View matrix (world -> camera space)
let viewMatrix := Mat4.lookAt cameraPos target up

-- Perspective projection
let fovY := Float.toRadians 60.0
let aspect := 16.0 / 9.0
let nearPlane := 0.1
let farPlane := 1000.0
let projMatrix := Mat4.perspective fovY aspect nearPlane farPlane

-- Orthographic projection (for 2D or UI)
let orthoMatrix := Mat4.orthographic
  (left := 0.0) (right := 1920.0)
  (bottom := 1080.0) (top := 0.0)  -- Y-down for screen coords
  (near := -1.0) (far := 1.0)

-- Combined view-projection matrix
let vpMatrix := projMatrix * viewMatrix
```

### 2D Transforms

For 2D games, use `Affine2D` (a 2x3 matrix):

```lean
-- Create a 2D transform
let transform := Affine2D.trs
  (translation := Vec2.mk 400.0 300.0)
  (rotation := Float.toRadians 45.0)
  (scale := Vec2.mk 2.0 2.0)

-- Transform a point
let localPos := Vec2.mk 10.0 0.0
let worldPos := transform.transformPoint localPos

-- Compose transforms
let parent := Affine2D.translation 100.0 100.0
let child := Affine2D.rotation Float.halfPi
let combined := parent.compose child
```

---

## Rotations with Quaternions

Quaternions avoid gimbal lock and interpolate smoothly.

### Creating Rotations

```lean
-- No rotation
let identity := Quat.identity

-- Rotate around an axis
let rotateY90 := Quat.fromAxisAngle Vec3.up Float.halfPi

-- From Euler angles (pitch, yaw, roll)
let euler := Quat.fromEuler
  (pitch := Float.toRadians 15.0)   -- Look up/down
  (yaw := Float.toRadians 45.0)     -- Turn left/right
  (roll := 0.0)                     -- Tilt

-- Look at a target
let lookAt := Quat.lookAt
  (forward := target.sub position |>.normalize)
  (up := Vec3.up)
```

### Applying Rotations

```lean
-- Rotate a vector
let direction := Vec3.forward
let rotation := Quat.fromAxisAngle Vec3.up Float.halfPi
let rotatedDir := rotation * direction  -- Uses HMul instance

-- Combine rotations (order matters!)
let pitch := Quat.fromAxisAngle Vec3.right (Float.toRadians 30.0)
let yaw := Quat.fromAxisAngle Vec3.up (Float.toRadians 45.0)
let combined := yaw * pitch  -- Yaw first, then pitch
```

### Smooth Rotation with Slerp

Spherical linear interpolation for smooth rotation:

```lean
-- Smoothly rotate from current to target orientation
let currentRotation := Quat.identity
let targetRotation := Quat.fromAxisAngle Vec3.up Float.pi
let t := 0.1  -- Interpolation factor (0-1)

let smoothedRotation := Quat.slerp currentRotation targetRotation t
```

### Converting Quaternions

```lean
-- To rotation matrix (for shaders)
let rotMatrix := rotation.toMat4

-- To Euler angles (for UI display)
let (pitch, yaw, roll) := rotation.toEuler
```

---

## Collision Detection

### Ray Casting

```lean
-- Create a ray from camera through mouse position
let rayOrigin := cameraPos
let rayDirection := screenToWorld mousePos |>.normalize
let ray := Ray.mk' rayOrigin rayDirection

-- Cast against a sphere
let sphere := Sphere.mk' (center := Vec3.zero) (radius := 5.0)
match Intersection.raySphere ray sphere with
| some hit =>
  let hitPoint := ray.pointAt hit.t
  let hitNormal := hit.normal
  IO.println s!"Hit at {hitPoint}"
| none =>
  IO.println "Miss"

-- Cast against an AABB
let box := AABB.fromCenterExtents Vec3.zero (Vec3.mk 2.0 2.0 2.0)
match Intersection.rayAABB ray box with
| some (tMin, tMax) =>
  let entryPoint := ray.pointAt tMin
  let exitPoint := ray.pointAt tMax
  IO.println s!"Entered at t={tMin}, exited at t={tMax}"
| none =>
  IO.println "Miss"

-- Cast against a triangle
let triangle := Triangle.mk
  (Vec3.mk 0.0 0.0 0.0)
  (Vec3.mk 5.0 0.0 0.0)
  (Vec3.mk 2.5 5.0 0.0)
match Intersection.rayTriangle ray triangle with
| some hit => IO.println s!"Hit triangle at t={hit.t}"
| none => IO.println "Miss"
```

### Primitive vs Primitive

```lean
-- Sphere-sphere collision
let sphere1 := Sphere.mk' Vec3.zero 2.0
let sphere2 := Sphere.mk' (Vec3.mk 3.0 0.0 0.0) 2.0
if Intersection.sphereSphere sphere1 sphere2 then
  IO.println "Spheres collide!"

-- AABB-AABB collision
let box1 := AABB.fromCenterExtents Vec3.zero Vec3.one
let box2 := AABB.fromCenterExtents (Vec3.mk 1.5 0.0 0.0) Vec3.one
if Intersection.aabbAABB box1 box2 then
  IO.println "Boxes overlap!"

-- Sphere-AABB collision
if Intersection.sphereAABB sphere1 box1 then
  IO.println "Sphere touches box!"
```

### Frustum Culling

Efficiently cull objects outside the camera view:

```lean
-- Build frustum from view-projection matrix
let frustum := Frustum.fromViewProjection vpMatrix

-- Test if objects are visible
let objectSphere := Sphere.mk' objectPos 5.0
let objectBox := AABB.fromCenterExtents objectPos (Vec3.mk 2.0 2.0 2.0)

if Intersection.frustumSphere frustum objectSphere then
  -- Object might be visible, render it
  render object

if Intersection.frustumAABB frustum objectBox then
  -- Object might be visible
  render object
```

### 2D Collision

```lean
-- Circle-circle
let circle1 := Circle.mk Vec2.zero 50.0
let circle2 := Circle.mk (Vec2.mk 60.0 0.0) 30.0
if Collision2D.circleCircle circle1 circle2 then
  IO.println "Circles overlap!"

-- Point in polygon
let polygon := Polygon2D.fromVertices #[
  Vec2.mk 0.0 0.0,
  Vec2.mk 100.0 0.0,
  Vec2.mk 100.0 100.0,
  Vec2.mk 0.0 100.0
]
let point := Vec2.mk 50.0 50.0
if polygon.containsPoint point then
  IO.println "Point is inside polygon!"

-- SAT collision with MTV (Minimum Translation Vector)
let poly1 := Polygon2D.regularPolygon Vec2.zero 50.0 6  -- Hexagon
let poly2 := Polygon2D.rectangle (Vec2.mk 40.0 0.0) 60.0 60.0
match Collision2D.satPolygonPolygon poly1 poly2 with
| some mtv =>
  -- mtv is the smallest vector to separate the shapes
  let separatedPos := poly1.centroid.add mtv
  IO.println s!"Collision! Separate by {mtv}"
| none =>
  IO.println "No collision"
```

---

## Animation and Easing

### Easing Functions

Make animations feel natural:

```lean
-- t goes from 0 to 1 over the animation duration
let t := elapsedTime / duration

-- Different easing curves
let linear := t
let smooth := Easing.smoothstep t           -- Smooth start and end
let quadIn := Easing.quadIn t               -- Slow start
let quadOut := Easing.quadOut t             -- Slow end
let quadInOut := Easing.quadInOut t         -- Slow start and end
let elastic := Easing.elasticOut t          -- Bouncy overshoot
let bounce := Easing.bounceOut t            -- Bouncing ball effect
let back := Easing.backOut t                -- Slight overshoot

-- Apply easing to position
let startPos := Vec3.mk 0.0 0.0 0.0
let endPos := Vec3.mk 100.0 0.0 0.0
let easedT := Easing.cubicInOut t
let animatedPos := startPos.lerp endPos easedT
```

### Spring Animation

Physics-based springy motion:

```lean
-- Spring parameters
let spring := Spring.mk
  (stiffness := 100.0)   -- Higher = faster oscillation
  (damping := 10.0)      -- Higher = less bounce

-- Update spring each frame
let (newValue, newVelocity) := Spring.update
  spring
  currentValue
  currentVelocity
  targetValue
  dt
```

### Smooth Damp

Unity-style smooth following:

```lean
-- Smooth camera follow
let (newPos, newVel) := Vec3.smoothDamp
  currentCameraPos
  currentVelocity
  targetPos
  (smoothTime := 0.3)    -- Approximate time to reach target
  (maxSpeed := 100.0)    -- Speed limit
  dt
```

---

## Curves and Paths

### Bezier Curves

```lean
-- Quadratic Bezier (3 control points)
let p0 := Vec3.mk 0.0 0.0 0.0
let p1 := Vec3.mk 5.0 10.0 0.0   -- Control point (curve bends toward this)
let p2 := Vec3.mk 10.0 0.0 0.0
let quadCurve := Bezier.quadratic p0 p1 p2

-- Sample points along the curve
for t in [0.0, 0.25, 0.5, 0.75, 1.0] do
  let point := quadCurve.evaluate t
  IO.println s!"t={t}: {point}"

-- Cubic Bezier (4 control points) - most common
let c0 := Vec3.mk 0.0 0.0 0.0
let c1 := Vec3.mk 3.0 10.0 0.0
let c2 := Vec3.mk 7.0 10.0 0.0
let c3 := Vec3.mk 10.0 0.0 0.0
let cubicCurve := Bezier.cubic c0 c1 c2 c3

-- Get tangent (direction) at a point
let tangent := cubicCurve.tangent 0.5 |>.normalize
```

### Catmull-Rom Splines

Smooth curves through control points:

```lean
-- Define waypoints
let waypoints := #[
  Vec3.mk 0.0 0.0 0.0,
  Vec3.mk 10.0 5.0 0.0,
  Vec3.mk 20.0 0.0 5.0,
  Vec3.mk 30.0 8.0 0.0
]

let spline := CatmullRom.fromPoints waypoints

-- Sample the full path
for i in [0:100] do
  let t := i.toFloat / 99.0
  let point := spline.evaluate t
  -- Draw or use point
```

### Arc-Length Parameterization

Move at constant speed along a curve:

```lean
-- Build arc-length table
let arcLengthCurve := cubicCurve.buildArcLengthTable 100

-- Move at constant speed
let totalLength := arcLengthCurve.totalLength
let speed := 5.0  -- Units per second
let distanceTraveled := speed * elapsedTime
let t := arcLengthCurve.tAtDistance distanceTraveled
let position := cubicCurve.evaluate t
```

---

## Procedural Generation with Noise

### Perlin Noise

```lean
-- 1D noise (for terrain height along a line)
let height := Noise.perlin1D x

-- 2D noise (for heightmaps, textures)
let value := Noise.perlin2D x y

-- 3D noise (for volumetric effects, 3D textures)
let density := Noise.perlin3D x y z
```

### Fractal Noise (FBM)

Layer multiple octaves for natural-looking results:

```lean
-- Fractal Brownian Motion
let terrain := Noise.fbm2D x y
  (octaves := 6)
  (lacunarity := 2.0)      -- Frequency multiplier per octave
  (persistence := 0.5)     -- Amplitude multiplier per octave

-- Ridged multifractal (for mountains)
let mountains := Noise.ridged2D x y
  (octaves := 6)
  (lacunarity := 2.0)
  (gain := 0.5)

-- Turbulence (absolute value creates billowy effect)
let clouds := Noise.turbulence2D x y (octaves := 4)
```

### Domain Warping

Distort coordinates for organic patterns:

```lean
-- Warp the input coordinates with noise
let warpedX := x + Noise.perlin2D (x * 0.5) (y * 0.5) * 2.0
let warpedY := y + Noise.perlin2D (x * 0.5 + 5.3) (y * 0.5 + 1.7) * 2.0
let value := Noise.perlin2D warpedX warpedY
```

### Worley (Cellular) Noise

For cell-like patterns:

```lean
let (f1, f2) := Noise.worley2D x y
-- f1 = distance to nearest point
-- f2 = distance to second nearest
let cellPattern := f2 - f1  -- Creates cell edges
```

---

## Spatial Data Structures

### Quadtree (2D)

For efficient 2D spatial queries:

```lean
-- Create a quadtree covering the game world
let bounds := AABB2D.fromMinMax Vec2.zero (Vec2.mk 1000.0 1000.0)
let quadtree := Quadtree.empty bounds (maxDepth := 8) (maxItems := 4)

-- Insert objects with bounding boxes
let quadtree := quadtree.insert objectId objectBounds

-- Query objects in a region
let queryRegion := AABB2D.fromCenterExtents playerPos (Vec2.mk 100.0 100.0)
let nearbyObjects := quadtree.query queryRegion

-- Find k nearest neighbors
let nearest := quadtree.kNearest playerPos 10
```

### Octree (3D)

For 3D spatial partitioning:

```lean
-- Create octree
let bounds := AABB.fromMinMax Vec3.zero (Vec3.mk 1000.0 1000.0 1000.0)
let octree := Octree.empty bounds (maxDepth := 8) (maxItems := 8)

-- Insert and query same as quadtree
let octree := octree.insert objectId objectBounds
let visible := octree.queryFrustum frustum
```

### BVH (Bounding Volume Hierarchy)

Optimal for ray tracing:

```lean
-- Build BVH from triangles
let triangles : Array (Triangle × UInt32) := ...
let bvh := BVH.build triangles

-- Fast ray intersection
match bvh.raycast ray with
| some (hit, triangleId) =>
  IO.println s!"Hit triangle {triangleId} at t={hit.t}"
| none =>
  IO.println "Miss"
```

### KD-Tree

For nearest neighbor queries:

```lean
-- Build from points
let points : Array Vec3 := ...
let kdtree := KDTree3D.build points

-- Find nearest point
let (nearestIdx, nearestDist) := kdtree.nearest queryPoint

-- Find all points within radius
let nearby := kdtree.radiusQuery queryPoint 10.0
```

---

## Physics Simulation

### Particle Physics

Simple point-mass physics:

```lean
-- Create a particle
let particle := Particle.create
  (position := Vec3.mk 0.0 10.0 0.0)
  (velocity := Vec3.mk 5.0 0.0 0.0)
  (mass := 1.0)

-- Apply gravity
let gravity := Vec3.mk 0.0 (-9.81) 0.0
let particle := particle.applyForce (gravity.scale particle.mass)

-- Integrate (multiple methods available)
let particle := Integration.semiImplicitEuler particle dt
-- Or: Integration.verlet, Integration.rk4
```

### Collision Response

Impulse-based collision resolution:

```lean
-- Two particles colliding
let contact := Contact.mk
  (point := collisionPoint)
  (normal := collisionNormal)
  (penetration := overlapAmount)

let (p1, p2) := CollisionResponse.resolveParticleCollision
  particle1 particle2 contact (restitution := 0.8)
```

### Rigid Body Physics

Full 3D rigid body with rotation:

```lean
-- Create a rigid body (box)
let mass := 10.0
let halfExtents := Vec3.mk 1.0 0.5 2.0
let inertiaTensor := InertiaTensor.solidBox mass halfExtents

let body := RigidBody.create
  (position := Vec3.mk 0.0 5.0 0.0)
  (mass := mass)
  (inertiaTensor := inertiaTensor)

-- Apply force at a point (creates torque)
let forcePoint := body.position.add (Vec3.mk 1.0 0.0 0.0)
let force := Vec3.mk 0.0 0.0 100.0
let body := body.applyForceAtPoint force forcePoint

-- Integrate
let body := Integration.integrateRigidBody body dt
```

### Swept Collision

Detect collisions during motion:

```lean
-- Sphere moving through space
let startPos := Vec3.mk 0.0 10.0 0.0
let endPos := Vec3.mk 0.0 0.0 0.0
let radius := 1.0

-- Check if it hits a plane (ground)
let ground := Plane.xz  -- Y = 0 plane
match SweptCollision.sphereVsPlane startPos endPos radius ground with
| some hit =>
  let hitPos := startPos.lerp endPos hit.t
  IO.println s!"Hit ground at t={hit.t}, position={hitPos}"
| none =>
  IO.println "No collision"
```

### Time of Impact

Find exact collision time:

```lean
-- Two spheres moving toward each other
let sphere1 := Sphere.mk' (Vec3.mk 0.0 0.0 0.0) 1.0
let vel1 := Vec3.mk 5.0 0.0 0.0

let sphere2 := Sphere.mk' (Vec3.mk 10.0 0.0 0.0) 1.0
let vel2 := Vec3.mk (-5.0) 0.0 0.0

match ContinuousCollision.sphereVsSphereTOI
  sphere1 vel1 sphere2 vel2 (maxTime := 1.0) with
| some toi =>
  IO.println s!"Collision at t={toi}"
  let pos1AtCollision := sphere1.center.add (vel1.scale toi)
  let pos2AtCollision := sphere2.center.add (vel2.scale toi)
| none =>
  IO.println "No collision within time window"
```

---

## Next Steps

- Browse the [API documentation](README.md#api-overview) for complete function listings
- Check out the [test suite](LinalgTests/) for more usage examples
- See [ROADMAP.md](ROADMAP.md) for planned features

Happy game programming!
