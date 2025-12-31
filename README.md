# Linalg

A Lean 4 library for linear algebra and game math, providing vectors, matrices, quaternions, and geometric primitives optimized for 2D and 3D game programming.

## Features

- **Vectors**: Vec2, Vec3, Vec4 with full algebra (add, sub, scale, dot, cross, normalize, lerp)
- **Matrices**: Mat2, Mat3, Mat4 in column-major order (GPU-compatible)
- **Quaternions**: Rotation representation with slerp interpolation
- **Geometry**: Ray, AABB, Sphere, Plane primitives
- **Intersections**: Ray-sphere, ray-AABB, ray-plane, sphere-sphere, AABB-AABB tests

## Installation

Add to your `lakefile.lean`:

```lean
require linalg from git "https://github.com/nathanial/linalg" @ "v0.0.1"
```

## Usage

```lean
import Linalg

open Linalg

-- Vectors
let v := Vec3.mk 1.0 2.0 3.0
let u := Vec3.mk 4.0 5.0 6.0
let dot := v.dot u
let cross := v.cross u
let normalized := v.normalize

-- Matrices
let translation := Mat4.translation 10.0 20.0 30.0
let rotation := Mat4.rotationY Float.halfPi
let transform := translation * rotation
let point := transform.transformPoint Vec3.zero

-- Quaternions
let q := Quat.fromAxisAngle Vec3.unitY Float.halfPi
let rotated := q * Vec3.unitX
let interpolated := Quat.slerp Quat.identity q 0.5

-- Geometry
let ray := Ray.mk' (Vec3.mk 0.0 0.0 (-10.0)) Vec3.unitZ
let sphere := Sphere.mk' Vec3.zero 1.0
match Intersection.raySphere ray sphere with
| some hit => IO.println s!"Hit at t={hit.t}"
| none => IO.println "Miss"
```

## API Overview

### Core (`Linalg.Core`)

Float utilities and constants:
- `Float.pi`, `Float.halfPi`, `Float.twoPi`
- `Float.epsilon`, `Float.deg2Rad`, `Float.rad2Deg`
- `Float.lerp`, `Float.clamp`, `Float.approxEq`

### Vec2

```lean
structure Vec2 where x : Float; y : Float

-- Constants
Vec2.zero, Vec2.one, Vec2.unitX, Vec2.unitY

-- Operations
add, sub, neg, scale, dot, length, lengthSquared
normalize, distance, lerp, perpendicular, cross
```

### Vec3

```lean
structure Vec3 where x : Float; y : Float; z : Float

-- Constants
Vec3.zero, Vec3.one, Vec3.unitX, Vec3.unitY, Vec3.unitZ

-- Operations
add, sub, neg, scale, dot, cross, length, lengthSquared
normalize, distance, lerp, reflect, project
```

### Vec4

```lean
structure Vec4 where x : Float; y : Float; z : Float; w : Float

-- Constructors
Vec4.fromPoint (v : Vec3)      -- w = 1
Vec4.fromDirection (v : Vec3)  -- w = 0

-- Operations
add, sub, neg, scale, dot, length, normalize
toVec3, toVec3Normalized (perspective divide)
```

### Mat2, Mat3, Mat4

Column-major matrices for GPU compatibility:

```lean
-- Common operations
identity, transpose, determinant, inverse
add, sub, mul (matrix × matrix, matrix × vector)

-- Mat4 transforms
Mat4.translation, Mat4.scaling
Mat4.rotationX, Mat4.rotationY, Mat4.rotationZ
Mat4.perspective, Mat4.orthographic, Mat4.lookAt
Mat4.transformPoint, Mat4.transformDirection
```

### Quat

```lean
structure Quat where x : Float; y : Float; z : Float; w : Float

-- Constructors
Quat.identity
Quat.fromAxisAngle (axis : Vec3) (angle : Float)
Quat.fromEuler (pitch yaw roll : Float)

-- Operations
mul, conjugate, inverse, normalize, length
slerp, toMat3, toMat4, rotateVec3
HMul Quat Vec3 Vec3  -- q * v rotates v by q
```

### Geometry

```lean
structure Ray where origin : Vec3; direction : Vec3
structure AABB where min : Vec3; max : Vec3
structure Sphere where center : Vec3; radius : Float
structure Plane where normal : Vec3; distance : Float

-- Ray
Ray.mk' (normalizes direction)
Ray.pointAt (t : Float)

-- AABB
AABB.fromMinMax, AABB.fromCenterExtents
AABB.center, AABB.extents, AABB.containsPoint, AABB.merge

-- Sphere
Sphere.mk' (ensures positive radius)
Sphere.containsPoint, Sphere.boundingBox

-- Plane
Plane.xy, Plane.xz, Plane.yz (standard planes)
Plane.signedDistance, Plane.distanceToPoint, Plane.projectPoint
```

### Intersection

```lean
-- Ray intersections (return Option with hit info)
Intersection.raySphere : Ray → Sphere → Option RayHit
Intersection.rayAABB : Ray → AABB → Option (Float × Float)
Intersection.rayPlane : Ray → Plane → Option RayHit

-- Primitive intersections (return Bool)
Intersection.sphereSphere : Sphere → Sphere → Bool
Intersection.aabbAABB : AABB → AABB → Bool
Intersection.sphereAABB : Sphere → AABB → Bool
```

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Matrix storage | Column-major | GPU compatibility (Metal/OpenGL) |
| Numeric type | Float | Practical game math |
| Failing operations | Option types | Safe matrix inverse, intersection tests |
| Normalization | Auto in constructors | Ray direction, Plane normal must be unit vectors |

## Building

```bash
lake build    # Build library
lake test     # Run 629 tests
```

## License

MIT License - see [LICENSE](LICENSE) for details.
