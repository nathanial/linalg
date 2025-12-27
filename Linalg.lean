/-
  Linalg - Game Math Library for Lean 4

  A comprehensive linear algebra library designed for game programming,
  including 2D and 3D vector math, matrices, quaternions, and geometric
  primitives with intersection tests.

  Features:
  - Vectors: Vec2, Vec3, Vec4 with full algebra operations
  - Matrices: Mat2, Mat3, Mat4 (column-major, GPU-compatible)
  - Quaternions: Full rotation support with slerp
  - Euler angles: Multiple rotation orders (XYZ, YXZ, ZYX, etc.)
  - Dual quaternions: Rigid transforms for skeletal animation blending
  - Transforms: Combined position, rotation, scale with hierarchy support
  - 2D transforms: Affine2D (Mat2x3) and Rotation2D
  - Geometric primitives: Ray, AABB, Sphere, Plane, Triangle, Frustum, OBB, Capsule
  - Intersection tests: Ray-primitive, primitive-primitive, frustum culling
  - Easing functions: Animation and interpolation utilities
  - Curves: Bezier, Catmull-Rom, B-splines, Bezier patches
  - Noise: Perlin, Simplex, FBM, ridged, turbulence, domain warping
-/

-- Core utilities
import Linalg.Core

-- Vector types
import Linalg.Vec2
import Linalg.Vec3
import Linalg.Vec4

-- Matrix types
import Linalg.Mat2
import Linalg.Mat3
import Linalg.Mat4

-- Quaternions
import Linalg.Quat

-- Euler angles
import Linalg.Euler

-- Dual Quaternions
import Linalg.DualQuat

-- Transform
import Linalg.Transform

-- 2D Affine Transforms
import Linalg.Affine2D
import Linalg.Rotation2D

-- Animation
import Linalg.Easing

-- Curves and Splines
import Linalg.Curves

-- Procedural Noise
import Linalg.Noise

-- Geometric primitives
import Linalg.Geometry.Ray
import Linalg.Geometry.AABB
import Linalg.Geometry.Sphere
import Linalg.Geometry.Plane
import Linalg.Geometry.Triangle
import Linalg.Geometry.Frustum
import Linalg.Geometry.OBB
import Linalg.Geometry.Capsule
import Linalg.Geometry.Circle
import Linalg.Geometry.Line2D
import Linalg.Geometry.Polygon2D
import Linalg.Geometry.Intersection
