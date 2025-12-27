/-
  Linalg - Game Math Library for Lean 4

  A comprehensive linear algebra library designed for game programming,
  including 2D and 3D vector math, matrices, quaternions, and geometric
  primitives with intersection tests.

  Features:
  - Vectors: Vec2, Vec3, Vec4 with full algebra operations
  - Matrices: Mat2, Mat3, Mat4 (column-major, GPU-compatible)
  - Quaternions: Full rotation support with slerp
  - Geometric primitives: Ray, AABB, Sphere, Plane
  - Intersection tests: Ray-primitive, primitive-primitive
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

-- Geometric primitives
import Linalg.Geometry.Ray
import Linalg.Geometry.AABB
import Linalg.Geometry.Sphere
import Linalg.Geometry.Plane
import Linalg.Geometry.Intersection
