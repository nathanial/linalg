/-
  Signed distance field (SDF) primitives and composition operators.
-/

import Linalg.Core
import Linalg.Vec2
import Linalg.Vec3

namespace Linalg
namespace SDF

private def max2 (a b : Float) : Float := Float.max a b
private def max3 (a b c : Float) : Float := Float.max a (Float.max b c)
private def min2 (a b : Float) : Float := Float.min a b

/-- Closest point on segment AB to point P (3D). -/
def closestPointOnSegment3 (a b p : Vec3) : Vec3 :=
  let ab := b - a
  let lenSq := ab.lengthSquared
  if lenSq < Float.epsilon then
    a
  else
    let t := Float.clamp ((p - a).dot ab / lenSq) 0.0 1.0
    a + ab.scale t

/-- Distance from point to segment AB (3D). -/
def segmentDistance3 (a b p : Vec3) : Float :=
  (p - closestPointOnSegment3 a b p).length

/-- Closest point on segment AB to point P (2D). -/
def closestPointOnSegment2 (a b p : Vec2) : Vec2 :=
  let ab := b - a
  let lenSq := ab.lengthSquared
  if lenSq < Float.epsilon then
    a
  else
    let t := Float.clamp ((p - a).dot ab / lenSq) 0.0 1.0
    a + ab.scale t

/-- Distance from point to segment AB (2D). -/
def segmentDistance2 (a b p : Vec2) : Float :=
  (p - closestPointOnSegment2 a b p).length

/-- Sphere SDF (3D). -/
def sphere (p center : Vec3) (radius : Float) : Float :=
  (p - center).length - radius

/-- Circle SDF (2D). -/
def circle (p center : Vec2) (radius : Float) : Float :=
  (p - center).length - radius

/-- Axis-aligned box SDF (3D). halfExtents is half-size. -/
def box (p center halfExtents : Vec3) : Float :=
  let q := Vec3.abs (p - center) - halfExtents
  let outside := Vec3.max q Vec3.zero
  let inside := min2 (max3 q.x q.y q.z) 0.0
  outside.length + inside

/-- Axis-aligned box SDF (2D). halfExtents is half-size. -/
def box2D (p center halfExtents : Vec2) : Float :=
  let q := Vec2.abs (p - center) - halfExtents
  let outside := Vec2.max q Vec2.zero
  let inside := min2 (max2 q.x q.y) 0.0
  outside.length + inside

/-- Rounded box SDF (3D). -/
def roundBox (p center halfExtents : Vec3) (radius : Float) : Float :=
  box p center halfExtents - radius

/-- Rounded box SDF (2D). -/
def roundBox2D (p center halfExtents : Vec2) (radius : Float) : Float :=
  box2D p center halfExtents - radius

/-- Capsule SDF (3D) defined by endpoints and radius. -/
def capsule (p a b : Vec3) (radius : Float) : Float :=
  segmentDistance3 a b p - radius

/-- Capsule SDF (2D) defined by endpoints and radius. -/
def capsule2D (p a b : Vec2) (radius : Float) : Float :=
  segmentDistance2 a b p - radius

/-- Plane SDF: dot(p, normal) + offset. Normal should be normalized. -/
def plane (p normal : Vec3) (offset : Float) : Float :=
  p.dot normal + offset

/-- Union of two SDFs. -/
def opUnion (d1 d2 : Float) : Float :=
  Float.min d1 d2

/-- Intersection of two SDFs. -/
def opIntersection (d1 d2 : Float) : Float :=
  Float.max d1 d2

/-- Subtract second SDF from first. -/
def opSubtract (d1 d2 : Float) : Float :=
  Float.max d1 (-d2)

private def smoothMin (a b k : Float) : Float :=
  if k <= 0.0 then
    Float.min a b
  else
    let h := Float.clamp (0.5 + 0.5 * (b - a) / k) 0.0 1.0
    Float.lerp b a h - k * h * (1.0 - h)

/-- Smooth union of two SDFs. -/
def smoothUnion (d1 d2 k : Float) : Float :=
  smoothMin d1 d2 k

/-- Smooth intersection of two SDFs. -/
def smoothIntersection (d1 d2 k : Float) : Float :=
  -smoothMin (-d1) (-d2) k

/-- Smooth subtraction of SDFs. -/
def smoothSubtract (d1 d2 k : Float) : Float :=
  smoothIntersection d1 (-d2) k

end SDF
end Linalg
