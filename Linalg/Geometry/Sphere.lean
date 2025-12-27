/-
  Sphere primitive for collision and intersection tests.
-/

import Linalg.Vec3
import Linalg.Core

namespace Linalg

/-- A sphere defined by center and radius. -/
structure Sphere where
  center : Vec3
  radius : Float
deriving Repr, BEq, Inhabited

namespace Sphere

/-- Create a sphere with non-negative radius. -/
def mk' (center : Vec3) (radius : Float) : Sphere :=
  ⟨center, Float.abs' radius⟩

/-- Unit sphere at origin. -/
def unit : Sphere := ⟨Vec3.zero, 1.0⟩

/-- Check if a point is contained within the sphere. -/
def containsPoint (s : Sphere) (p : Vec3) : Bool :=
  s.center.distanceSquared p <= s.radius * s.radius

/-- Check if a point is on the surface of the sphere. -/
def onSurface (s : Sphere) (p : Vec3) (eps : Float := Float.epsilon) : Bool :=
  Float.approxEq (s.center.distance p) s.radius eps

/-- Get the surface area of the sphere. -/
def surfaceArea (s : Sphere) : Float :=
  4.0 * Float.pi * s.radius * s.radius

/-- Get the volume of the sphere. -/
def volume (s : Sphere) : Float :=
  (4.0 / 3.0) * Float.pi * s.radius * s.radius * s.radius

/-- Get the diameter of the sphere. -/
def diameter (s : Sphere) : Float := s.radius * 2.0

/-- Get the closest point on the sphere surface to a given point. -/
def closestPoint (s : Sphere) (p : Vec3) : Vec3 :=
  let dir := (p.sub s.center).normalize
  s.center.add (dir.scale s.radius)

/-- Get the distance from a point to the sphere surface. -/
def distanceToSurface (s : Sphere) (p : Vec3) : Float :=
  Float.abs' (s.center.distance p - s.radius)

/-- Expand the sphere to contain a point. -/
def expand (s : Sphere) (p : Vec3) : Sphere :=
  let dist := s.center.distance p
  if dist <= s.radius then s
  else ⟨s.center, dist⟩

/-- Merge two spheres into a bounding sphere. -/
def merge (a b : Sphere) : Sphere :=
  let dist := a.center.distance b.center
  -- Check if one sphere contains the other
  if dist + b.radius <= a.radius then a
  else if dist + a.radius <= b.radius then b
  else
    -- Calculate new bounding sphere
    let newRadius := (dist + a.radius + b.radius) / 2.0
    let dir := (b.center.sub a.center).normalize
    let offset := newRadius - a.radius
    let newCenter := a.center.add (dir.scale offset)
    ⟨newCenter, newRadius⟩

/-- Check if two spheres are approximately equal. -/
def approxEq (a b : Sphere) (eps : Float := Float.epsilon) : Bool :=
  a.center.approxEq b.center eps && Float.approxEq a.radius b.radius eps

end Sphere

end Linalg
