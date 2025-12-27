/-
  Plane primitive using normal + distance representation.
-/

import Linalg.Vec3
import Linalg.Core

namespace Linalg

/-- A plane defined by normal vector and distance from origin.
    The plane equation is: normal.dot(point) = distance -/
structure Plane where
  normal : Vec3  -- Should be normalized
  distance : Float
deriving Repr, BEq, Inhabited

namespace Plane

/-- Create a plane from normal and distance, normalizing the normal. -/
def fromNormalDistance (normal : Vec3) (distance : Float) : Plane :=
  ⟨normal.normalize, distance⟩

/-- Create a plane from a normal and a point on the plane. -/
def fromNormalPoint (normal : Vec3) (point : Vec3) : Plane :=
  let n := normal.normalize
  ⟨n, n.dot point⟩

/-- Create a plane from three non-collinear points (CCW winding). -/
def fromPoints (a b c : Vec3) : Plane :=
  let normal := (b.sub a).cross (c.sub a) |>.normalize
  ⟨normal, normal.dot a⟩

/-- Standard planes. -/
def xy : Plane := ⟨Vec3.unitZ, 0.0⟩  -- Z = 0
def xz : Plane := ⟨Vec3.unitY, 0.0⟩  -- Y = 0
def yz : Plane := ⟨Vec3.unitX, 0.0⟩  -- X = 0

/-- Get the signed distance from a point to the plane.
    Positive = in front (same side as normal), negative = behind. -/
def signedDistance (p : Plane) (point : Vec3) : Float :=
  p.normal.dot point - p.distance

/-- Get the absolute distance from a point to the plane. -/
def distanceToPoint (p : Plane) (point : Vec3) : Float :=
  Float.abs' (p.signedDistance point)

/-- Project a point onto the plane. -/
def projectPoint (p : Plane) (point : Vec3) : Vec3 :=
  point.sub (p.normal.scale (p.signedDistance point))

/-- Check if a point is on the plane (within epsilon). -/
def containsPoint (p : Plane) (point : Vec3) (eps : Float := Float.epsilon) : Bool :=
  Float.abs' (p.signedDistance point) < eps

/-- Check if a point is in front of the plane (positive side of normal). -/
def isInFront (p : Plane) (point : Vec3) : Bool :=
  p.signedDistance point > 0.0

/-- Check if a point is behind the plane (negative side of normal). -/
def isBehind (p : Plane) (point : Vec3) : Bool :=
  p.signedDistance point < 0.0

/-- Flip the plane (reverse normal and negate distance). -/
def flip (p : Plane) : Plane :=
  ⟨p.normal.neg, -p.distance⟩

/-- Get a point on the plane (the point closest to origin). -/
def origin (p : Plane) : Vec3 :=
  p.normal.scale p.distance

/-- Reflect a point across the plane. -/
def reflectPoint (p : Plane) (point : Vec3) : Vec3 :=
  point.sub (p.normal.scale (2.0 * p.signedDistance point))

/-- Reflect a direction vector across the plane. -/
def reflectDirection (p : Plane) (dir : Vec3) : Vec3 :=
  dir.reflect p.normal

/-- Normalize the plane (ensure normal is unit length). -/
def normalized (p : Plane) : Plane :=
  let len := p.normal.length
  if len > Float.epsilon then ⟨p.normal.scale (1.0 / len), p.distance / len⟩
  else p

/-- Check if two planes are approximately equal. -/
def approxEq (a b : Plane) (eps : Float := Float.epsilon) : Bool :=
  a.normal.approxEq b.normal eps && Float.approxEq a.distance b.distance eps

/-- Check if two planes are approximately parallel. -/
def isParallel (a b : Plane) (eps : Float := Float.epsilon) : Bool :=
  Float.approxEq (Float.abs' (a.normal.dot b.normal)) 1.0 eps

end Plane

end Linalg
