/-
  Triangle primitive with barycentric coordinate support.
-/

import Linalg.Vec3

namespace Linalg

/-- A triangle defined by three vertices. -/
structure Triangle where
  v0 : Vec3
  v1 : Vec3
  v2 : Vec3
  deriving Repr

namespace Triangle

/-- Create a triangle from three points. -/
@[inline]
def mk' (a b c : Vec3) : Triangle :=
  { v0 := a, v1 := b, v2 := c }

/-- Edge from v0 to v1. -/
@[inline]
def edge01 (t : Triangle) : Vec3 := t.v1.sub t.v0

/-- Edge from v0 to v2. -/
@[inline]
def edge02 (t : Triangle) : Vec3 := t.v2.sub t.v0

/-- Edge from v1 to v2. -/
@[inline]
def edge12 (t : Triangle) : Vec3 := t.v2.sub t.v1

/-- Compute the (unnormalized) normal of the triangle using right-hand rule. -/
@[inline]
def normal (t : Triangle) : Vec3 :=
  t.edge01.cross t.edge02

/-- Compute the unit normal of the triangle. -/
@[inline]
def unitNormal (t : Triangle) : Vec3 :=
  t.normal.normalize

/-- Compute the area of the triangle. -/
@[inline]
def area (t : Triangle) : Float :=
  t.normal.length * 0.5

/-- Compute the centroid (center of mass) of the triangle. -/
@[inline]
def centroid (t : Triangle) : Vec3 :=
  Vec3.mk
    ((t.v0.x + t.v1.x + t.v2.x) / 3.0)
    ((t.v0.y + t.v1.y + t.v2.y) / 3.0)
    ((t.v0.z + t.v1.z + t.v2.z) / 3.0)

/-- Barycentric coordinates (u, v, w) where point = u*v0 + v*v1 + w*v2. -/
structure BarycentricCoords where
  u : Float
  v : Float
  w : Float
  deriving Repr

namespace BarycentricCoords

/-- Check if barycentric coordinates represent a point inside the triangle. -/
@[inline]
def isInside (bc : BarycentricCoords) : Bool :=
  bc.u >= 0.0 && bc.v >= 0.0 && bc.w >= 0.0

/-- Check if barycentric coordinates are valid (sum to 1). -/
@[inline]
def isValid (bc : BarycentricCoords) (epsilon : Float := 0.0001) : Bool :=
  Float.abs (bc.u + bc.v + bc.w - 1.0) < epsilon

end BarycentricCoords

/-- Compute barycentric coordinates of a point with respect to this triangle.
    Uses the method from "Real-Time Collision Detection" by Christer Ericson. -/
def barycentric (t : Triangle) (p : Vec3) : BarycentricCoords :=
  let v0v1 := t.edge01
  let v0v2 := t.edge02
  let v0p := p.sub t.v0

  let d00 := v0v1.dot v0v1
  let d01 := v0v1.dot v0v2
  let d11 := v0v2.dot v0v2
  let d20 := v0p.dot v0v1
  let d21 := v0p.dot v0v2

  let denom := d00 * d11 - d01 * d01
  if denom.abs < 1e-10 then
    -- Degenerate triangle
    { u := 1.0, v := 0.0, w := 0.0 }
  else
    let invDenom := 1.0 / denom
    let v := (d11 * d20 - d01 * d21) * invDenom
    let w := (d00 * d21 - d01 * d20) * invDenom
    let u := 1.0 - v - w
    { u := u, v := v, w := w }

/-- Convert barycentric coordinates back to a 3D point. -/
@[inline]
def fromBarycentric (t : Triangle) (bc : BarycentricCoords) : Vec3 :=
  (t.v0.scale bc.u).add ((t.v1.scale bc.v).add (t.v2.scale bc.w))

/-- Check if a point lies inside the triangle (in 3D, projects onto triangle plane). -/
def containsPoint (t : Triangle) (p : Vec3) : Bool :=
  let bc := t.barycentric p
  bc.isInside

/-- Compute the closest point on the triangle to a given point. -/
def closestPoint (t : Triangle) (p : Vec3) : Vec3 := Id.run do
  -- Check if P is in vertex region outside V0
  let ab := t.edge01
  let ac := t.edge02
  let ap := p.sub t.v0

  let d1 := ab.dot ap
  let d2 := ac.dot ap
  if d1 <= 0.0 && d2 <= 0.0 then
    return t.v0  -- barycentric coordinates (1,0,0)

  -- Check if P is in vertex region outside V1
  let bp := p.sub t.v1
  let d3 := ab.dot bp
  let d4 := ac.dot bp
  if d3 >= 0.0 && d4 <= d3 then
    return t.v1  -- barycentric coordinates (0,1,0)

  -- Check if P is in edge region of AB
  let vc := d1 * d4 - d3 * d2
  if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 then
    let v := d1 / (d1 - d3)
    return t.v0.add (ab.scale v)  -- barycentric coordinates (1-v,v,0)

  -- Check if P is in vertex region outside V2
  let cp := p.sub t.v2
  let d5 := ab.dot cp
  let d6 := ac.dot cp
  if d6 >= 0.0 && d5 <= d6 then
    return t.v2  -- barycentric coordinates (0,0,1)

  -- Check if P is in edge region of AC
  let vb := d5 * d2 - d1 * d6
  if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 then
    let w := d2 / (d2 - d6)
    return t.v0.add (ac.scale w)  -- barycentric coordinates (1-w,0,w)

  -- Check if P is in edge region of BC
  let va := d3 * d6 - d5 * d4
  if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 then
    let w := (d4 - d3) / ((d4 - d3) + (d5 - d6))
    return t.v1.add ((t.v2.sub t.v1).scale w)  -- barycentric coordinates (0,1-w,w)

  -- P is inside the triangle
  let denom := va + vb + vc
  let v := vb / denom
  let w := vc / denom
  return t.v0.add ((ab.scale v).add (ac.scale w))

/-- Compute bounding box of the triangle. -/
def boundingBox (t : Triangle) : Vec3 × Vec3 :=
  let minX := Float.min t.v0.x (Float.min t.v1.x t.v2.x)
  let minY := Float.min t.v0.y (Float.min t.v1.y t.v2.y)
  let minZ := Float.min t.v0.z (Float.min t.v1.z t.v2.z)
  let maxX := Float.max t.v0.x (Float.max t.v1.x t.v2.x)
  let maxY := Float.max t.v0.y (Float.max t.v1.y t.v2.y)
  let maxZ := Float.max t.v0.z (Float.max t.v1.z t.v2.z)
  (Vec3.mk minX minY minZ, Vec3.mk maxX maxY maxZ)

end Triangle

end Linalg
