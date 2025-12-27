/-
  Oriented Bounding Box primitive.

  An OBB is a box that can be rotated to any orientation, defined by:
  - center: the center point
  - halfExtents: half-size along each local axis
  - orientation: quaternion rotation from local to world space
-/

import Linalg.Vec3
import Linalg.Quat
import Linalg.Geometry.AABB

namespace Linalg

/-- Oriented Bounding Box defined by center, half-extents, and orientation. -/
structure OBB where
  center : Vec3
  halfExtents : Vec3
  orientation : Quat := Quat.identity
deriving Repr, BEq, Inhabited

namespace OBB

/-- Create an OBB from center and half-extents (axis-aligned). -/
def fromCenterExtents (center halfExtents : Vec3) : OBB :=
  { center := center, halfExtents := halfExtents, orientation := Quat.identity }

/-- Create an OBB from an AABB (no rotation). -/
def fromAABB (aabb : AABB) : OBB :=
  { center := aabb.center, halfExtents := aabb.extents, orientation := Quat.identity }

/-- Create an OBB from center, half-extents, and axis-angle rotation. -/
def fromAxisAngle (center halfExtents axis : Vec3) (angle : Float) : OBB :=
  { center := center, halfExtents := halfExtents, orientation := Quat.fromAxisAngle axis angle }

/-- Get the local X axis in world space. -/
def axisX (obb : OBB) : Vec3 :=
  obb.orientation.rotateVec3 Vec3.unitX

/-- Get the local Y axis in world space. -/
def axisY (obb : OBB) : Vec3 :=
  obb.orientation.rotateVec3 Vec3.unitY

/-- Get the local Z axis in world space. -/
def axisZ (obb : OBB) : Vec3 :=
  obb.orientation.rotateVec3 Vec3.unitZ

/-- Get all three local axes in world space. -/
def axes (obb : OBB) : Vec3 × Vec3 × Vec3 :=
  (obb.axisX, obb.axisY, obb.axisZ)

/-- Transform a world-space point to local OBB space. -/
def worldToLocal (obb : OBB) (p : Vec3) : Vec3 :=
  -- Translate to OBB origin, then rotate by inverse orientation
  let relative := p.sub obb.center
  obb.orientation.inverse.rotateVec3 relative

/-- Transform a local OBB-space point to world space. -/
def localToWorld (obb : OBB) (p : Vec3) : Vec3 :=
  -- Rotate by orientation, then translate
  (obb.orientation.rotateVec3 p).add obb.center

/-- Check if a point is contained within the OBB. -/
def containsPoint (obb : OBB) (p : Vec3) : Bool :=
  let localPt := obb.worldToLocal p
  Float.abs' localPt.x <= obb.halfExtents.x &&
  Float.abs' localPt.y <= obb.halfExtents.y &&
  Float.abs' localPt.z <= obb.halfExtents.z

/-- Get the closest point on the OBB surface to a given point. -/
def closestPoint (obb : OBB) (p : Vec3) : Vec3 :=
  -- Transform to local space
  let localPt := obb.worldToLocal p
  -- Clamp to box in local space
  let clamped := Vec3.mk
    (Float.clamp localPt.x (-obb.halfExtents.x) obb.halfExtents.x)
    (Float.clamp localPt.y (-obb.halfExtents.y) obb.halfExtents.y)
    (Float.clamp localPt.z (-obb.halfExtents.z) obb.halfExtents.z)
  -- Transform back to world space
  obb.localToWorld clamped

/-- Get the squared distance from a point to the OBB. -/
def distanceSquared (obb : OBB) (p : Vec3) : Float :=
  let closest := obb.closestPoint p
  p.distanceSquared closest

/-- Get the distance from a point to the OBB. -/
def distance (obb : OBB) (p : Vec3) : Float :=
  Float.sqrt (obb.distanceSquared p)

/-- Get all 8 corners of the OBB in world space. -/
def corners (obb : OBB) : Array Vec3 :=
  let hx := obb.halfExtents.x
  let hy := obb.halfExtents.y
  let hz := obb.halfExtents.z
  #[
    obb.localToWorld (Vec3.mk (-hx) (-hy) (-hz)),
    obb.localToWorld (Vec3.mk   hx  (-hy) (-hz)),
    obb.localToWorld (Vec3.mk (-hx)   hy  (-hz)),
    obb.localToWorld (Vec3.mk   hx    hy  (-hz)),
    obb.localToWorld (Vec3.mk (-hx) (-hy)   hz),
    obb.localToWorld (Vec3.mk   hx  (-hy)   hz),
    obb.localToWorld (Vec3.mk (-hx)   hy    hz),
    obb.localToWorld (Vec3.mk   hx    hy    hz)
  ]

/-- Get the AABB that encloses this OBB. -/
def toAABB (obb : OBB) : AABB := Id.run do
  let cs := obb.corners
  let mut minPt := cs.getD 0 Vec3.zero
  let mut maxPt := minPt
  for i in [1:8] do
    let c := cs.getD i Vec3.zero
    minPt := Vec3.min minPt c
    maxPt := Vec3.max maxPt c
  return AABB.fromMinMax minPt maxPt

/-- Get the volume of the OBB. -/
def volume (obb : OBB) : Float :=
  8.0 * obb.halfExtents.x * obb.halfExtents.y * obb.halfExtents.z

/-- Get the surface area of the OBB. -/
def surfaceArea (obb : OBB) : Float :=
  let hx := obb.halfExtents.x
  let hy := obb.halfExtents.y
  let hz := obb.halfExtents.z
  8.0 * (hx * hy + hy * hz + hz * hx)

/-- Rotate the OBB by an additional quaternion. -/
def rotate (obb : OBB) (q : Quat) : OBB :=
  { obb with
    center := q.rotateVec3 obb.center
    orientation := q.multiply obb.orientation }

/-- Translate the OBB by a vector. -/
def translate (obb : OBB) (v : Vec3) : OBB :=
  { obb with center := obb.center.add v }

/-- Scale the OBB uniformly. -/
def scale (obb : OBB) (s : Float) : OBB :=
  { obb with
    center := obb.center.scale s
    halfExtents := obb.halfExtents.scale s }

/-- Scale the OBB non-uniformly (in local space). -/
def scaleNonUniform (obb : OBB) (s : Vec3) : OBB :=
  { obb with
    halfExtents := Vec3.mk (obb.halfExtents.x * s.x) (obb.halfExtents.y * s.y) (obb.halfExtents.z * s.z) }

/-- Check if two OBBs approximately equal. -/
def approxEq (a b : OBB) (eps : Float := Float.epsilon) : Bool :=
  a.center.approxEq b.center eps &&
  a.halfExtents.approxEq b.halfExtents eps &&
  -- Quaternions can represent same rotation with opposite signs
  (a.orientation.dot b.orientation |> Float.abs' |> fun d => d > 1.0 - eps)

/-- Project the OBB onto an axis and get the half-width of the projection. -/
def projectOntoAxis (obb : OBB) (axis : Vec3) : Float :=
  let ax := obb.axisX
  let ay := obb.axisY
  let az := obb.axisZ
  obb.halfExtents.x * Float.abs' (ax.dot axis) +
  obb.halfExtents.y * Float.abs' (ay.dot axis) +
  obb.halfExtents.z * Float.abs' (az.dot axis)

/-- Test if two OBBs intersect using the Separating Axis Theorem. -/
def intersectsOBB (a b : OBB) : Bool := Id.run do
  -- Get axes of both OBBs
  let (ax1, ay1, az1) := a.axes
  let (ax2, ay2, az2) := b.axes

  -- Vector between centers
  let d := b.center.sub a.center

  -- Test the 6 face axes (3 from each OBB)
  let faceAxes := #[ax1, ay1, az1, ax2, ay2, az2]
  for axis in faceAxes do
    let ra := a.projectOntoAxis axis
    let rb := b.projectOntoAxis axis
    let dist := Float.abs' (d.dot axis)
    if dist > ra + rb then return false

  -- Test the 9 edge-edge cross product axes
  let aAxes := #[ax1, ay1, az1]
  let bAxes := #[ax2, ay2, az2]
  for aAxis in aAxes do
    for bAxis in bAxes do
      let cross := aAxis.cross bAxis
      -- Skip near-parallel axes
      if cross.lengthSquared > Float.epsilon then
        let axis := cross.normalize
        let ra := a.projectOntoAxis axis
        let rb := b.projectOntoAxis axis
        let dist := Float.abs' (d.dot axis)
        if dist > ra + rb then return false

  return true

/-- Test if the OBB intersects a sphere. -/
def intersectsSphere (obb : OBB) (center : Vec3) (radius : Float) : Bool :=
  obb.distanceSquared center <= radius * radius

/-- Test if the OBB intersects an AABB. -/
def intersectsAABB (obb : OBB) (aabb : AABB) : Bool :=
  obb.intersectsOBB (OBB.fromAABB aabb)

end OBB

end Linalg
