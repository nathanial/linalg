/-
  View frustum for culling operations.
-/

import Linalg.Vec3
import Linalg.Vec4
import Linalg.Mat4
import Linalg.Geometry.Plane
import Linalg.Geometry.AABB
import Linalg.Geometry.Sphere

namespace Linalg

/-- Frustum plane indices. -/
inductive FrustumPlane where
  | near
  | far
  | left
  | right
  | top
  | bottom
  deriving Repr, BEq

/-- A view frustum defined by 6 planes.
    Plane normals point inward (toward the visible region). -/
structure Frustum where
  near   : Plane
  far    : Plane
  left   : Plane
  right  : Plane
  top    : Plane
  bottom : Plane
  deriving Repr

namespace Frustum

/-- Get a plane by index. -/
def getPlane (f : Frustum) (p : FrustumPlane) : Plane :=
  match p with
  | .near   => f.near
  | .far    => f.far
  | .left   => f.left
  | .right  => f.right
  | .top    => f.top
  | .bottom => f.bottom

/-- Get all 6 planes as an array. -/
def planes (f : Frustum) : Array Plane :=
  #[f.near, f.far, f.left, f.right, f.top, f.bottom]

/-- Extract frustum planes from a view-projection matrix.
    Uses the Gribb/Hartmann method for plane extraction.
    The matrix should be viewMatrix * projectionMatrix. -/
def fromViewProjection (vp : Mat4) : Frustum :=
  -- Extract rows of the matrix
  let row0 := Vec4.mk (vp.get 0 0) (vp.get 0 1) (vp.get 0 2) (vp.get 0 3)
  let row1 := Vec4.mk (vp.get 1 0) (vp.get 1 1) (vp.get 1 2) (vp.get 1 3)
  let row2 := Vec4.mk (vp.get 2 0) (vp.get 2 1) (vp.get 2 2) (vp.get 2 3)
  let row3 := Vec4.mk (vp.get 3 0) (vp.get 3 1) (vp.get 3 2) (vp.get 3 3)

  -- Extract planes: plane = row3 ± rowN
  -- Left:   row3 + row0
  -- Right:  row3 - row0
  -- Bottom: row3 + row1
  -- Top:    row3 - row1
  -- Near:   row3 + row2
  -- Far:    row3 - row2

  let leftVec := row3.add row0
  let rightVec := row3.sub row0
  let bottomVec := row3.add row1
  let topVec := row3.sub row1
  let nearVec := row3.add row2
  let farVec := row3.sub row2

  -- Normalize planes
  let normalizePlane (v : Vec4) : Plane :=
    let normal := Vec3.mk v.x v.y v.z
    let len := normal.length
    if len < Float.epsilon then
      Plane.mk Vec3.unitZ 0.0
    else
      let invLen := 1.0 / len
      Plane.mk (normal.scale invLen) (v.w * invLen)

  {
    left   := normalizePlane leftVec
    right  := normalizePlane rightVec
    bottom := normalizePlane bottomVec
    top    := normalizePlane topVec
    near   := normalizePlane nearVec
    far    := normalizePlane farVec
  }

/-- Create a frustum from individual planes. -/
def mk' (near far left right top bottom : Plane) : Frustum :=
  { near := near, far := far, left := left, right := right, top := top, bottom := bottom }

/-- Result of a containment test. -/
inductive Containment where
  | outside    -- Completely outside the frustum
  | intersects -- Partially inside
  | inside     -- Completely inside
  deriving Repr, BEq

/-- Test if a point is inside the frustum. -/
def containsPoint (f : Frustum) (p : Vec3) : Bool :=
  f.near.signedDistance p >= 0.0 &&
  f.far.signedDistance p >= 0.0 &&
  f.left.signedDistance p >= 0.0 &&
  f.right.signedDistance p >= 0.0 &&
  f.top.signedDistance p >= 0.0 &&
  f.bottom.signedDistance p >= 0.0

/-- Test sphere against frustum. Returns containment status. -/
def testSphere (f : Frustum) (s : Sphere) : Containment := Id.run do
  let planes := f.planes
  let mut allInside := true

  for plane in planes do
    let dist := plane.signedDistance s.center
    if dist < -s.radius then
      return .outside
    if dist < s.radius then
      allInside := false

  if allInside then .inside else .intersects

/-- Quick sphere visibility test (is the sphere at least partially visible?). -/
def isSphereVisible (f : Frustum) (s : Sphere) : Bool :=
  f.testSphere s != .outside

/-- Test AABB against frustum. Returns containment status.
    Uses the p-vertex/n-vertex optimization. -/
def testAABB (f : Frustum) (aabb : AABB) : Containment := Id.run do
  let planes := f.planes
  let mut allInside := true

  for plane in planes do
    -- Compute positive vertex (p-vertex): corner furthest in direction of normal
    let px := if plane.normal.x >= 0.0 then aabb.max.x else aabb.min.x
    let py := if plane.normal.y >= 0.0 then aabb.max.y else aabb.min.y
    let pz := if plane.normal.z >= 0.0 then aabb.max.z else aabb.min.z
    let pVertex := Vec3.mk px py pz

    -- If p-vertex is outside, entire AABB is outside
    if plane.signedDistance pVertex < 0.0 then
      return .outside

    -- Compute negative vertex (n-vertex): corner closest in direction of normal
    let nx := if plane.normal.x >= 0.0 then aabb.min.x else aabb.max.x
    let ny := if plane.normal.y >= 0.0 then aabb.min.y else aabb.max.y
    let nz := if plane.normal.z >= 0.0 then aabb.min.z else aabb.max.z
    let nVertex := Vec3.mk nx ny nz

    -- If n-vertex is outside, AABB intersects the plane
    if plane.signedDistance nVertex < 0.0 then
      allInside := false

  if allInside then .inside else .intersects

/-- Quick AABB visibility test (is the AABB at least partially visible?). -/
def isAABBVisible (f : Frustum) (aabb : AABB) : Bool :=
  f.testAABB aabb != .outside

/-- Get the 8 corner points of the frustum (in world space).
    Requires the inverse of the view-projection matrix. -/
def corners (invViewProj : Mat4) : Array Vec3 :=
  -- NDC corners
  let ndcCorners := #[
    Vec4.mk (-1.0) (-1.0) (-1.0) 1.0,  -- near bottom left
    Vec4.mk   1.0  (-1.0) (-1.0) 1.0,  -- near bottom right
    Vec4.mk   1.0    1.0  (-1.0) 1.0,  -- near top right
    Vec4.mk (-1.0)   1.0  (-1.0) 1.0,  -- near top left
    Vec4.mk (-1.0) (-1.0)   1.0  1.0,  -- far bottom left
    Vec4.mk   1.0  (-1.0)   1.0  1.0,  -- far bottom right
    Vec4.mk   1.0    1.0    1.0  1.0,  -- far top right
    Vec4.mk (-1.0)   1.0    1.0  1.0   -- far top left
  ]

  ndcCorners.map fun ndc =>
    let clip : Vec4 := invViewProj * ndc
    clip.toVec3Normalized

end Frustum

end Linalg
