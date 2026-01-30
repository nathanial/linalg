/-
  Quaternion type and operations for rotation representation.
-/

import Linalg.Core
import Linalg.Vec3
import Linalg.Mat3
import Linalg.Mat4

namespace Linalg

/-- Quaternion for rotation representation (w + xi + yj + zk). -/
structure Quat where
  x : Float
  y : Float
  z : Float
  w : Float
deriving Repr, BEq, Inhabited

namespace Quat

/-- Identity quaternion (no rotation). -/
def identity : Quat := ⟨0.0, 0.0, 0.0, 1.0⟩

/-- Create quaternion from axis and angle (angle in radians). -/
def fromAxisAngle (axis : Vec3) (angle : Float) : Quat :=
  let halfAngle := angle / 2.0
  let s := Float.sin halfAngle
  let n := axis.normalize
  ⟨n.x * s, n.y * s, n.z * s, Float.cos halfAngle⟩

/-- Create quaternion from Euler angles (pitch, yaw, roll in radians, XYZ order). -/
def fromEuler (pitch yaw roll : Float) : Quat :=
  let cy := Float.cos (yaw / 2.0)
  let sy := Float.sin (yaw / 2.0)
  let cp := Float.cos (pitch / 2.0)
  let sp := Float.sin (pitch / 2.0)
  let cr := Float.cos (roll / 2.0)
  let sr := Float.sin (roll / 2.0)
  ⟨sr * cp * cy - cr * sp * sy,
   cr * sp * cy + sr * cp * sy,
   cr * cp * sy - sr * sp * cy,
   cr * cp * cy + sr * sp * sy⟩

/-- Create quaternion from rotation matrix. -/
def fromMat3 (m : Mat3) : Quat :=
  let trace := m.get 0 0 + m.get 1 1 + m.get 2 2
  if trace > 0.0 then
    let s := Float.sqrt (trace + 1.0) * 2.0
    ⟨(m.get 2 1 - m.get 1 2) / s,
     (m.get 0 2 - m.get 2 0) / s,
     (m.get 1 0 - m.get 0 1) / s,
     s / 4.0⟩
  else if m.get 0 0 > m.get 1 1 && m.get 0 0 > m.get 2 2 then
    let s := Float.sqrt (1.0 + m.get 0 0 - m.get 1 1 - m.get 2 2) * 2.0
    ⟨s / 4.0,
     (m.get 0 1 + m.get 1 0) / s,
     (m.get 0 2 + m.get 2 0) / s,
     (m.get 2 1 - m.get 1 2) / s⟩
  else if m.get 1 1 > m.get 2 2 then
    let s := Float.sqrt (1.0 + m.get 1 1 - m.get 0 0 - m.get 2 2) * 2.0
    ⟨(m.get 0 1 + m.get 1 0) / s,
     s / 4.0,
     (m.get 1 2 + m.get 2 1) / s,
     (m.get 0 2 - m.get 2 0) / s⟩
  else
    let s := Float.sqrt (1.0 + m.get 2 2 - m.get 0 0 - m.get 1 1) * 2.0
    ⟨(m.get 0 2 + m.get 2 0) / s,
     (m.get 1 2 + m.get 2 1) / s,
     s / 4.0,
     (m.get 1 0 - m.get 0 1) / s⟩

/-- Create quaternion from 4x4 rotation matrix (uses upper-left 3x3). -/
def fromMat4 (m : Mat4) : Quat := fromMat3 m.toMat3

/-- Squared length of the quaternion. -/
@[inline]
def lengthSquared (q : Quat) : Float := q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w

/-- Length of the quaternion. -/
@[inline]
def length (q : Quat) : Float := Float.sqrt q.lengthSquared

/-- Normalize the quaternion to unit length. -/
def normalize (q : Quat) : Quat :=
  let len := q.length
  if len > Float.epsilon then ⟨q.x / len, q.y / len, q.z / len, q.w / len⟩
  else identity

/-- Check if quaternion is approximately unit length. -/
def isNormalized (q : Quat) (eps : Float := Float.epsilon) : Bool :=
  Float.approxEq q.lengthSquared 1.0 eps

/-- Conjugate of the quaternion. -/
def conjugate (q : Quat) : Quat := ⟨-q.x, -q.y, -q.z, q.w⟩

/-- Inverse of the quaternion. -/
def inverse (q : Quat) : Quat :=
  let lenSq := q.lengthSquared
  if lenSq > Float.epsilon then
    let invLenSq := 1.0 / lenSq
    ⟨-q.x * invLenSq, -q.y * invLenSq, -q.z * invLenSq, q.w * invLenSq⟩
  else identity

/-- Negate the quaternion. -/
def neg (q : Quat) : Quat := ⟨-q.x, -q.y, -q.z, -q.w⟩

/-- Quaternion multiplication (Hamilton product). -/
def multiply (a b : Quat) : Quat :=
  ⟨a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
   a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
   a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
   a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z⟩

/-- Dot product of two quaternions. -/
@[inline]
def dot (a b : Quat) : Float := a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w

/-- Rotate a vector by this quaternion. Uses the formula: q * v * q^(-1) -/
def rotateVec3 (q : Quat) (v : Vec3) : Vec3 :=
  -- Optimized rotation: v' = v + 2w(qv × v) + 2(qv × (qv × v))
  let qv := Vec3.mk q.x q.y q.z
  let uv := qv.cross v
  let uuv := qv.cross uv
  v + (uv.scale (2.0 * q.w)) + (uuv.scale 2.0)

/-- Linear interpolation (not normalized). -/
def lerpUnnormalized (a b : Quat) (t : Float) : Quat :=
  ⟨Float.lerp a.x b.x t, Float.lerp a.y b.y t,
   Float.lerp a.z b.z t, Float.lerp a.w b.w t⟩

/-- Linear interpolation (normalized). -/
def lerp (a b : Quat) (t : Float) : Quat :=
  (lerpUnnormalized a b t).normalize

/-- Spherical linear interpolation. -/
def slerp (a b : Quat) (t : Float) : Quat :=
  let dotProd := a.dot b
  -- Handle negative dot (shortest path)
  let (b', dotProd') := if dotProd < 0.0 then (b.neg, -dotProd) else (b, dotProd)
  -- Use linear interpolation for very close quaternions
  if dotProd' > 0.9995 then
    lerp a b' t
  else
    let theta := Float.acos (Float.clamp dotProd' (-1.0) 1.0)
    let sinTheta := Float.sin theta
    if Float.abs' sinTheta < Float.epsilon then
      lerp a b' t
    else
      let wa := Float.sin ((1.0 - t) * theta) / sinTheta
      let wb := Float.sin (t * theta) / sinTheta
      ⟨a.x * wa + b'.x * wb, a.y * wa + b'.y * wb,
       a.z * wa + b'.z * wb, a.w * wa + b'.w * wb⟩

/-- Convert to 3x3 rotation matrix. -/
def toMat3 (q : Quat) : Mat3 :=
  let x2 := q.x * 2.0; let y2 := q.y * 2.0; let z2 := q.z * 2.0
  let xx := q.x * x2; let xy := q.x * y2; let xz := q.x * z2
  let yy := q.y * y2; let yz := q.y * z2; let zz := q.z * z2
  let wx := q.w * x2; let wy := q.w * y2; let wz := q.w * z2
  Mat3.fromColumns
    ⟨1.0 - (yy + zz), xy + wz, xz - wy⟩
    ⟨xy - wz, 1.0 - (xx + zz), yz + wx⟩
    ⟨xz + wy, yz - wx, 1.0 - (xx + yy)⟩

/-- Convert to 4x4 rotation matrix. -/
def toMat4 (q : Quat) : Mat4 :=
  let x2 := q.x * 2.0; let y2 := q.y * 2.0; let z2 := q.z * 2.0
  let xx := q.x * x2; let xy := q.x * y2; let xz := q.x * z2
  let yy := q.y * y2; let yz := q.y * z2; let zz := q.z * z2
  let wx := q.w * x2; let wy := q.w * y2; let wz := q.w * z2
  { data := #[
    1.0 - (yy + zz), xy + wz, xz - wy, 0,
    xy - wz, 1.0 - (xx + zz), yz + wx, 0,
    xz + wy, yz - wx, 1.0 - (xx + yy), 0,
    0, 0, 0, 1
  ]}

/-- Get the rotation axis (normalized). -/
def getAxis (q : Quat) : Vec3 :=
  let sinHalfAngleSq := 1.0 - q.w * q.w
  if sinHalfAngleSq < Float.epsilon then
    Vec3.unitY  -- Default axis for identity quaternion
  else
    let invSinHalfAngle := 1.0 / Float.sqrt sinHalfAngleSq
    ⟨q.x * invSinHalfAngle, q.y * invSinHalfAngle, q.z * invSinHalfAngle⟩

/-- Get the rotation angle in radians. -/
def getAngle (q : Quat) : Float :=
  2.0 * Float.acos (Float.clamp q.w (-1.0) 1.0)

/-- Get axis and angle as a tuple. -/
def toAxisAngle (q : Quat) : Vec3 × Float :=
  (q.getAxis, q.getAngle)

/-- Convert to Euler angles (pitch, yaw, roll) in radians. -/
def toEuler (q : Quat) : Vec3 :=
  -- Roll (x-axis rotation)
  let sinr_cosp := 2.0 * (q.w * q.x + q.y * q.z)
  let cosr_cosp := 1.0 - 2.0 * (q.x * q.x + q.y * q.y)
  let roll := Float.atan2 sinr_cosp cosr_cosp

  -- Pitch (y-axis rotation)
  let sinp := 2.0 * (q.w * q.y - q.z * q.x)
  let pitch := if Float.abs' sinp >= 1.0 then
    Float.halfPi * (if sinp >= 0.0 then 1.0 else -1.0)  -- Use 90 degrees if out of range
  else
    Float.asin sinp

  -- Yaw (z-axis rotation)
  let siny_cosp := 2.0 * (q.w * q.z + q.x * q.y)
  let cosy_cosp := 1.0 - 2.0 * (q.y * q.y + q.z * q.z)
  let yaw := Float.atan2 siny_cosp cosy_cosp

  ⟨pitch, yaw, roll⟩

/-- Check if two quaternions are approximately equal. -/
def approxEq (a b : Quat) (eps : Float := Float.epsilon) : Bool :=
  Float.approxEq a.x b.x eps && Float.approxEq a.y b.y eps &&
  Float.approxEq a.z b.z eps && Float.approxEq a.w b.w eps

/-- Check if two quaternions represent the same rotation (accounts for q = -q). -/
def sameRotation (a b : Quat) (eps : Float := Float.epsilon) : Bool :=
  a.approxEq b eps || a.approxEq b.neg eps

/-- Create a quaternion that looks in the given direction.
    forward: The direction to look at (will be normalized)
    worldUp: The world up vector (default Y-up) -/
def lookRotation (forward : Vec3) (worldUp : Vec3 := Vec3.unitY) : Quat :=
  let fwd := forward.normalize
  -- Handle degenerate case where forward is parallel to up
  let dot := fwd.dot worldUp
  let up := if Float.abs' dot > 0.999 then
    -- Forward is parallel to up, use a different up vector
    if Float.abs' fwd.y > 0.999 then Vec3.unitZ else Vec3.unitY
  else
    worldUp

  let right := up.cross fwd |>.normalize
  let realUp := fwd.cross right

  -- Build rotation matrix and convert to quaternion
  let m := Mat3.fromColumns right realUp fwd
  fromMat3 m

/-- Create a quaternion that rotates from one direction to another. -/
def fromToRotation (from_ to : Vec3) : Quat :=
  let fromN := from_.normalize
  let toN := to.normalize
  let dot := fromN.dot toN

  if dot > 0.999999 then
    -- Vectors are nearly parallel
    identity
  else if dot < -0.999999 then
    -- Vectors are nearly opposite, find a perpendicular axis
    let axis := if Float.abs' fromN.x < 0.9 then
      Vec3.unitX.cross fromN |>.normalize
    else
      Vec3.unitY.cross fromN |>.normalize
    fromAxisAngle axis Float.pi
  else
    let axis := fromN.cross toN
    let s := Float.sqrt ((1.0 + dot) * 2.0)
    let invS := 1.0 / s
    ⟨axis.x * invS, axis.y * invS, axis.z * invS, s * 0.5⟩

/-- Quaternion natural logarithm. Returns a pure quaternion (w=0).
    For unit quaternions: log(q) = (0, θ * axis) where θ is half-angle. -/
def log (q : Quat) : Quat :=
  let len := Float.sqrt (q.x * q.x + q.y * q.y + q.z * q.z)
  if len < Float.epsilon then
    ⟨0.0, 0.0, 0.0, 0.0⟩
  else
    let theta := Float.atan2 len q.w
    let scale := theta / len
    ⟨q.x * scale, q.y * scale, q.z * scale, 0.0⟩

/-- Quaternion exponential. Takes a pure quaternion (w=0) and returns unit quaternion. -/
def exp (q : Quat) : Quat :=
  let theta := Float.sqrt (q.x * q.x + q.y * q.y + q.z * q.z)
  if theta < Float.epsilon then
    ⟨0.0, 0.0, 0.0, 1.0⟩
  else
    let sinTheta := Float.sin theta
    let cosTheta := Float.cos theta
    let scale := sinTheta / theta
    ⟨q.x * scale, q.y * scale, q.z * scale, cosTheta⟩

/-- Compute SQUAD intermediate control point for smooth spline interpolation.
    Given three consecutive quaternions (prev, curr, next), computes the
    tangent control point at curr for use in squad interpolation. -/
def squadIntermediate (prev curr next : Quat) : Quat :=
  -- Ensure shortest path
  let prev' := if prev.dot curr < 0.0 then prev.neg else prev
  let next' := if next.dot curr < 0.0 then next.neg else next
  -- s_i = q_i * exp(-0.25 * (log(q_i^-1 * q_{i-1}) + log(q_i^-1 * q_{i+1})))
  let currInv := curr.inverse
  let logPrev := (currInv.multiply prev').log
  let logNext := (currInv.multiply next').log
  let sum : Quat := ⟨-0.25 * (logPrev.x + logNext.x),
                     -0.25 * (logPrev.y + logNext.y),
                     -0.25 * (logPrev.z + logNext.z),
                     0.0⟩
  (curr.multiply sum.exp).normalize

/-- SQUAD (Spherical Quadrangle Interpolation) for smooth quaternion spline.
    Interpolates between q1 and q2 using control points s1 and s2.
    Use `squadIntermediate` to compute s1 and s2 from neighboring quaternions.
    Parameter t should be in [0, 1]. -/
def squad (q1 q2 s1 s2 : Quat) (t : Float) : Quat :=
  let slerpQ := slerp q1 q2 t
  let slerpS := slerp s1 s2 t
  slerp slerpQ slerpS (2.0 * t * (1.0 - t))

/-- Interpolate through a sequence of quaternion keyframes using SQUAD.
    Takes an array of quaternions and a global parameter t in [0, n-1]
    where n is the number of quaternions. Returns interpolated rotation. -/
def squadPath (keyframes : Array Quat) (t : Float) : Quat :=
  if keyframes.size < 2 then
    keyframes.getD 0 identity
  else if keyframes.size == 2 then
    slerp (keyframes.getD 0 identity) (keyframes.getD 1 identity) t
  else
    let n := keyframes.size
    let tClamped := Float.clamp t 0.0 (Float.ofNat (n - 1) - Float.epsilon)
    let i := tClamped.toUInt32.toNat
    let localT := tClamped - Float.ofNat i

    let q0 := keyframes.getD (if i == 0 then 0 else i - 1) identity
    let q1 := keyframes.getD i identity
    let q2 := keyframes.getD (min (i + 1) (n - 1)) identity
    let q3 := keyframes.getD (min (i + 2) (n - 1)) identity

    let s1 := squadIntermediate q0 q1 q2
    let s2 := squadIntermediate q1 q2 q3
    squad q1 q2 s1 s2 localT

instance : Neg Quat := ⟨neg⟩
instance : HMul Quat Quat Quat := ⟨multiply⟩
instance : HMul Quat Vec3 Vec3 := ⟨rotateVec3⟩

end Quat

/-- Coerce Quat to Mat3 rotation matrix. -/
instance : Coe Quat Mat3 := ⟨Quat.toMat3⟩

/-- Coerce Quat to Mat4 rotation matrix. -/
instance : Coe Quat Mat4 := ⟨Quat.toMat4⟩

end Linalg
