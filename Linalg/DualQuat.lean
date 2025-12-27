/-
  Dual Quaternion type for rigid body transforms.
  Combines rotation and translation in a single representation.
  Useful for skeletal animation blending and screw motion interpolation.
-/

import Linalg.Core
import Linalg.Vec3
import Linalg.Quat
import Linalg.Mat4

namespace Linalg

/-- Dual quaternion representing a rigid body transform (rotation + translation).
    Stored as real part (rotation) and dual part (translation encoding).

    A dual quaternion is: q = q_r + ε * q_d
    where ε² = 0 (dual number property), q_r is rotation, q_d encodes translation.

    For a rigid transform with rotation r and translation t:
    - q_r = r (unit quaternion)
    - q_d = 0.5 * t_quat * r, where t_quat = (t.x, t.y, t.z, 0) -/
structure DualQuat where
  /-- Real part: rotation quaternion (should be unit length). -/
  real : Quat
  /-- Dual part: encodes translation relative to rotation. -/
  dual : Quat
  deriving Repr, BEq, Inhabited

namespace DualQuat

/-- Identity dual quaternion (no rotation, no translation). -/
def identity : DualQuat :=
  { real := Quat.identity
    dual := ⟨0, 0, 0, 0⟩ }

/-- Create dual quaternion from rotation quaternion only (no translation). -/
def fromRotation (rot : Quat) : DualQuat :=
  { real := rot.normalize
    dual := ⟨0, 0, 0, 0⟩ }

/-- Create dual quaternion from translation only (no rotation). -/
def fromTranslation (t : Vec3) : DualQuat :=
  { real := Quat.identity
    dual := ⟨t.x * 0.5, t.y * 0.5, t.z * 0.5, 0⟩ }

/-- Create dual quaternion from rotation and translation.
    This is the primary constructor for rigid transforms. -/
def fromRotationTranslation (rot : Quat) (t : Vec3) : DualQuat :=
  let r := rot.normalize
  -- q_d = 0.5 * t_quat * r, where t_quat = (t.x, t.y, t.z, 0)
  let t_quat : Quat := ⟨t.x, t.y, t.z, 0⟩
  let dual := Quat.multiply t_quat r
  { real := r
    dual := ⟨dual.x * 0.5, dual.y * 0.5, dual.z * 0.5, dual.w * 0.5⟩ }

/-- Create dual quaternion from 4x4 transformation matrix. -/
def fromMat4 (m : Mat4) : DualQuat :=
  let rot := Quat.fromMat4 m
  let t := m.getTranslation
  fromRotationTranslation rot t

/-- Get the rotation component as a quaternion. -/
def getRotation (dq : DualQuat) : Quat := dq.real

/-- Get the translation component as a vector. -/
def getTranslation (dq : DualQuat) : Vec3 :=
  -- t = 2 * q_d * q_r^(-1) = 2 * q_d * conjugate(q_r) for unit quaternions
  let conj := dq.real.conjugate
  let t_quat := Quat.multiply (⟨dq.dual.x * 2, dq.dual.y * 2, dq.dual.z * 2, dq.dual.w * 2⟩) conj
  ⟨t_quat.x, t_quat.y, t_quat.z⟩

/-- Convert to 4x4 transformation matrix. -/
def toMat4 (dq : DualQuat) : Mat4 :=
  let rotMat := dq.real.toMat4
  let t := dq.getTranslation
  -- Set translation in last column
  { data := #[
    rotMat.get 0 0, rotMat.get 1 0, rotMat.get 2 0, 0,
    rotMat.get 0 1, rotMat.get 1 1, rotMat.get 2 1, 0,
    rotMat.get 0 2, rotMat.get 1 2, rotMat.get 2 2, 0,
    t.x, t.y, t.z, 1
  ]}

/-- Squared magnitude of the dual quaternion. -/
def lengthSquared (dq : DualQuat) : Float := dq.real.lengthSquared

/-- Magnitude of the dual quaternion. -/
def length (dq : DualQuat) : Float := dq.real.length

/-- Normalize the dual quaternion (makes real part unit length). -/
def normalize (dq : DualQuat) : DualQuat :=
  let len := dq.real.length
  if len > Float.epsilon then
    let invLen := 1.0 / len
    { real := ⟨dq.real.x * invLen, dq.real.y * invLen, dq.real.z * invLen, dq.real.w * invLen⟩
      dual := ⟨dq.dual.x * invLen, dq.dual.y * invLen, dq.dual.z * invLen, dq.dual.w * invLen⟩ }
  else
    identity

/-- Check if dual quaternion is normalized. -/
def isNormalized (dq : DualQuat) (eps : Float := 0.0001) : Bool :=
  Float.approxEq dq.real.lengthSquared 1.0 eps

/-- Conjugate of a dual quaternion.
    For unit dual quaternion, this gives the inverse transform. -/
def conjugate (dq : DualQuat) : DualQuat :=
  { real := dq.real.conjugate
    dual := dq.dual.conjugate }

/-- Dual conjugate (negates dual part only). -/
def dualConjugate (dq : DualQuat) : DualQuat :=
  { real := dq.real
    dual := ⟨-dq.dual.x, -dq.dual.y, -dq.dual.z, -dq.dual.w⟩ }

/-- Full conjugate (both quaternion and dual conjugate). -/
def fullConjugate (dq : DualQuat) : DualQuat :=
  { real := dq.real.conjugate
    dual := ⟨-dq.dual.conjugate.x, -dq.dual.conjugate.y, -dq.dual.conjugate.z, -dq.dual.conjugate.w⟩ }

/-- Negate the dual quaternion. -/
def neg (dq : DualQuat) : DualQuat :=
  { real := dq.real.neg
    dual := dq.dual.neg }

/-- Add two dual quaternions (component-wise). -/
def add (a b : DualQuat) : DualQuat :=
  { real := ⟨a.real.x + b.real.x, a.real.y + b.real.y, a.real.z + b.real.z, a.real.w + b.real.w⟩
    dual := ⟨a.dual.x + b.dual.x, a.dual.y + b.dual.y, a.dual.z + b.dual.z, a.dual.w + b.dual.w⟩ }

/-- Scale a dual quaternion by a scalar. -/
def scale (dq : DualQuat) (s : Float) : DualQuat :=
  { real := ⟨dq.real.x * s, dq.real.y * s, dq.real.z * s, dq.real.w * s⟩
    dual := ⟨dq.dual.x * s, dq.dual.y * s, dq.dual.z * s, dq.dual.w * s⟩ }

/-- Multiply two dual quaternions (composition of transforms).
    Result represents applying b first, then a. -/
def multiply (a b : DualQuat) : DualQuat :=
  -- (a_r + ε*a_d) * (b_r + ε*b_d) = a_r*b_r + ε*(a_r*b_d + a_d*b_r)
  let real := Quat.multiply a.real b.real
  let dual1 := Quat.multiply a.real b.dual
  let dual2 := Quat.multiply a.dual b.real
  { real := real
    dual := ⟨dual1.x + dual2.x, dual1.y + dual2.y, dual1.z + dual2.z, dual1.w + dual2.w⟩ }

/-- Dot product of two dual quaternions. -/
def dot (a b : DualQuat) : Float :=
  a.real.dot b.real

/-- Transform a point by the dual quaternion. -/
def transformPoint (dq : DualQuat) (p : Vec3) : Vec3 :=
  -- First rotate, then translate
  let rotated := dq.real.rotateVec3 p
  let t := dq.getTranslation
  rotated.add t

/-- Transform a direction by the dual quaternion (rotation only, no translation). -/
def transformVector (dq : DualQuat) (v : Vec3) : Vec3 :=
  dq.real.rotateVec3 v

/-- Inverse of a dual quaternion.
    For unit dual quaternions, this equals the conjugate. -/
def inverse (dq : DualQuat) : DualQuat :=
  -- For unit dual quaternion, inverse = conjugate
  if dq.isNormalized then
    dq.conjugate
  else
    -- General inverse: more complex computation
    let lenSq := dq.real.lengthSquared
    if lenSq > Float.epsilon then
      let invLenSq := 1.0 / lenSq
      let realConj := dq.real.conjugate
      let dualConj := dq.dual.conjugate
      { real := ⟨realConj.x * invLenSq, realConj.y * invLenSq, realConj.z * invLenSq, realConj.w * invLenSq⟩
        dual := ⟨dualConj.x * invLenSq, dualConj.y * invLenSq, dualConj.z * invLenSq, dualConj.w * invLenSq⟩ }
    else
      identity

/-- Dual quaternion Linear Blending (DLB).
    Blends multiple dual quaternions with weights.
    Better than matrix blending for skinning as it preserves rigid transforms. -/
def blend (dqs : Array DualQuat) (weights : Array Float) : DualQuat :=
  if dqs.size == 0 || weights.size == 0 then
    identity
  else
    let n := min dqs.size weights.size
    -- Ensure consistent hemisphere (avoid flipping)
    let pivot := dqs[0]!
    let zero : Quat := ⟨0, 0, 0, 0⟩
    let (sumReal, sumDual) := (List.range n).foldl (init := (zero, zero)) fun (accReal, accDual) i =>
      match dqs[i]?, weights[i]? with
      | some dq, some w =>
        -- Check if on same hemisphere
        let d := if pivot.real.dot dq.real < 0 then dq.neg else dq
        let newReal : Quat := ⟨accReal.x + d.real.x * w, accReal.y + d.real.y * w,
                        accReal.z + d.real.z * w, accReal.w + d.real.w * w⟩
        let newDual : Quat := ⟨accDual.x + d.dual.x * w, accDual.y + d.dual.y * w,
                        accDual.z + d.dual.z * w, accDual.w + d.dual.w * w⟩
        (newReal, newDual)
      | _, _ => (accReal, accDual)
    (DualQuat.mk sumReal sumDual).normalize

/-- Linear interpolation between two dual quaternions (not normalized). -/
def lerpUnnormalized (a b : DualQuat) (t : Float) : DualQuat :=
  { real := Quat.lerpUnnormalized a.real b.real t
    dual := ⟨Float.lerp a.dual.x b.dual.x t, Float.lerp a.dual.y b.dual.y t,
             Float.lerp a.dual.z b.dual.z t, Float.lerp a.dual.w b.dual.w t⟩ }

/-- Linear interpolation between two dual quaternions (normalized).
    Uses shortest path by checking hemisphere. -/
def lerp (a b : DualQuat) (t : Float) : DualQuat :=
  let b' := if a.real.dot b.real < 0 then b.neg else b
  (lerpUnnormalized a b' t).normalize

/-- Screw Linear Interpolation (ScLERP).
    Interpolates along the screw axis, producing constant velocity motion.
    More accurate than DLB for interpolating between two transforms. -/
def sclerp (a b : DualQuat) (t : Float) : DualQuat :=
  -- Handle near-identical transforms
  if (a.real.dot b.real).abs > 0.9999 then
    lerp a b t
  else
    -- Compute relative transform: delta = a^(-1) * b
    let aInv := a.inverse
    let delta := multiply aInv b
    -- Ensure positive rotation
    let delta' := if delta.real.w < 0 then delta.neg else delta
    -- Extract screw parameters
    let realLen := delta'.real.length
    if realLen < Float.epsilon then
      lerp a b t
    else
      -- Get rotation angle
      let halfAngle := Float.acos (Float.clamp delta'.real.w (-1.0) 1.0)
      if halfAngle.abs < Float.epsilon then
        lerp a b t
      else
        -- Interpolate angle
        let halfAngleT := halfAngle * t
        let s := Float.sin halfAngleT
        let c := Float.cos halfAngleT
        let sinHalfAngle := Float.sin halfAngle
        let invSinHalfAngle := if sinHalfAngle.abs > Float.epsilon then 1.0 / sinHalfAngle else 0.0

        -- Scale the axis
        let axis := Vec3.mk delta'.real.x delta'.real.y delta'.real.z
        let axisScaled := axis.scale (s * invSinHalfAngle)

        -- Compute interpolated real part
        let realInterp : Quat := ⟨axisScaled.x, axisScaled.y, axisScaled.z, c⟩

        -- Interpolate dual part (simplified for common case)
        let dualInterp := ⟨delta'.dual.x * t, delta'.dual.y * t, delta'.dual.z * t, delta'.dual.w * t⟩

        -- Compose: result = a * delta^t
        multiply a { real := realInterp, dual := dualInterp }

/-- Check approximate equality. -/
def approxEq (a b : DualQuat) (eps : Float := 0.0001) : Bool :=
  (a.real.approxEq b.real eps && a.dual.approxEq b.dual eps) ||
  (a.real.approxEq b.real.neg eps && a.dual.approxEq b.dual.neg eps)

/-- Check if two dual quaternions represent the same transform. -/
def sameTransform (a b : DualQuat) (eps : Float := 0.0001) : Bool :=
  let ta := a.getTranslation
  let tb := b.getTranslation
  a.real.sameRotation b.real eps &&
  ta.distance tb < eps

instance : Neg DualQuat := ⟨neg⟩
instance : HAdd DualQuat DualQuat DualQuat := ⟨add⟩
instance : HMul DualQuat DualQuat DualQuat := ⟨multiply⟩
instance : HMul DualQuat Vec3 Vec3 := ⟨transformPoint⟩

end DualQuat

end Linalg
