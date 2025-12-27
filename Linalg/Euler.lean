/-
  Euler angle utilities with different rotation orders.
  Supports both intrinsic (rotating axes) and extrinsic (fixed axes) conventions.
-/

import Linalg.Core
import Linalg.Vec3
import Linalg.Mat3
import Linalg.Quat

namespace Linalg

/-- Euler rotation order (intrinsic rotation sequence).
    For example, XYZ means rotate around X, then around the new Y, then around the new Z. -/
inductive EulerOrder
  | XYZ | XZY | YXZ | YZX | ZXY | ZYX
  deriving Repr, BEq, Inhabited

namespace EulerOrder

/-- Convert to string representation. -/
def toString : EulerOrder → String
  | XYZ => "XYZ"
  | XZY => "XZY"
  | YXZ => "YXZ"
  | YZX => "YZX"
  | ZXY => "ZXY"
  | ZYX => "ZYX"

instance : ToString EulerOrder := ⟨toString⟩

end EulerOrder

/-- Euler angles with explicit rotation order. -/
structure Euler where
  /-- Rotation angle around first axis (radians). -/
  a1 : Float
  /-- Rotation angle around second axis (radians). -/
  a2 : Float
  /-- Rotation angle around third axis (radians). -/
  a3 : Float
  /-- Rotation order. -/
  order : EulerOrder
  deriving Repr, BEq, Inhabited

namespace Euler

/-- Create Euler angles with XYZ order (pitch around X, yaw around Y, roll around Z). -/
def xyz (x y z : Float) : Euler := ⟨x, y, z, .XYZ⟩

/-- Create Euler angles with XZY order. -/
def xzy (x z y : Float) : Euler := ⟨x, z, y, .XZY⟩

/-- Create Euler angles with YXZ order (common in games: yaw, pitch, roll). -/
def yxz (y x z : Float) : Euler := ⟨y, x, z, .YXZ⟩

/-- Create Euler angles with YZX order. -/
def yzx (y z x : Float) : Euler := ⟨y, z, x, .YZX⟩

/-- Create Euler angles with ZXY order. -/
def zxy (z x y : Float) : Euler := ⟨z, x, y, .ZXY⟩

/-- Create Euler angles with ZYX order (aerospace/nautical: yaw, pitch, roll). -/
def zyx (z y x : Float) : Euler := ⟨z, y, x, .ZYX⟩

/-- Create zero Euler angles with given order. -/
def zero (order : EulerOrder := .XYZ) : Euler := ⟨0, 0, 0, order⟩

/-- Create Euler angles from degrees. -/
def fromDegrees (a1 a2 a3 : Float) (order : EulerOrder := .XYZ) : Euler :=
  ⟨a1 * Float.pi / 180.0, a2 * Float.pi / 180.0, a3 * Float.pi / 180.0, order⟩

/-- Get angles as degrees. -/
def toDegrees (e : Euler) : Vec3 :=
  let toDeg := 180.0 / Float.pi
  Vec3.mk (e.a1 * toDeg) (e.a2 * toDeg) (e.a3 * toDeg)

/-- Get angles as radians vector. -/
def toRadians (e : Euler) : Vec3 := Vec3.mk e.a1 e.a2 e.a3

/-- Rotation matrix around X axis. -/
private def rotX (angle : Float) : Mat3 :=
  let c := Float.cos angle
  let s := Float.sin angle
  Mat3.fromColumns
    (Vec3.mk 1 0 0)
    (Vec3.mk 0 c s)
    (Vec3.mk 0 (-s) c)

/-- Rotation matrix around Y axis. -/
private def rotY (angle : Float) : Mat3 :=
  let c := Float.cos angle
  let s := Float.sin angle
  Mat3.fromColumns
    (Vec3.mk c 0 (-s))
    (Vec3.mk 0 1 0)
    (Vec3.mk s 0 c)

/-- Rotation matrix around Z axis. -/
private def rotZ (angle : Float) : Mat3 :=
  let c := Float.cos angle
  let s := Float.sin angle
  Mat3.fromColumns
    (Vec3.mk c s 0)
    (Vec3.mk (-s) c 0)
    (Vec3.mk 0 0 1)

/-- Convert Euler angles to rotation matrix.
    Uses extrinsic rotation convention (fixed axes, multiply in order). -/
def toMat3 (e : Euler) : Mat3 :=
  -- Extrinsic rotations: multiply in order (first rotation is leftmost)
  -- E.g., for XYZ: result = Rx * Ry * Rz
  match e.order with
  | .XYZ => rotX e.a1 * rotY e.a2 * rotZ e.a3
  | .XZY => rotX e.a1 * rotZ e.a2 * rotY e.a3
  | .YXZ => rotY e.a1 * rotX e.a2 * rotZ e.a3
  | .YZX => rotY e.a1 * rotZ e.a2 * rotX e.a3
  | .ZXY => rotZ e.a1 * rotX e.a2 * rotY e.a3
  | .ZYX => rotZ e.a1 * rotY e.a2 * rotX e.a3

/-- Convert Euler angles to quaternion. -/
def toQuat (e : Euler) : Quat :=
  Quat.fromMat3 e.toMat3

/-- Extract Euler angles from rotation matrix with given order.
    Note: May have gimbal lock issues when the middle axis is at ±90°. -/
def fromMat3 (m : Mat3) (order : EulerOrder := .XYZ) : Euler :=
  -- Extract based on rotation order
  -- Using standard rotation matrix decomposition formulas
  match order with
  | .XYZ =>
    -- For XYZ: R = Rz * Ry * Rx
    let sy := m.get 0 2
    if Float.abs sy > 0.99999 then
      -- Gimbal lock: y ≈ ±90°
      let a1 := Float.atan2 (-m.get 1 0) (m.get 1 1)
      let a2 := if sy > 0 then Float.halfPi else -Float.halfPi
      let a3 := 0.0
      ⟨a1, a2, a3, .XYZ⟩
    else
      let a1 := Float.atan2 (-m.get 1 2) (m.get 2 2)
      let a2 := Float.asin sy
      let a3 := Float.atan2 (-m.get 0 1) (m.get 0 0)
      ⟨a1, a2, a3, .XYZ⟩

  | .XZY =>
    let sz := -m.get 0 1
    if Float.abs sz > 0.99999 then
      let a1 := Float.atan2 (m.get 2 0) (m.get 2 2)
      let a2 := if sz > 0 then Float.halfPi else -Float.halfPi
      let a3 := 0.0
      ⟨a1, a2, a3, .XZY⟩
    else
      let a1 := Float.atan2 (m.get 2 1) (m.get 1 1)
      let a2 := Float.asin sz
      let a3 := Float.atan2 (m.get 0 2) (m.get 0 0)
      ⟨a1, a2, a3, .XZY⟩

  | .YXZ =>
    let sx := -m.get 1 2
    if Float.abs sx > 0.99999 then
      let a1 := Float.atan2 (-m.get 0 1) (m.get 0 0)
      let a2 := if sx > 0 then Float.halfPi else -Float.halfPi
      let a3 := 0.0
      ⟨a1, a2, a3, .YXZ⟩
    else
      let a1 := Float.atan2 (m.get 0 2) (m.get 2 2)
      let a2 := Float.asin sx
      let a3 := Float.atan2 (m.get 1 0) (m.get 1 1)
      ⟨a1, a2, a3, .YXZ⟩

  | .YZX =>
    let sz := m.get 1 0
    if Float.abs sz > 0.99999 then
      let a1 := Float.atan2 (m.get 0 2) (m.get 0 0)
      let a2 := if sz > 0 then Float.halfPi else -Float.halfPi
      let a3 := 0.0
      ⟨a1, a2, a3, .YZX⟩
    else
      let a1 := Float.atan2 (-m.get 2 0) (m.get 0 0)
      let a2 := Float.asin sz
      let a3 := Float.atan2 (-m.get 1 2) (m.get 1 1)
      ⟨a1, a2, a3, .YZX⟩

  | .ZXY =>
    let sx := m.get 2 1
    if Float.abs sx > 0.99999 then
      let a1 := Float.atan2 (m.get 1 0) (m.get 1 1)
      let a2 := if sx > 0 then Float.halfPi else -Float.halfPi
      let a3 := 0.0
      ⟨a1, a2, a3, .ZXY⟩
    else
      let a1 := Float.atan2 (-m.get 0 1) (m.get 1 1)
      let a2 := Float.asin sx
      let a3 := Float.atan2 (-m.get 2 0) (m.get 2 2)
      ⟨a1, a2, a3, .ZXY⟩

  | .ZYX =>
    let sy := -m.get 2 0
    if Float.abs sy > 0.99999 then
      let a1 := Float.atan2 (-m.get 1 2) (m.get 1 1)
      let a2 := if sy > 0 then Float.halfPi else -Float.halfPi
      let a3 := 0.0
      ⟨a1, a2, a3, .ZYX⟩
    else
      let a1 := Float.atan2 (m.get 1 0) (m.get 0 0)
      let a2 := Float.asin sy
      let a3 := Float.atan2 (m.get 2 1) (m.get 2 2)
      ⟨a1, a2, a3, .ZYX⟩

/-- Extract Euler angles from quaternion with given order. -/
def fromQuat (q : Quat) (order : EulerOrder := .XYZ) : Euler :=
  fromMat3 q.toMat3 order

/-- Check approximate equality. -/
def approxEq (a b : Euler) (eps : Float := 0.0001) : Bool :=
  a.order == b.order &&
  Float.approxEq a.a1 b.a1 eps &&
  Float.approxEq a.a2 b.a2 eps &&
  Float.approxEq a.a3 b.a3 eps

/-- Clamp angles to a range (useful for limiting rotations). -/
def clamp (e : Euler) (minAngles maxAngles : Vec3) : Euler :=
  { e with
    a1 := Float.clamp e.a1 minAngles.x maxAngles.x
    a2 := Float.clamp e.a2 minAngles.y maxAngles.y
    a3 := Float.clamp e.a3 minAngles.z maxAngles.z }

/-- Wrap angles to [-π, π] range. -/
def normalize (e : Euler) : Euler :=
  let wrap (a : Float) : Float :=
    let twoPi := 2.0 * Float.pi
    -- Compute a mod 2π using floor
    let a' := a - twoPi * Float.floor (a / twoPi)
    if a' > Float.pi then a' - twoPi
    else if a' < -Float.pi then a' + twoPi
    else a'
  { e with a1 := wrap e.a1, a2 := wrap e.a2, a3 := wrap e.a3 }

/-- Linear interpolation between Euler angles.
    Note: May not produce shortest path rotation. Use Quat.slerp for that. -/
def lerp (a b : Euler) (t : Float) : Euler :=
  if a.order != b.order then a  -- Can't interpolate different orders
  else
    { a1 := Float.lerp a.a1 b.a1 t
      a2 := Float.lerp a.a2 b.a2 t
      a3 := Float.lerp a.a3 b.a3 t
      order := a.order }

end Euler

-- Extend Quat with Euler utilities
namespace Quat

/-- Create quaternion from Euler angles with explicit order. -/
def fromEulerOrdered (e : Euler) : Quat := e.toQuat

/-- Create quaternion from Euler angles in XYZ order. -/
def fromEulerXYZ (x y z : Float) : Quat := (Euler.xyz x y z).toQuat

/-- Create quaternion from Euler angles in YXZ order (common for games). -/
def fromEulerYXZ (y x z : Float) : Quat := (Euler.yxz y x z).toQuat

/-- Create quaternion from Euler angles in ZYX order (aerospace). -/
def fromEulerZYX (z y x : Float) : Quat := (Euler.zyx z y x).toQuat

/-- Extract Euler angles with given order. -/
def toEulerOrdered (q : Quat) (order : EulerOrder := .XYZ) : Euler :=
  Euler.fromQuat q order

/-- Extract Euler angles in XYZ order. -/
def toEulerXYZ (q : Quat) : Vec3 := (Euler.fromQuat q .XYZ).toRadians

/-- Extract Euler angles in YXZ order. -/
def toEulerYXZ (q : Quat) : Vec3 := (Euler.fromQuat q .YXZ).toRadians

/-- Extract Euler angles in ZYX order. -/
def toEulerZYX (q : Quat) : Vec3 := (Euler.fromQuat q .ZYX).toRadians

end Quat

end Linalg
