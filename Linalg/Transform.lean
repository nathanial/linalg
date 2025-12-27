/-
  Transform type combining position, rotation, and scale.
  Commonly used in game engines for scene graph hierarchies.
-/

import Linalg.Vec3
import Linalg.Quat
import Linalg.Mat4

namespace Linalg

/-- A transform representing position, rotation, and scale.
    The order of operations is: Scale → Rotate → Translate (SRT). -/
structure Transform where
  position : Vec3
  rotation : Quat
  scale : Vec3
  deriving Repr

namespace Transform

/-- Identity transform: no translation, no rotation, uniform scale of 1. -/
def identity : Transform :=
  { position := Vec3.zero
    rotation := Quat.identity
    scale := Vec3.one }

/-- Create a transform with just a position. -/
def fromPosition (pos : Vec3) : Transform :=
  { identity with position := pos }

/-- Create a transform with just a rotation. -/
def fromRotation (rot : Quat) : Transform :=
  { identity with rotation := rot }

/-- Create a transform with just a scale. -/
def fromScale (s : Vec3) : Transform :=
  { identity with scale := s }

/-- Create a transform with uniform scale. -/
def fromUniformScale (s : Float) : Transform :=
  { identity with scale := Vec3.mk s s s }

/-- Create a transform from position and rotation. -/
def fromPositionRotation (pos : Vec3) (rot : Quat) : Transform :=
  { position := pos, rotation := rot, scale := Vec3.one }

/-- Create a transform from all components. -/
def mk' (pos : Vec3) (rot : Quat) (scale : Vec3) : Transform :=
  { position := pos, rotation := rot, scale := scale }

/-- Get the forward direction (negative Z in local space). -/
def forward (t : Transform) : Vec3 :=
  t.rotation * Vec3.forward

/-- Get the right direction (positive X in local space). -/
def right (t : Transform) : Vec3 :=
  t.rotation * Vec3.unitX

/-- Get the up direction (positive Y in local space). -/
def up (t : Transform) : Vec3 :=
  t.rotation * Vec3.unitY

/-- Get the backward direction (positive Z in local space). -/
def backward (t : Transform) : Vec3 :=
  t.rotation * Vec3.unitZ

/-- Get the left direction (negative X in local space). -/
def left (t : Transform) : Vec3 :=
  t.rotation * Vec3.left

/-- Get the down direction (negative Y in local space). -/
def down (t : Transform) : Vec3 :=
  t.rotation * Vec3.down

/-- Convert transform to a 4x4 matrix (TRS order). -/
def toMat4 (t : Transform) : Mat4 :=
  let translationMat := Mat4.translation t.position.x t.position.y t.position.z
  let rotationMat := t.rotation.toMat4
  let scaleMat := Mat4.scaling t.scale.x t.scale.y t.scale.z
  translationMat * rotationMat * scaleMat

/-- Convert transform to a 4x4 matrix without scale. -/
def toMat4NoScale (t : Transform) : Mat4 :=
  let translationMat := Mat4.translation t.position.x t.position.y t.position.z
  let rotationMat := t.rotation.toMat4
  translationMat * rotationMat

/-- Extract transform from a 4x4 matrix.
    Note: This assumes the matrix represents a valid TRS transform with no shear. -/
def fromMat4 (m : Mat4) : Transform :=
  -- Extract translation from last column
  let position := m.getTranslation

  -- Extract scale from column magnitudes
  let col0 := Vec3.mk (m.get 0 0) (m.get 1 0) (m.get 2 0)
  let col1 := Vec3.mk (m.get 0 1) (m.get 1 1) (m.get 2 1)
  let col2 := Vec3.mk (m.get 0 2) (m.get 1 2) (m.get 2 2)
  let scaleX := col0.length
  let scaleY := col1.length
  let scaleZ := col2.length
  let scale := Vec3.mk scaleX scaleY scaleZ

  -- Normalize columns to get rotation matrix
  let invScaleX := if scaleX > Float.epsilon then 1.0 / scaleX else 0.0
  let invScaleY := if scaleY > Float.epsilon then 1.0 / scaleY else 0.0
  let invScaleZ := if scaleZ > Float.epsilon then 1.0 / scaleZ else 0.0

  let rotMat := Mat3.fromColumns
    (col0.scale invScaleX)
    (col1.scale invScaleY)
    (col2.scale invScaleZ)

  let rotation := Quat.fromMat3 rotMat

  { position := position, rotation := rotation, scale := scale }

/-- Transform a point from local to world space. -/
def transformPoint (t : Transform) (p : Vec3) : Vec3 :=
  -- Scale, then rotate, then translate
  let scaled := Vec3.mk (p.x * t.scale.x) (p.y * t.scale.y) (p.z * t.scale.z)
  let rotated := t.rotation * scaled
  rotated.add t.position

/-- Transform a direction from local to world space (ignores translation). -/
def transformDirection (t : Transform) (d : Vec3) : Vec3 :=
  -- Scale, then rotate
  let scaled := Vec3.mk (d.x * t.scale.x) (d.y * t.scale.y) (d.z * t.scale.z)
  t.rotation * scaled

/-- Transform a direction without scale (useful for normals). -/
def transformDirectionNoScale (t : Transform) (d : Vec3) : Vec3 :=
  t.rotation * d

/-- Inverse transform a point from world to local space. -/
def inverseTransformPoint (t : Transform) (p : Vec3) : Vec3 :=
  -- Inverse: translate, then inverse rotate, then inverse scale
  let translated := p.sub t.position
  let rotated := t.rotation.conjugate * translated
  let invScaleX := if t.scale.x.abs > Float.epsilon then 1.0 / t.scale.x else 0.0
  let invScaleY := if t.scale.y.abs > Float.epsilon then 1.0 / t.scale.y else 0.0
  let invScaleZ := if t.scale.z.abs > Float.epsilon then 1.0 / t.scale.z else 0.0
  Vec3.mk (rotated.x * invScaleX) (rotated.y * invScaleY) (rotated.z * invScaleZ)

/-- Inverse transform a direction from world to local space. -/
def inverseTransformDirection (t : Transform) (d : Vec3) : Vec3 :=
  let rotated := t.rotation.conjugate * d
  let invScaleX := if t.scale.x.abs > Float.epsilon then 1.0 / t.scale.x else 0.0
  let invScaleY := if t.scale.y.abs > Float.epsilon then 1.0 / t.scale.y else 0.0
  let invScaleZ := if t.scale.z.abs > Float.epsilon then 1.0 / t.scale.z else 0.0
  Vec3.mk (rotated.x * invScaleX) (rotated.y * invScaleY) (rotated.z * invScaleZ)

/-- Compute the inverse transform. -/
def inverse (t : Transform) : Transform :=
  let invRot := t.rotation.conjugate
  let invScaleX := if t.scale.x.abs > Float.epsilon then 1.0 / t.scale.x else 0.0
  let invScaleY := if t.scale.y.abs > Float.epsilon then 1.0 / t.scale.y else 0.0
  let invScaleZ := if t.scale.z.abs > Float.epsilon then 1.0 / t.scale.z else 0.0
  let invScale := Vec3.mk invScaleX invScaleY invScaleZ
  -- Inverse position: -inverse(R*S) * position = -S^-1 * R^-1 * position
  let rotatedPos := invRot * t.position
  let invPos := Vec3.mk
    (-rotatedPos.x * invScaleX)
    (-rotatedPos.y * invScaleY)
    (-rotatedPos.z * invScaleZ)
  { position := invPos, rotation := invRot, scale := invScale }

/-- Combine two transforms: child transform in parent's space.
    Result represents first applying child, then parent. -/
def mul (parent child : Transform) : Transform :=
  -- Combined position: parent.pos + parent.rot * (parent.scale * child.pos)
  let scaledChildPos := Vec3.mk
    (child.position.x * parent.scale.x)
    (child.position.y * parent.scale.y)
    (child.position.z * parent.scale.z)
  let position := parent.position.add (parent.rotation * scaledChildPos)

  -- Combined rotation: parent.rot * child.rot
  let rotation := parent.rotation * child.rotation

  -- Combined scale: parent.scale * child.scale
  let scale := Vec3.mk
    (parent.scale.x * child.scale.x)
    (parent.scale.y * child.scale.y)
    (parent.scale.z * child.scale.z)

  { position := position, rotation := rotation, scale := scale }

instance : HMul Transform Transform Transform where
  hMul := mul

/-- Interpolate between two transforms. -/
def lerp (a b : Transform) (t : Float) : Transform :=
  { position := a.position.lerp b.position t
    rotation := Quat.slerp a.rotation b.rotation t
    scale := a.scale.lerp b.scale t }

/-- Look at a target point, setting rotation to face it. -/
def lookAt (t : Transform) (target : Vec3) (worldUp : Vec3 := Vec3.unitY) : Transform :=
  let direction := (target.sub t.position).normalize
  let rotation := Quat.lookRotation direction worldUp
  { t with rotation := rotation }

/-- Translate by a world-space offset. -/
def translate (t : Transform) (offset : Vec3) : Transform :=
  { t with position := t.position.add offset }

/-- Translate by a local-space offset. -/
def translateLocal (t : Transform) (offset : Vec3) : Transform :=
  let worldOffset := t.rotation * offset
  { t with position := t.position.add worldOffset }

/-- Rotate by a quaternion (in world space). -/
def rotate (t : Transform) (rot : Quat) : Transform :=
  { t with rotation := rot * t.rotation }

/-- Rotate by a quaternion (in local space). -/
def rotateLocal (t : Transform) (rot : Quat) : Transform :=
  { t with rotation := t.rotation * rot }

/-- Rotate around an axis by an angle (in radians, world space). -/
def rotateAround (t : Transform) (axis : Vec3) (angle : Float) : Transform :=
  t.rotate (Quat.fromAxisAngle axis angle)

/-- Scale uniformly. -/
def scaleUniform (t : Transform) (s : Float) : Transform :=
  { t with scale := t.scale.scale s }

/-- Scale non-uniformly. -/
def scaleBy (t : Transform) (s : Vec3) : Transform :=
  { t with scale := Vec3.mk (t.scale.x * s.x) (t.scale.y * s.y) (t.scale.z * s.z) }

/-- Check if two transforms are approximately equal. -/
def approxEq (a b : Transform) (epsilon : Float := 0.0001) : Bool :=
  a.position.distance b.position < epsilon &&
  (a.rotation.sameRotation b.rotation) &&
  a.scale.distance b.scale < epsilon

end Transform

/-! ## Transform Hierarchy

    Utilities for managing parent-child transform relationships.
    Supports computing world transforms from local transforms and vice versa.
-/

/-- A node in a transform hierarchy with optional parent reference.
    The transform is in local space relative to the parent. -/
structure TransformNode where
  /-- Local transform relative to parent (or world if no parent). -/
  local_ : Transform
  /-- Optional parent node index (for external hierarchy storage). -/
  parentIdx : Option Nat
  deriving Repr

instance : Inhabited TransformNode where
  default := { local_ := Transform.identity, parentIdx := none }

namespace TransformNode

/-- Create a root node (no parent). -/
def root (t : Transform) : TransformNode :=
  { local_ := t, parentIdx := none }

/-- Create a child node with parent index. -/
def child (t : Transform) (parent : Nat) : TransformNode :=
  { local_ := t, parentIdx := some parent }

/-- Check if this is a root node. -/
def isRoot (n : TransformNode) : Bool := n.parentIdx.isNone

end TransformNode

/-- A simple transform hierarchy storing nodes with parent references.
    Nodes are stored in an array, with parent indices referring to earlier entries.
    Index 0 is typically the root. -/
structure TransformHierarchy where
  nodes : Array TransformNode
  deriving Repr, Inhabited

namespace TransformHierarchy

/-- Create an empty hierarchy. -/
def empty : TransformHierarchy := { nodes := #[] }

/-- Create a hierarchy with a single root node. -/
def singleton (t : Transform) : TransformHierarchy :=
  { nodes := #[TransformNode.root t] }

/-- Add a root node to the hierarchy. Returns (hierarchy, node index). -/
def addRoot (h : TransformHierarchy) (t : Transform) : TransformHierarchy × Nat :=
  let idx := h.nodes.size
  ({ nodes := h.nodes.push (TransformNode.root t) }, idx)

/-- Add a child node to the hierarchy. Returns (hierarchy, node index).
    Fails if parent index is invalid. -/
def addChild (h : TransformHierarchy) (t : Transform) (parentIdx : Nat) : Option (TransformHierarchy × Nat) :=
  if parentIdx < h.nodes.size then
    let idx := h.nodes.size
    some ({ nodes := h.nodes.push (TransformNode.child t parentIdx) }, idx)
  else
    none

/-- Get a node by index. -/
def getNode (h : TransformHierarchy) (idx : Nat) : Option TransformNode :=
  h.nodes[idx]?

/-- Get the local transform of a node. -/
def getLocal (h : TransformHierarchy) (idx : Nat) : Option Transform :=
  h.nodes[idx]?.map (·.local_)

/-- Set the local transform of a node. -/
def setLocal (h : TransformHierarchy) (idx : Nat) (t : Transform) : TransformHierarchy :=
  match h.nodes[idx]? with
  | some node => { nodes := h.nodes.set! idx { node with local_ := t } }
  | none => h

/-- Compute the world transform for a node by walking up the parent chain. -/
def getWorld (h : TransformHierarchy) (idx : Nat) : Option Transform :=
  if idx >= h.nodes.size then none
  else
    -- Walk up parent chain, accumulating transforms
    let rec go (i : Nat) (acc : Transform) (fuel : Nat) : Transform :=
      if fuel == 0 then acc
      else
        match h.nodes[i]? with
        | none => acc
        | some node =>
          let worldAcc := node.local_ * acc
          match node.parentIdx with
          | none => worldAcc
          | some pi =>
            if pi < i then go pi worldAcc (fuel - 1)
            else worldAcc  -- Invalid parent (must be earlier in array)
    some (go idx Transform.identity h.nodes.size)

/-- Convert a world-space transform to local-space given parent's world transform. -/
def worldToLocal (parentWorld : Transform) (worldT : Transform) : Transform :=
  parentWorld.inverse * worldT

/-- Convert a local-space transform to world-space given parent's world transform. -/
def localToWorld (parentWorld : Transform) (localT : Transform) : Transform :=
  parentWorld * localT

/-- Get all children of a node. -/
def getChildren (h : TransformHierarchy) (parentIdx : Nat) : Array Nat :=
  (List.range h.nodes.size).foldl (init := #[]) fun acc i =>
    match h.nodes[i]? with
    | some node =>
      if node.parentIdx == some parentIdx then acc.push i else acc
    | none => acc

/-- Get the depth of a node (0 for root). -/
def getDepth (h : TransformHierarchy) (idx : Nat) : Nat :=
  let rec go (i : Nat) (depth : Nat) (fuel : Nat) : Nat :=
    if fuel == 0 then depth
    else
      match h.nodes[i]? with
      | none => depth
      | some node =>
        match node.parentIdx with
        | none => depth
        | some pi =>
          if pi < i then go pi (depth + 1) (fuel - 1)
          else depth
  go idx 0 h.nodes.size

/-- Number of nodes in the hierarchy. -/
def size (h : TransformHierarchy) : Nat := h.nodes.size

/-- Check if hierarchy is empty. -/
def isEmpty (h : TransformHierarchy) : Bool := h.nodes.isEmpty

end TransformHierarchy

/-! ## Transform Chain Utilities

    Functions for working with chains of transforms without a full hierarchy structure.
-/

namespace Transform

/-- Compute the world transform from a chain of local transforms (root first). -/
def chainToWorld (chain : Array Transform) : Transform :=
  chain.foldl (· * ·) Transform.identity

/-- Given a parent's world transform and a desired world transform,
    compute the required local transform for a child. -/
def computeLocalTransform (parentWorld desiredWorld : Transform) : Transform :=
  parentWorld.inverse * desiredWorld

/-- Given the local transforms of ancestors (root first) and a local transform,
    compute the world transform. -/
def localToWorldChain (ancestors : Array Transform) (localT : Transform) : Transform :=
  let parentWorld := ancestors.foldl (· * ·) Transform.identity
  parentWorld * localT

/-- Reparent a transform: given its current world transform and new parent's world transform,
    compute the new local transform to maintain the same world position. -/
def reparent (currentWorld newParentWorld : Transform) : Transform :=
  newParentWorld.inverse * currentWorld

end Transform

end Linalg
