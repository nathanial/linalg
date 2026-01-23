/-
  Tests for Transform type.
-/

import Linalg
import Crucible

namespace LinalgTests.TransformTests

open Crucible
open Linalg

testSuite "Transform Basic"

test "identity transform" := do
  let t := Transform.identity
  ensure (floatNear t.position.x 0.0 0.0001) "position x should be 0"
  ensure (floatNear t.position.y 0.0 0.0001) "position y should be 0"
  ensure (floatNear t.position.z 0.0 0.0001) "position z should be 0"
  ensure (floatNear t.scale.x 1.0 0.0001) "scale x should be 1"
  ensure (t.rotation.sameRotation Quat.identity) "rotation should be identity"

test "fromPosition creates correct transform" := do
  let t := Transform.fromPosition (Vec3.mk 5.0 10.0 15.0)
  ensure (floatNear t.position.x 5.0 0.0001) "position x"
  ensure (floatNear t.position.y 10.0 0.0001) "position y"
  ensure (floatNear t.position.z 15.0 0.0001) "position z"

test "fromUniformScale creates correct transform" := do
  let t := Transform.fromUniformScale 2.0
  ensure (floatNear t.scale.x 2.0 0.0001) "scale x"
  ensure (floatNear t.scale.y 2.0 0.0001) "scale y"
  ensure (floatNear t.scale.z 2.0 0.0001) "scale z"

testSuite "Transform Point"

test "identity transform leaves point unchanged" := do
  let t := Transform.identity
  let p := Vec3.mk 1.0 2.0 3.0
  let result := t.transformPoint p
  ensure (floatNear result.x 1.0 0.0001) "x"
  ensure (floatNear result.y 2.0 0.0001) "y"
  ensure (floatNear result.z 3.0 0.0001) "z"

test "translation moves point" := do
  let t := Transform.fromPosition (Vec3.mk 10.0 0.0 0.0)
  let p := Vec3.zero
  let result := t.transformPoint p
  ensure (floatNear result.x 10.0 0.0001) "x should be 10"

test "scale scales point" := do
  let t := Transform.fromUniformScale 2.0
  let p := Vec3.mk 1.0 1.0 1.0
  let result := t.transformPoint p
  ensure (floatNear result.x 2.0 0.0001) "x should be 2"
  ensure (floatNear result.y 2.0 0.0001) "y should be 2"
  ensure (floatNear result.z 2.0 0.0001) "z should be 2"

test "rotation rotates point" := do
  let t := Transform.fromRotation (Quat.fromAxisAngle Vec3.unitY Float.halfPi)
  let p := Vec3.unitX
  let result := t.transformPoint p
  ensure (floatNear result.x 0.0 0.0001) "x should be 0"
  ensure (floatNear result.z (-1.0) 0.0001) "z should be -1"

testSuite "Transform Inverse"

test "inverse of identity is identity" := do
  let t := Transform.identity
  let inv := t.inverse
  ensure (inv.approxEq Transform.identity) "inverse of identity should be identity"

test "inverse undoes transform" := do
  let t := Transform.mk' (Vec3.mk 5.0 10.0 15.0) Quat.identity (Vec3.mk 2.0 2.0 2.0)
  let inv := t.inverse
  let p := Vec3.mk 1.0 2.0 3.0
  let transformed := t.transformPoint p
  let result := inv.transformPoint transformed
  ensure (floatNear result.x p.x 0.001) "x should return to original"
  ensure (floatNear result.y p.y 0.001) "y should return to original"
  ensure (floatNear result.z p.z 0.001) "z should return to original"

testSuite "Transform Composition"

test "composing with identity gives same transform" := do
  let t := Transform.fromPosition (Vec3.mk 5.0 0.0 0.0)
  let result := t * Transform.identity
  ensure (result.approxEq t) "composing with identity should give same"

test "translation composition" := do
  let t1 := Transform.fromPosition (Vec3.mk 5.0 0.0 0.0)
  let t2 := Transform.fromPosition (Vec3.mk 0.0 10.0 0.0)
  let result := t1 * t2
  ensure (floatNear result.position.x 5.0 0.0001) "x should be 5"
  ensure (floatNear result.position.y 10.0 0.0001) "y should be 10"

testSuite "Transform Lerp"

test "lerp at 0 returns first transform" := do
  let a := Transform.fromPosition Vec3.zero
  let b := Transform.fromPosition (Vec3.mk 10.0 0.0 0.0)
  let result := Transform.lerp a b 0.0
  ensure (floatNear result.position.x 0.0 0.0001) "should be at first position"

test "lerp at 1 returns second transform" := do
  let a := Transform.fromPosition Vec3.zero
  let b := Transform.fromPosition (Vec3.mk 10.0 0.0 0.0)
  let result := Transform.lerp a b 1.0
  ensure (floatNear result.position.x 10.0 0.0001) "should be at second position"

test "lerp at 0.5 returns midpoint" := do
  let a := Transform.fromPosition Vec3.zero
  let b := Transform.fromPosition (Vec3.mk 10.0 0.0 0.0)
  let result := Transform.lerp a b 0.5
  ensure (floatNear result.position.x 5.0 0.0001) "should be at midpoint"

testSuite "Transform Hierarchy"

test "empty hierarchy" := do
  let h := TransformHierarchy.empty
  ensure (h.size == 0) "should be empty"
  ensure h.isEmpty "isEmpty should be true"

test "singleton hierarchy" := do
  let h := TransformHierarchy.singleton Transform.identity
  ensure (h.size == 1) "should have 1 node"
  ensure (!h.isEmpty) "isEmpty should be false"

test "add root node" := do
  let h := TransformHierarchy.empty
  let (h2, idx) := h.addRoot (Transform.fromPosition (Vec3.mk 5.0 0.0 0.0))
  ensure (idx == 0) "first node should be index 0"
  ensure (h2.size == 1) "should have 1 node"

test "add child node" := do
  let h := TransformHierarchy.singleton Transform.identity
  match h.addChild (Transform.fromPosition (Vec3.mk 1.0 0.0 0.0)) 0 with
  | some (h2, idx) =>
    ensure (idx == 1) "child should be index 1"
    ensure (h2.size == 2) "should have 2 nodes"
  | none => ensure false "addChild should succeed"

test "add child with invalid parent fails" := do
  let h := TransformHierarchy.singleton Transform.identity
  match h.addChild Transform.identity 99 with
  | some _ => ensure false "should fail with invalid parent"
  | none => pure ()

test "getLocal retrieves transform" := do
  let t := Transform.fromPosition (Vec3.mk 3.0 4.0 5.0)
  let h := TransformHierarchy.singleton t
  match h.getLocal 0 with
  | some t2 => ensure (t.approxEq t2) "should match original"
  | none => ensure false "getLocal should succeed"

test "setLocal updates transform" := do
  let h := TransformHierarchy.singleton Transform.identity
  let newT := Transform.fromPosition (Vec3.mk 10.0 0.0 0.0)
  let h2 := h.setLocal 0 newT
  match h2.getLocal 0 with
  | some t => ensure (t.approxEq newT) "should be updated"
  | none => ensure false "getLocal should succeed"

test "getWorld for root returns local" := do
  let t := Transform.fromPosition (Vec3.mk 5.0 0.0 0.0)
  let h := TransformHierarchy.singleton t
  match h.getWorld 0 with
  | some world => ensure (world.approxEq t) "world should equal local for root"
  | none => ensure false "getWorld should succeed"

test "getWorld combines parent and child" := do
  let parent := Transform.fromPosition (Vec3.mk 10.0 0.0 0.0)
  let child := Transform.fromPosition (Vec3.mk 5.0 0.0 0.0)
  let (h, _) := TransformHierarchy.empty.addRoot parent
  match h.addChild child 0 with
  | some (h2, childIdx) =>
    match h2.getWorld childIdx with
    | some world =>
      -- World position should be 10 + 5 = 15
      ensure (floatNear world.position.x 15.0 0.001) "x should be 15"
    | none => ensure false "getWorld should succeed"
  | none => ensure false "addChild should succeed"

test "getChildren finds children" := do
  let (h, rootIdx) := TransformHierarchy.empty.addRoot Transform.identity
  match h.addChild Transform.identity rootIdx with
  | some (h2, _) =>
    match h2.addChild Transform.identity rootIdx with
    | some (h3, _) =>
      let children := h3.getChildren rootIdx
      ensure (children.size == 2) "should have 2 children"
    | none => ensure false "addChild should succeed"
  | none => ensure false "addChild should succeed"

test "getDepth returns correct depth" := do
  let (h, _) := TransformHierarchy.empty.addRoot Transform.identity
  match h.addChild Transform.identity 0 with
  | some (h2, idx1) =>
    match h2.addChild Transform.identity idx1 with
    | some (h3, idx2) =>
      ensure (h3.getDepth 0 == 0) "root depth should be 0"
      ensure (h3.getDepth idx1 == 1) "child depth should be 1"
      ensure (h3.getDepth idx2 == 2) "grandchild depth should be 2"
    | none => ensure false "addChild should succeed"
  | none => ensure false "addChild should succeed"

testSuite "Transform Chain Utilities"

test "chainToWorld combines transforms" := do
  let t1 := Transform.fromPosition (Vec3.mk 5.0 0.0 0.0)
  let t2 := Transform.fromPosition (Vec3.mk 3.0 0.0 0.0)
  let result := Transform.chainToWorld #[t1, t2]
  ensure (floatNear result.position.x 8.0 0.001) "x should be 5+3=8"

test "computeLocalTransform gives correct local" := do
  let parentWorld := Transform.fromPosition (Vec3.mk 10.0 0.0 0.0)
  let desiredWorld := Transform.fromPosition (Vec3.mk 15.0 0.0 0.0)
  let local_ := Transform.computeLocalTransform parentWorld desiredWorld
  ensure (floatNear local_.position.x 5.0 0.001) "local x should be 5"

test "localToWorldChain computes world" := do
  let ancestors := #[
    Transform.fromPosition (Vec3.mk 10.0 0.0 0.0),
    Transform.fromPosition (Vec3.mk 5.0 0.0 0.0)
  ]
  let local_ := Transform.fromPosition (Vec3.mk 2.0 0.0 0.0)
  let world := Transform.localToWorldChain ancestors local_
  ensure (floatNear world.position.x 17.0 0.001) "x should be 10+5+2=17"

test "reparent maintains world position" := do
  let currentWorld := Transform.fromPosition (Vec3.mk 20.0 0.0 0.0)
  let newParentWorld := Transform.fromPosition (Vec3.mk 5.0 0.0 0.0)
  let newLocal := Transform.reparent currentWorld newParentWorld
  -- The new local should be such that newParent * newLocal = currentWorld
  let recomputed := newParentWorld * newLocal
  ensure (floatNear recomputed.position.x 20.0 0.001) "world should be preserved"



end LinalgTests.TransformTests
