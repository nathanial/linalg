/-
  Tests for Vec2 operations.
-/

import Linalg
import Crucible

namespace LinalgTests.Vec2Tests

open Crucible
open Linalg

testSuite "Vec2 Basic"

test "zero vector has all zero components" := do
  Vec2.zero.x ≡ 0.0
  Vec2.zero.y ≡ 0.0

test "one vector has all one components" := do
  Vec2.one.x ≡ 1.0
  Vec2.one.y ≡ 1.0

test "addition works correctly" := do
  let a := Vec2.mk 1.0 2.0
  let b := Vec2.mk 3.0 4.0
  let c := a + b
  c.x ≡ 4.0
  c.y ≡ 6.0

test "subtraction works correctly" := do
  let a := Vec2.mk 5.0 7.0
  let b := Vec2.mk 2.0 3.0
  let c := a - b
  c.x ≡ 3.0
  c.y ≡ 4.0

test "negation works correctly" := do
  let a := Vec2.mk 1.0 (-2.0)
  let b := -a
  b.x ≡ (-1.0)
  b.y ≡ 2.0

test "scalar multiplication works" := do
  let a := Vec2.mk 2.0 3.0
  let b := a * 2.0
  b.x ≡ 4.0
  b.y ≡ 6.0

testSuite "Vec2 Math"

test "dot product works correctly" := do
  let a := Vec2.mk 1.0 2.0
  let b := Vec2.mk 3.0 4.0
  a.dot b ≡ 11.0

test "length calculation is correct" := do
  let v := Vec2.mk 3.0 4.0
  ensure (floatNear v.length 5.0 0.0001) "length should be 5"

test "normalize produces unit length" := do
  let v := Vec2.mk 3.0 4.0
  let n := v.normalize
  ensure (floatNear n.length 1.0 0.0001) "normalized length should be 1"

test "perpendicular is orthogonal" := do
  let v := Vec2.mk 3.0 4.0
  let p := v.perpendicular
  ensure (floatNear (v.dot p) 0.0 0.0001) "perpendicular should have 0 dot product"

test "lerp at 0 returns first vector" := do
  let a := Vec2.mk 0.0 0.0
  let b := Vec2.mk 10.0 10.0
  let c := Vec2.lerp a b 0.0
  c.x ≡ 0.0
  c.y ≡ 0.0

test "lerp at 1 returns second vector" := do
  let a := Vec2.mk 0.0 0.0
  let b := Vec2.mk 10.0 10.0
  let c := Vec2.lerp a b 1.0
  c.x ≡ 10.0
  c.y ≡ 10.0

test "lerp at 0.5 returns midpoint" := do
  let a := Vec2.mk 0.0 0.0
  let b := Vec2.mk 10.0 10.0
  let c := Vec2.lerp a b 0.5
  c.x ≡ 5.0
  c.y ≡ 5.0



end LinalgTests.Vec2Tests
