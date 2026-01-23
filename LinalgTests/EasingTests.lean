/-
  Tests for Easing functions.
-/

import Linalg
import Crucible

namespace LinalgTests.EasingTests

open Crucible
open Linalg

testSuite "Easing Boundaries"

test "linear at 0 returns 0" := do
  ensure (floatNear (Easing.linear 0.0) 0.0 0.0001) "should be 0"

test "linear at 1 returns 1" := do
  ensure (floatNear (Easing.linear 1.0) 1.0 0.0001) "should be 1"

test "quadIn at 0 returns 0" := do
  ensure (floatNear (Easing.quadIn 0.0) 0.0 0.0001) "should be 0"

test "quadIn at 1 returns 1" := do
  ensure (floatNear (Easing.quadIn 1.0) 1.0 0.0001) "should be 1"

test "quadOut at 0 returns 0" := do
  ensure (floatNear (Easing.quadOut 0.0) 0.0 0.0001) "should be 0"

test "quadOut at 1 returns 1" := do
  ensure (floatNear (Easing.quadOut 1.0) 1.0 0.0001) "should be 1"

test "sineIn at 0 returns 0" := do
  ensure (floatNear (Easing.sineIn 0.0) 0.0 0.0001) "should be 0"

test "sineIn at 1 returns 1" := do
  ensure (floatNear (Easing.sineIn 1.0) 1.0 0.0001) "should be 1"

testSuite "Easing Midpoints"

test "linear at 0.5 returns 0.5" := do
  ensure (floatNear (Easing.linear 0.5) 0.5 0.0001) "should be 0.5"

test "quadInOut at 0.5 returns 0.5" := do
  ensure (floatNear (Easing.quadInOut 0.5) 0.5 0.0001) "should be 0.5"

test "cubicInOut at 0.5 returns 0.5" := do
  ensure (floatNear (Easing.cubicInOut 0.5) 0.5 0.0001) "should be 0.5"

test "sineInOut at 0.5 returns 0.5" := do
  ensure (floatNear (Easing.sineInOut 0.5) 0.5 0.0001) "should be 0.5"

testSuite "Easing Properties"

test "quadIn is slower than linear at start" := do
  let t := 0.3
  ensure (Easing.quadIn t < Easing.linear t) "quadIn should be slower at start"

test "quadOut is faster than linear at start" := do
  let t := 0.3
  ensure (Easing.quadOut t > Easing.linear t) "quadOut should be faster at start"

test "smoothStep at 0.5 returns 0.5" := do
  ensure (floatNear (Easing.smoothStep 0.5) 0.5 0.0001) "should be 0.5"

test "smootherStep at 0.5 returns 0.5" := do
  ensure (floatNear (Easing.smootherStep 0.5) 0.5 0.0001) "should be 0.5"

testSuite "Easing Apply"

test "apply interpolates between values" := do
  let result := Easing.apply Easing.linear 10.0 20.0 0.5
  ensure (floatNear result 15.0 0.0001) "should interpolate to 15"

test "apply with quadIn at midpoint" := do
  let result := Easing.apply Easing.quadIn 0.0 100.0 0.5
  ensure (floatNear result 25.0 0.0001) "quadIn at 0.5 should give 25"

testSuite "SmoothDamp"

test "smoothdamp moves toward target" := do
  let state := SmoothDamp.init 0.0
  let (newVal, _) := SmoothDamp.step state 10.0 0.5 0.1
  ensure (newVal > 0.0) "should move toward target"
  ensure (newVal < 10.0) "should not overshoot"

test "smoothdamp reaches target over time" := do
  let mut state := SmoothDamp.init 0.0
  -- Run for 3 seconds (enough to settle with smoothTime=0.3)
  for _ in [:200] do
    let (_, newState) := SmoothDamp.step state 10.0 0.3 0.016
    state := newState
  ensure (floatNear state.current 10.0 0.5) "should be near target after many steps"

test "smoothdamp vec2 moves toward target" := do
  let state := SmoothDampState2.init Vec2.zero
  let (newVal, _) := state.step (Vec2.mk 10.0 10.0) 0.5 0.1
  ensure (newVal.x > 0.0) "x should move toward target"
  ensure (newVal.y > 0.0) "y should move toward target"

test "smoothdamp vec3 moves toward target" := do
  let state := SmoothDampState3.init Vec3.zero
  let (newVal, _) := state.step (Vec3.mk 10.0 10.0 10.0) 0.5 0.1
  ensure (newVal.x > 0.0) "x should move toward target"
  ensure (newVal.z > 0.0) "z should move toward target"

testSuite "Hermite Interpolation"

test "hermite at t=0 returns p0" := do
  let result := Hermite.interpolate 0.0 1.0 10.0 1.0 0.0
  ensure (floatNear result 0.0 0.0001) "should return p0"

test "hermite at t=1 returns p1" := do
  let result := Hermite.interpolate 0.0 1.0 10.0 1.0 1.0
  ensure (floatNear result 10.0 0.0001) "should return p1"

test "hermite at t=0.5 with zero tangents" := do
  -- With zero tangents, hermite should pass through midpoint smoothly
  let result := Hermite.interpolate 0.0 0.0 10.0 0.0 0.5
  ensure (floatNear result 5.0 0.0001) "should be at midpoint"

test "hermite vec2 interpolation" := do
  let p0 := Vec2.zero
  let p1 := Vec2.mk 10.0 10.0
  let m0 := Vec2.mk 5.0 0.0
  let m1 := Vec2.mk 5.0 0.0
  let result := Hermite.interpolateVec2 p0 m0 p1 m1 0.5
  ensure (result.x > 0.0 && result.x < 10.0) "x should be in range"
  ensure (result.y > 0.0 && result.y < 10.0) "y should be in range"

test "hermite vec3 interpolation" := do
  let p0 := Vec3.zero
  let p1 := Vec3.mk 10.0 10.0 10.0
  let result := Hermite.interpolateVec3 p0 Vec3.zero p1 Vec3.zero 0.5
  ensure (floatNear result.x 5.0 0.5) "x should be near midpoint"

test "hermite derivative at endpoints" := do
  -- With tangent m0=2 at t=0, derivative should be 2
  let d := Hermite.derivative 0.0 2.0 10.0 0.0 0.0
  ensure (floatNear d 2.0 0.0001) "derivative at t=0 should equal m0"

test "cardinal spline with tension 0 is catmull-rom" := do
  let result := Hermite.cardinal (-1.0) 0.0 1.0 2.0 0.5 0.0
  -- Should smoothly interpolate through control points
  ensure (result > 0.0 && result < 1.0) "should be between p1 and p2"

testSuite "Animation Curves"

test "empty curve returns 0" := do
  let curve := AnimationCurve.empty
  ensure (floatNear (curve.evaluate 0.5) 0.0 0.0001) "empty curve should return 0"

test "constant curve returns value" := do
  let curve := AnimationCurve.constant 5.0
  ensure (floatNear (curve.evaluate 0.0) 5.0 0.0001) "should return 5"
  ensure (floatNear (curve.evaluate 100.0) 5.0 0.0001) "should return 5 at any time"

test "linear curve interpolates" := do
  let curve := AnimationCurve.linear 0.0 0.0 1.0 10.0
  ensure (floatNear (curve.evaluate 0.0) 0.0 0.0001) "should be 0 at start"
  ensure (floatNear (curve.evaluate 1.0) 10.0 0.0001) "should be 10 at end"
  ensure (floatNear (curve.evaluate 0.5) 5.0 0.0001) "should be 5 at midpoint"

test "curve from points" := do
  let curve := AnimationCurve.fromPoints #[(0.0, 0.0), (1.0, 10.0), (2.0, 5.0)]
  ensure (floatNear (curve.evaluate 0.0) 0.0 0.0001) "should be 0 at t=0"
  ensure (floatNear (curve.evaluate 2.0) 5.0 0.0001) "should be 5 at t=2"

test "curve duration" := do
  let curve := AnimationCurve.linear 1.0 0.0 5.0 10.0
  ensure (floatNear curve.duration 4.0 0.0001) "duration should be 4"
  ensure (floatNear curve.startTime 1.0 0.0001) "start should be 1"
  ensure (floatNear curve.endTime 5.0 0.0001) "end should be 5"

test "curve clamp wrap mode" := do
  let curve := AnimationCurve.linear 0.0 0.0 1.0 10.0
  ensure (floatNear (curve.evaluate (-1.0)) 0.0 0.0001) "before start should clamp"
  ensure (floatNear (curve.evaluate 2.0) 10.0 0.0001) "after end should clamp"

test "curve loop wrap mode" := do
  let curve := (AnimationCurve.linear 0.0 0.0 1.0 10.0).withWrapMode WrapMode.loop WrapMode.loop
  let v1 := curve.evaluate 0.5
  let v2 := curve.evaluate 1.5
  ensure (floatNear v1 v2 0.1) "looped value should repeat"

test "curve add and remove key" := do
  let curve := AnimationCurve.linear 0.0 0.0 2.0 20.0
  let withMid := curve.addKey (Keyframe.smooth 1.0 15.0)
  ensure (withMid.keyframes.size == 3) "should have 3 keys"
  let removed := withMid.removeKey 1
  ensure (removed.keyframes.size == 2) "should have 2 keys after remove"

test "curve sample" := do
  let curve := AnimationCurve.linear 0.0 0.0 1.0 10.0
  let samples := curve.sample 5
  ensure (samples.size == 5) "should have 5 samples"
  ensure (floatNear samples[0]! 0.0 0.0001) "first sample should be 0"
  ensure (floatNear samples[4]! 10.0 0.0001) "last sample should be 10"

test "eased curve" := do
  let curve := AnimationCurve.eased 0.0 0.0 1.0 10.0 Easing.quadIn 5
  ensure (curve.keyframes.size == 6) "should have 6 keyframes"
  let mid := curve.evaluate 0.5
  -- quadIn at 0.5 gives 0.25, so value should be 2.5
  ensure (mid < 5.0) "eased curve should be slower at start"



end LinalgTests.EasingTests
