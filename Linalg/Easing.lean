/-
  Easing functions for animation and interpolation.

  All easing functions take a parameter t in [0, 1] and return a value
  that represents the eased progression. The functions follow the naming
  convention: easeIn (slow start), easeOut (slow end), easeInOut (slow both).
-/

import Linalg.Core
import Linalg.Vec2
import Linalg.Vec3

namespace Linalg

namespace Easing

/-- Linear interpolation (no easing). -/
@[inline]
def linear (t : Float) : Float := t

-- ============================================================================
-- Quadratic (t²)
-- ============================================================================

/-- Quadratic ease in: starts slow, accelerates. -/
@[inline]
def quadIn (t : Float) : Float := t * t

/-- Quadratic ease out: starts fast, decelerates. -/
@[inline]
def quadOut (t : Float) : Float :=
  let inv := 1.0 - t
  1.0 - inv * inv

/-- Quadratic ease in-out: slow start and end. -/
@[inline]
def quadInOut (t : Float) : Float :=
  if t < 0.5 then
    2.0 * t * t
  else
    let inv := -2.0 * t + 2.0
    1.0 - inv * inv / 2.0

-- ============================================================================
-- Cubic (t³)
-- ============================================================================

/-- Cubic ease in. -/
@[inline]
def cubicIn (t : Float) : Float := t * t * t

/-- Cubic ease out. -/
@[inline]
def cubicOut (t : Float) : Float :=
  let inv := 1.0 - t
  1.0 - inv * inv * inv

/-- Cubic ease in-out. -/
@[inline]
def cubicInOut (t : Float) : Float :=
  if t < 0.5 then
    4.0 * t * t * t
  else
    let inv := -2.0 * t + 2.0
    1.0 - inv * inv * inv / 2.0

-- ============================================================================
-- Quartic (t⁴)
-- ============================================================================

/-- Quartic ease in. -/
@[inline]
def quartIn (t : Float) : Float := t * t * t * t

/-- Quartic ease out. -/
@[inline]
def quartOut (t : Float) : Float :=
  let inv := 1.0 - t
  1.0 - inv * inv * inv * inv

/-- Quartic ease in-out. -/
@[inline]
def quartInOut (t : Float) : Float :=
  if t < 0.5 then
    8.0 * t * t * t * t
  else
    let inv := -2.0 * t + 2.0
    1.0 - inv * inv * inv * inv / 2.0

-- ============================================================================
-- Quintic (t⁵)
-- ============================================================================

/-- Quintic ease in. -/
@[inline]
def quintIn (t : Float) : Float := t * t * t * t * t

/-- Quintic ease out. -/
@[inline]
def quintOut (t : Float) : Float :=
  let inv := 1.0 - t
  1.0 - inv * inv * inv * inv * inv

/-- Quintic ease in-out. -/
@[inline]
def quintInOut (t : Float) : Float :=
  if t < 0.5 then
    16.0 * t * t * t * t * t
  else
    let inv := -2.0 * t + 2.0
    1.0 - inv * inv * inv * inv * inv / 2.0

-- ============================================================================
-- Sinusoidal
-- ============================================================================

/-- Sine ease in. -/
@[inline]
def sineIn (t : Float) : Float :=
  1.0 - Float.cos (t * Float.halfPi)

/-- Sine ease out. -/
@[inline]
def sineOut (t : Float) : Float :=
  Float.sin (t * Float.halfPi)

/-- Sine ease in-out. -/
@[inline]
def sineInOut (t : Float) : Float :=
  -(Float.cos (Float.pi * t) - 1.0) / 2.0

-- ============================================================================
-- Exponential
-- ============================================================================

/-- Exponential ease in. -/
@[inline]
def expoIn (t : Float) : Float :=
  if t == 0.0 then 0.0 else Float.pow 2.0 (10.0 * t - 10.0)

/-- Exponential ease out. -/
@[inline]
def expoOut (t : Float) : Float :=
  if t == 1.0 then 1.0 else 1.0 - Float.pow 2.0 (-10.0 * t)

/-- Exponential ease in-out. -/
@[inline]
def expoInOut (t : Float) : Float :=
  if t == 0.0 then 0.0
  else if t == 1.0 then 1.0
  else if t < 0.5 then
    Float.pow 2.0 (20.0 * t - 10.0) / 2.0
  else
    (2.0 - Float.pow 2.0 (-20.0 * t + 10.0)) / 2.0

-- ============================================================================
-- Circular
-- ============================================================================

/-- Circular ease in. -/
@[inline]
def circIn (t : Float) : Float :=
  1.0 - Float.sqrt (1.0 - t * t)

/-- Circular ease out. -/
@[inline]
def circOut (t : Float) : Float :=
  let tm1 := t - 1.0
  Float.sqrt (1.0 - tm1 * tm1)

/-- Circular ease in-out. -/
@[inline]
def circInOut (t : Float) : Float :=
  if t < 0.5 then
    (1.0 - Float.sqrt (1.0 - 4.0 * t * t)) / 2.0
  else
    let v := -2.0 * t + 2.0
    (Float.sqrt (1.0 - v * v) + 1.0) / 2.0

-- ============================================================================
-- Back (overshoots)
-- ============================================================================

private def backC1 : Float := 1.70158
private def backC2 : Float := backC1 * 1.525
private def backC3 : Float := backC1 + 1.0

/-- Back ease in: overshoots at the start. -/
@[inline]
def backIn (t : Float) : Float :=
  backC3 * t * t * t - backC1 * t * t

/-- Back ease out: overshoots at the end. -/
@[inline]
def backOut (t : Float) : Float :=
  let tm1 := t - 1.0
  1.0 + backC3 * tm1 * tm1 * tm1 + backC1 * tm1 * tm1

/-- Back ease in-out: overshoots at both ends. -/
@[inline]
def backInOut (t : Float) : Float :=
  if t < 0.5 then
    let v := 2.0 * t
    (v * v * ((backC2 + 1.0) * v - backC2)) / 2.0
  else
    let v := 2.0 * t - 2.0
    (v * v * ((backC2 + 1.0) * v + backC2) + 2.0) / 2.0

-- ============================================================================
-- Elastic (springy)
-- ============================================================================

private def elasticC4 : Float := Float.twoPi / 3.0
private def elasticC5 : Float := Float.twoPi / 4.5

/-- Elastic ease in: springy start. -/
@[inline]
def elasticIn (t : Float) : Float :=
  if t == 0.0 then 0.0
  else if t == 1.0 then 1.0
  else
    -Float.pow 2.0 (10.0 * t - 10.0) * Float.sin ((t * 10.0 - 10.75) * elasticC4)

/-- Elastic ease out: springy end. -/
@[inline]
def elasticOut (t : Float) : Float :=
  if t == 0.0 then 0.0
  else if t == 1.0 then 1.0
  else
    Float.pow 2.0 (-10.0 * t) * Float.sin ((t * 10.0 - 0.75) * elasticC4) + 1.0

/-- Elastic ease in-out: springy both ends. -/
@[inline]
def elasticInOut (t : Float) : Float :=
  if t == 0.0 then 0.0
  else if t == 1.0 then 1.0
  else if t < 0.5 then
    -(Float.pow 2.0 (20.0 * t - 10.0) * Float.sin ((20.0 * t - 11.125) * elasticC5)) / 2.0
  else
    Float.pow 2.0 (-20.0 * t + 10.0) * Float.sin ((20.0 * t - 11.125) * elasticC5) / 2.0 + 1.0

-- ============================================================================
-- Bounce
-- ============================================================================

private def bounceN1 : Float := 7.5625
private def bounceD1 : Float := 2.75

/-- Bounce ease out: bouncing at the end. -/
def bounceOut (t : Float) : Float :=
  if t < 1.0 / bounceD1 then
    bounceN1 * t * t
  else if t < 2.0 / bounceD1 then
    let t' := t - 1.5 / bounceD1
    bounceN1 * t' * t' + 0.75
  else if t < 2.5 / bounceD1 then
    let t' := t - 2.25 / bounceD1
    bounceN1 * t' * t' + 0.9375
  else
    let t' := t - 2.625 / bounceD1
    bounceN1 * t' * t' + 0.984375

/-- Bounce ease in: bouncing at the start. -/
@[inline]
def bounceIn (t : Float) : Float :=
  1.0 - bounceOut (1.0 - t)

/-- Bounce ease in-out: bouncing at both ends. -/
@[inline]
def bounceInOut (t : Float) : Float :=
  if t < 0.5 then
    (1.0 - bounceOut (1.0 - 2.0 * t)) / 2.0
  else
    (1.0 + bounceOut (2.0 * t - 1.0)) / 2.0

-- ============================================================================
-- Smooth step variants
-- ============================================================================

/-- Hermite smooth step (3t² - 2t³). -/
@[inline]
def smoothStep (t : Float) : Float :=
  t * t * (3.0 - 2.0 * t)

/-- Perlin's improved smooth step (6t⁵ - 15t⁴ + 10t³). -/
@[inline]
def smootherStep (t : Float) : Float :=
  t * t * t * (t * (t * 6.0 - 15.0) + 10.0)

-- ============================================================================
-- Spring interpolation
-- ============================================================================

/-- Spring interpolation with configurable damping and frequency.
    Returns a value that oscillates around the target before settling. -/
def spring (t : Float) (damping : Float := 0.5) (frequency : Float := 10.0) : Float :=
  if t <= 0.0 then 0.0
  else if t >= 1.0 then 1.0
  else
    let decay := Float.exp (-damping * frequency * t)
    let oscillation := Float.cos (frequency * Float.sqrt (1.0 - damping * damping) * t)
    1.0 - decay * oscillation

-- ============================================================================
-- Utility functions
-- ============================================================================

/-- Apply an easing function to interpolate between two values. -/
@[inline]
def apply (ease : Float → Float) (start finish t : Float) : Float :=
  start + (finish - start) * ease t

/-- Clamp t to [0, 1] before applying easing. -/
@[inline]
def clampedApply (ease : Float → Float) (start finish t : Float) : Float :=
  apply ease start finish (Float.clamp t 0.0 1.0)

/-- Reverse an easing function (ease out becomes ease in). -/
@[inline]
def reverse (ease : Float → Float) (t : Float) : Float :=
  1.0 - ease (1.0 - t)

/-- Mirror an easing function to create in-out effect. -/
@[inline]
def mirror (ease : Float → Float) (t : Float) : Float :=
  if t < 0.5 then
    ease (2.0 * t) / 2.0
  else
    1.0 - ease (2.0 - 2.0 * t) / 2.0

end Easing

-- ============================================================================
-- SmoothDamp (Unity-style smooth following)
-- ============================================================================

/-- State for SmoothDamp interpolation.
    Tracks current value and velocity for smooth following. -/
structure SmoothDampState where
  current : Float
  velocity : Float
  deriving Repr, Inhabited

namespace SmoothDamp

/-- Smoothly interpolate a value toward a target using a critically damped spring.
    This is the Unity-style SmoothDamp algorithm.

    Parameters:
    - state: Current value and velocity
    - target: Target value to approach
    - smoothTime: Approximate time to reach target (seconds)
    - deltaTime: Time since last update (seconds)
    - maxSpeed: Maximum speed (optional, use Float.infinity for unlimited)

    Returns: (newValue, newState) -/
def step (state : SmoothDampState) (target : Float) (smoothTime : Float)
    (deltaTime : Float) (maxSpeed : Float := Float.infinity) : Float × SmoothDampState :=
  -- Clamp smoothTime to avoid division issues
  let smoothTime' := Float.max 0.0001 smoothTime
  let omega := 2.0 / smoothTime'

  let x := omega * deltaTime
  -- Approximate exp(-x) using Padé approximant for stability
  let exp := 1.0 / (1.0 + x + 0.48 * x * x + 0.235 * x * x * x)

  let change := state.current - target
  let originalTarget := target

  -- Clamp maximum change
  let maxChange := maxSpeed * smoothTime'
  let change' := Float.clamp change (-maxChange) maxChange
  let target' := state.current - change'

  let temp := (state.velocity + omega * change') * deltaTime
  let newVelocity := (state.velocity - omega * temp) * exp
  let newValue := target' + (change' + temp) * exp

  -- Prevent overshooting
  if (originalTarget - state.current > 0.0) == (newValue > originalTarget) then
    (originalTarget, { current := originalTarget, velocity := 0.0 })
  else
    (newValue, { current := newValue, velocity := newVelocity })

/-- Create initial state at a given value. -/
def init (value : Float) : SmoothDampState :=
  { current := value, velocity := 0.0 }

/-- Convenience function for one-shot smoothing (without state tracking).
    Less efficient than using state, but simpler for quick usage. -/
def value (current target smoothTime deltaTime : Float)
    (maxSpeed : Float := Float.infinity) : Float :=
  let state := { current := current, velocity := 0.0 : SmoothDampState }
  (step state target smoothTime deltaTime maxSpeed).1

end SmoothDamp

/-- SmoothDamp state for Vec2. -/
structure SmoothDampState2 where
  current : Vec2
  velocity : Vec2
  deriving Repr, Inhabited

namespace SmoothDampState2

/-- Create initial state at a given position. -/
def init (value : Vec2) : SmoothDampState2 :=
  { current := value, velocity := Vec2.zero }

/-- Step the smooth damp for Vec2. -/
def step (state : SmoothDampState2) (target : Vec2) (smoothTime : Float)
    (deltaTime : Float) (maxSpeed : Float := Float.infinity) : Vec2 × SmoothDampState2 :=
  let (newX, stateX) := SmoothDamp.step
    { current := state.current.x, velocity := state.velocity.x }
    target.x smoothTime deltaTime maxSpeed
  let (newY, stateY) := SmoothDamp.step
    { current := state.current.y, velocity := state.velocity.y }
    target.y smoothTime deltaTime maxSpeed
  let newValue := Vec2.mk newX newY
  let newVelocity := Vec2.mk stateX.velocity stateY.velocity
  (newValue, { current := newValue, velocity := newVelocity })

end SmoothDampState2

/-- SmoothDamp state for Vec3. -/
structure SmoothDampState3 where
  current : Vec3
  velocity : Vec3
  deriving Repr, Inhabited

namespace SmoothDampState3

/-- Create initial state at a given position. -/
def init (value : Vec3) : SmoothDampState3 :=
  { current := value, velocity := Vec3.zero }

/-- Step the smooth damp for Vec3. -/
def step (state : SmoothDampState3) (target : Vec3) (smoothTime : Float)
    (deltaTime : Float) (maxSpeed : Float := Float.infinity) : Vec3 × SmoothDampState3 :=
  let (newX, stateX) := SmoothDamp.step
    { current := state.current.x, velocity := state.velocity.x }
    target.x smoothTime deltaTime maxSpeed
  let (newY, stateY) := SmoothDamp.step
    { current := state.current.y, velocity := state.velocity.y }
    target.y smoothTime deltaTime maxSpeed
  let (newZ, stateZ) := SmoothDamp.step
    { current := state.current.z, velocity := state.velocity.z }
    target.z smoothTime deltaTime maxSpeed
  let newValue := Vec3.mk newX newY newZ
  let newVelocity := Vec3.mk stateX.velocity stateY.velocity stateZ.velocity
  (newValue, { current := newValue, velocity := newVelocity })

end SmoothDampState3

-- ============================================================================
-- Hermite Interpolation
-- ============================================================================

namespace Hermite

/-- Cubic Hermite interpolation between two points with tangents.
    Given points p0, p1 and tangents m0, m1, interpolates at parameter t ∈ [0, 1].

    H(t) = (2t³ - 3t² + 1)p0 + (t³ - 2t² + t)m0 + (-2t³ + 3t²)p1 + (t³ - t²)m1 -/
def interpolate (p0 m0 p1 m1 t : Float) : Float :=
  let t2 := t * t
  let t3 := t2 * t
  let h00 := 2.0 * t3 - 3.0 * t2 + 1.0  -- Basis for p0
  let h10 := t3 - 2.0 * t2 + t          -- Basis for m0
  let h01 := -2.0 * t3 + 3.0 * t2       -- Basis for p1
  let h11 := t3 - t2                     -- Basis for m1
  h00 * p0 + h10 * m0 + h01 * p1 + h11 * m1

/-- Cubic Hermite interpolation for Vec2. -/
def interpolateVec2 (p0 m0 p1 m1 : Vec2) (t : Float) : Vec2 :=
  Vec2.mk
    (interpolate p0.x m0.x p1.x m1.x t)
    (interpolate p0.y m0.y p1.y m1.y t)

/-- Cubic Hermite interpolation for Vec3. -/
def interpolateVec3 (p0 m0 p1 m1 : Vec3) (t : Float) : Vec3 :=
  Vec3.mk
    (interpolate p0.x m0.x p1.x m1.x t)
    (interpolate p0.y m0.y p1.y m1.y t)
    (interpolate p0.z m0.z p1.z m1.z t)

/-- Compute the derivative of cubic Hermite at parameter t.
    H'(t) = (6t² - 6t)p0 + (3t² - 4t + 1)m0 + (-6t² + 6t)p1 + (3t² - 2t)m1 -/
def derivative (p0 m0 p1 m1 t : Float) : Float :=
  let t2 := t * t
  let h00' := 6.0 * t2 - 6.0 * t
  let h10' := 3.0 * t2 - 4.0 * t + 1.0
  let h01' := -6.0 * t2 + 6.0 * t
  let h11' := 3.0 * t2 - 2.0 * t
  h00' * p0 + h10' * m0 + h01' * p1 + h11' * m1

/-- Compute the derivative for Vec2. -/
def derivativeVec2 (p0 m0 p1 m1 : Vec2) (t : Float) : Vec2 :=
  Vec2.mk
    (derivative p0.x m0.x p1.x m1.x t)
    (derivative p0.y m0.y p1.y m1.y t)

/-- Compute the derivative for Vec3. -/
def derivativeVec3 (p0 m0 p1 m1 : Vec3) (t : Float) : Vec3 :=
  Vec3.mk
    (derivative p0.x m0.x p1.x m1.x t)
    (derivative p0.y m0.y p1.y m1.y t)
    (derivative p0.z m0.z p1.z m1.z t)

/-- Cardinal spline interpolation (Catmull-Rom is a special case with tension=0).
    Automatically computes tangents from neighboring points.
    tension: 0 = Catmull-Rom, 0.5 = tighter curve, 1 = linear -/
def cardinal (p0 p1 p2 p3 t tension : Float) : Float :=
  let m1 := (1.0 - tension) * (p2 - p0) / 2.0
  let m2 := (1.0 - tension) * (p3 - p1) / 2.0
  interpolate p1 m1 p2 m2 t

/-- Cardinal spline for Vec2. -/
def cardinalVec2 (p0 p1 p2 p3 : Vec2) (t tension : Float) : Vec2 :=
  let m1 := (p2 - p0).scale ((1.0 - tension) / 2.0)
  let m2 := (p3 - p1).scale ((1.0 - tension) / 2.0)
  interpolateVec2 p1 m1 p2 m2 t

/-- Cardinal spline for Vec3. -/
def cardinalVec3 (p0 p1 p2 p3 : Vec3) (t tension : Float) : Vec3 :=
  let m1 := (p2 - p0).scale ((1.0 - tension) / 2.0)
  let m2 := (p3 - p1).scale ((1.0 - tension) / 2.0)
  interpolateVec3 p1 m1 p2 m2 t

end Hermite

-- ============================================================================
-- Animation Curves (Keyframe-based)
-- ============================================================================

/-- Tangent mode for animation curve keyframes. -/
inductive TangentMode
  | linear      -- Linear interpolation to/from this key
  | constant    -- Step function (hold value until next key)
  | smooth      -- Auto-compute smooth tangent
  | custom      -- User-specified tangent value
  deriving Repr, BEq, Inhabited

/-- A keyframe in an animation curve. -/
structure Keyframe where
  time : Float           -- Time position
  value : Float          -- Value at this time
  inTangent : Float      -- Incoming tangent (for Hermite)
  outTangent : Float     -- Outgoing tangent (for Hermite)
  inMode : TangentMode   -- How to handle incoming interpolation
  outMode : TangentMode  -- How to handle outgoing interpolation
  deriving Repr, Inhabited

namespace Keyframe

/-- Create a keyframe with automatic smooth tangents. -/
def smooth (time value : Float) : Keyframe :=
  { time := time
  , value := value
  , inTangent := 0.0
  , outTangent := 0.0
  , inMode := TangentMode.smooth
  , outMode := TangentMode.smooth }

/-- Create a keyframe with linear tangents. -/
def linear (time value : Float) : Keyframe :=
  { time := time
  , value := value
  , inTangent := 0.0
  , outTangent := 0.0
  , inMode := TangentMode.linear
  , outMode := TangentMode.linear }

/-- Create a keyframe with constant (step) interpolation. -/
def constant (time value : Float) : Keyframe :=
  { time := time
  , value := value
  , inTangent := 0.0
  , outTangent := 0.0
  , inMode := TangentMode.constant
  , outMode := TangentMode.constant }

/-- Create a keyframe with custom tangents. -/
def custom (time value inTan outTan : Float) : Keyframe :=
  { time := time
  , value := value
  , inTangent := inTan
  , outTangent := outTan
  , inMode := TangentMode.custom
  , outMode := TangentMode.custom }

end Keyframe

/-- Wrap mode for animation curves. -/
inductive WrapMode
  | clamp      -- Clamp to first/last keyframe value
  | loop       -- Loop back to start
  | pingPong   -- Alternate direction
  deriving Repr, BEq, Inhabited

/-- An animation curve defined by keyframes with Hermite interpolation. -/
structure AnimationCurve where
  keyframes : Array Keyframe
  preWrapMode : WrapMode
  postWrapMode : WrapMode
  deriving Repr, Inhabited

namespace AnimationCurve

/-- Create an empty animation curve. -/
def empty : AnimationCurve :=
  { keyframes := #[]
  , preWrapMode := WrapMode.clamp
  , postWrapMode := WrapMode.clamp }

/-- Create a curve from an array of (time, value) pairs with smooth tangents. -/
def fromPoints (points : Array (Float × Float)) : AnimationCurve :=
  let keys := points.map fun (t, v) => Keyframe.smooth t v
  let sorted := keys.qsort (fun a b => a.time < b.time)
  { keyframes := computeSmoothTangents sorted
  , preWrapMode := WrapMode.clamp
  , postWrapMode := WrapMode.clamp }
where
  computeSmoothTangents (keys : Array Keyframe) : Array Keyframe :=
    if keys.size < 2 then keys
    else Id.run do
      let mut result := keys
      for i in [:keys.size] do
        let key := keys[i]!
        if key.inMode == TangentMode.smooth || key.outMode == TangentMode.smooth then
          let tangent :=
            if i == 0 then
              -- First key: use forward difference
              let next := keys[1]!
              let dt := next.time - key.time
              if dt > Float.epsilon then (next.value - key.value) / dt else 0.0
            else if i == keys.size - 1 then
              -- Last key: use backward difference
              let prev := keys[i - 1]!
              let dt := key.time - prev.time
              if dt > Float.epsilon then (key.value - prev.value) / dt else 0.0
            else
              -- Middle key: use central difference
              let prev := keys[i - 1]!
              let next := keys[i + 1]!
              let dt := next.time - prev.time
              if dt > Float.epsilon then (next.value - prev.value) / dt else 0.0
          let newKey := { key with
            inTangent := if key.inMode == TangentMode.smooth then tangent else key.inTangent
            outTangent := if key.outMode == TangentMode.smooth then tangent else key.outTangent }
          result := result.set! i newKey
      return result

/-- Create a linear curve from start to end. -/
def linear (startTime startValue endTime endValue : Float) : AnimationCurve :=
  { keyframes := #[Keyframe.linear startTime startValue, Keyframe.linear endTime endValue]
  , preWrapMode := WrapMode.clamp
  , postWrapMode := WrapMode.clamp }

/-- Create a constant curve (always returns the same value). -/
def constant (value : Float) : AnimationCurve :=
  { keyframes := #[Keyframe.constant 0.0 value]
  , preWrapMode := WrapMode.clamp
  , postWrapMode := WrapMode.clamp }

/-- Create an eased curve using an easing function. -/
def eased (startTime startValue endTime endValue : Float)
    (ease : Float → Float) (samples : Nat := 10) : AnimationCurve :=
  let points := Id.run do
    let mut arr : Array (Float × Float) := #[]
    for i in [:samples + 1] do
      let t := i.toFloat / samples.toFloat
      let time := Float.lerp startTime endTime t
      let value := Float.lerp startValue endValue (ease t)
      arr := arr.push (time, value)
    return arr
  fromPoints points

/-- Get the duration of the curve (time of last keyframe - time of first keyframe). -/
def duration (curve : AnimationCurve) : Float :=
  if curve.keyframes.isEmpty then 0.0
  else
    let first := curve.keyframes[0]!
    let last := curve.keyframes[curve.keyframes.size - 1]!
    last.time - first.time

/-- Get the start time of the curve. -/
def startTime (curve : AnimationCurve) : Float :=
  if curve.keyframes.isEmpty then 0.0
  else curve.keyframes[0]!.time

/-- Get the end time of the curve. -/
def endTime (curve : AnimationCurve) : Float :=
  if curve.keyframes.isEmpty then 0.0
  else curve.keyframes[curve.keyframes.size - 1]!.time

/-- Apply wrap mode to a time value. -/
private def wrapTime (curve : AnimationCurve) (time : Float) : Float :=
  if curve.keyframes.isEmpty then time
  else
    let start := curve.startTime
    let finish := curve.endTime
    let dur := finish - start

    if dur <= Float.epsilon then start
    else if time < start then
      match curve.preWrapMode with
      | WrapMode.clamp => start
      | WrapMode.loop =>
        let cycles := Float.ceil ((start - time) / dur)
        time + cycles * dur
      | WrapMode.pingPong =>
        let offset := start - time
        let cycles := Float.floor (offset / dur)
        let remainder := offset - cycles * dur
        if cycles.toUInt64 % 2 == 0 then start + remainder
        else finish - remainder
    else if time > finish then
      match curve.postWrapMode with
      | WrapMode.clamp => finish
      | WrapMode.loop =>
        let cycles := Float.floor ((time - start) / dur)
        time - cycles * dur
      | WrapMode.pingPong =>
        let offset := time - start
        let cycles := Float.floor (offset / dur)
        let remainder := offset - cycles * dur
        if cycles.toUInt64 % 2 == 0 then start + remainder
        else finish - remainder
    else time

/-- Find the keyframe segment containing a given time.
    Returns (index of key before, local t within segment). -/
private def findSegment (curve : AnimationCurve) (time : Float) : Nat × Float :=
  if curve.keyframes.size < 2 then (0, 0.0)
  else Id.run do
    -- Binary search for the segment
    let mut lo := 0
    let mut hi := curve.keyframes.size - 1
    while lo < hi - 1 do
      let mid := (lo + hi) / 2
      if curve.keyframes[mid]!.time <= time then
        lo := mid
      else
        hi := mid

    let k0 := curve.keyframes[lo]!
    let k1 := curve.keyframes[lo + 1]!
    let dt := k1.time - k0.time
    let localT := if dt > Float.epsilon then (time - k0.time) / dt else 0.0
    return (lo, Float.clamp localT 0.0 1.0)

/-- Evaluate the animation curve at a given time. -/
def evaluate (curve : AnimationCurve) (time : Float) : Float :=
  if curve.keyframes.isEmpty then 0.0
  else if curve.keyframes.size == 1 then curve.keyframes[0]!.value
  else
    let wrappedTime := curve.wrapTime time
    let (idx, localT) := curve.findSegment wrappedTime
    let k0 := curve.keyframes[idx]!
    let k1 := curve.keyframes[idx + 1]!

    -- Handle constant mode
    if k0.outMode == TangentMode.constant then k0.value
    else if k1.inMode == TangentMode.constant then k0.value
    -- Handle linear mode
    else if k0.outMode == TangentMode.linear && k1.inMode == TangentMode.linear then
      Float.lerp k0.value k1.value localT
    -- Default: Hermite interpolation
    else
      let dt := k1.time - k0.time
      -- Scale tangents by time interval for proper Hermite
      let m0 := k0.outTangent * dt
      let m1 := k1.inTangent * dt
      Hermite.interpolate k0.value m0 k1.value m1 localT

/-- Add a keyframe to the curve (maintains sorted order). -/
def addKey (curve : AnimationCurve) (key : Keyframe) : AnimationCurve :=
  let idx := Id.run do
    for i in [:curve.keyframes.size] do
      if curve.keyframes[i]!.time > key.time then return i
    return curve.keyframes.size
  -- Insert at position idx
  let before := curve.keyframes.toSubarray 0 idx |>.toArray
  let after := curve.keyframes.toSubarray idx curve.keyframes.size |>.toArray
  { curve with keyframes := before ++ #[key] ++ after }

/-- Remove the keyframe at a given index. -/
def removeKey (curve : AnimationCurve) (idx : Nat) : AnimationCurve :=
  if _ : idx < curve.keyframes.size then
    let before := curve.keyframes.toSubarray 0 idx |>.toArray
    let after := curve.keyframes.toSubarray (idx + 1) curve.keyframes.size |>.toArray
    { curve with keyframes := before ++ after }
  else curve

/-- Set the wrap modes. -/
def withWrapMode (curve : AnimationCurve) (pre post : WrapMode) : AnimationCurve :=
  { curve with preWrapMode := pre, postWrapMode := post }

/-- Sample the curve at regular intervals. -/
def sample (curve : AnimationCurve) (numSamples : Nat) : Array Float := Id.run do
  if numSamples < 2 then return #[curve.evaluate curve.startTime]
  let start := curve.startTime
  let finish := curve.endTime
  let mut result : Array Float := #[]
  for i in [:numSamples] do
    let t := Float.lerp start finish (i.toFloat / (numSamples - 1).toFloat)
    result := result.push (curve.evaluate t)
  return result

end AnimationCurve

end Linalg
