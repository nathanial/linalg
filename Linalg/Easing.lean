/-
  Easing functions for animation and interpolation.

  All easing functions take a parameter t in [0, 1] and return a value
  that represents the eased progression. The functions follow the naming
  convention: easeIn (slow start), easeOut (slow end), easeInOut (slow both).
-/

import Linalg.Core

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

end Linalg
