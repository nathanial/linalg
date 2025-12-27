/-
  Linalg Core Utilities
  Mathematical constants and Float helper functions.
-/

namespace Linalg

/-- Pi constant. -/
def Float.pi : Float := 3.14159265358979323846

/-- Two times pi. -/
def Float.twoPi : Float := 6.28318530717958647692

/-- Pi divided by two. -/
def Float.halfPi : Float := 1.57079632679489661923

/-- Degrees to radians conversion factor. -/
def Float.degToRad : Float := Float.pi / 180.0

/-- Radians to degrees conversion factor. -/
def Float.radToDeg : Float := 180.0 / Float.pi

/-- Small epsilon for floating point comparisons. -/
def Float.epsilon : Float := 0.00001

/-- Positive infinity. -/
def Float.infinity : Float := 1.0 / 0.0

/-- Negative infinity. -/
def Float.negInfinity : Float := -1.0 / 0.0

/-- Maximum of two floats. -/
@[inline]
def Float.max (a b : Float) : Float := if a >= b then a else b

/-- Minimum of two floats. -/
@[inline]
def Float.min (a b : Float) : Float := if a <= b then a else b

/-- Absolute value of a float. -/
@[inline]
def Float.abs' (x : Float) : Float := if x >= 0.0 then x else -x

/-- Clamp a float to a range. -/
@[inline]
def Float.clamp (x lo hi : Float) : Float :=
  Float.max lo (Float.min hi x)

/-- Linear interpolation between two values. -/
@[inline]
def Float.lerp (a b t : Float) : Float :=
  a + (b - a) * t

/-- Check if two floats are approximately equal. -/
@[inline]
def Float.approxEq (a b : Float) (eps : Float := Float.epsilon) : Bool :=
  Float.abs' (a - b) < eps

/-- Convert degrees to radians. -/
@[inline]
def Float.toRadians (degrees : Float) : Float :=
  degrees * Float.degToRad

/-- Convert radians to degrees. -/
@[inline]
def Float.toDegrees (radians : Float) : Float :=
  radians * Float.radToDeg

end Linalg
