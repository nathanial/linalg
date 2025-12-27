/-
  Curves and splines for animation, paths, and procedural generation.

  Includes:
  - Bezier curves (quadratic and cubic)
  - Catmull-Rom splines
  - Arc-length parameterization
-/

import Linalg.Core
import Linalg.Vec2
import Linalg.Vec3

namespace Linalg

-- ============================================================================
-- Quadratic Bezier Curve
-- ============================================================================

/-- A quadratic Bezier curve defined by 3 control points.
    B(t) = (1-t)²P₀ + 2(1-t)tP₁ + t²P₂ -/
structure Bezier2 (V : Type) where
  p0 : V  -- Start point
  p1 : V  -- Control point
  p2 : V  -- End point
  deriving Repr

namespace Bezier2

/-- Evaluate using Vec2. -/
def evalVec2 (b : Bezier2 Vec2) (t : Float) : Vec2 :=
  let t2 := t * t
  let mt := 1.0 - t
  let mt2 := mt * mt
  let c0 := mt2
  let c1 := 2.0 * mt * t
  let c2 := t2
  Vec2.mk
    (c0 * b.p0.x + c1 * b.p1.x + c2 * b.p2.x)
    (c0 * b.p0.y + c1 * b.p1.y + c2 * b.p2.y)

/-- Evaluate using Vec3. -/
def evalVec3 (b : Bezier2 Vec3) (t : Float) : Vec3 :=
  let t2 := t * t
  let mt := 1.0 - t
  let mt2 := mt * mt
  let c0 := mt2
  let c1 := 2.0 * mt * t
  let c2 := t2
  Vec3.mk
    (c0 * b.p0.x + c1 * b.p1.x + c2 * b.p2.x)
    (c0 * b.p0.y + c1 * b.p1.y + c2 * b.p2.y)
    (c0 * b.p0.z + c1 * b.p1.z + c2 * b.p2.z)

/-- Compute the derivative (tangent) at parameter t.
    B'(t) = 2(1-t)(P₁-P₀) + 2t(P₂-P₁) -/
def derivativeVec2 (b : Bezier2 Vec2) (t : Float) : Vec2 :=
  let mt := 1.0 - t
  let d0 := b.p1.sub b.p0
  let d1 := b.p2.sub b.p1
  (d0.scale (2.0 * mt)).add (d1.scale (2.0 * t))

/-- Compute the derivative (tangent) at parameter t for Vec3. -/
def derivativeVec3 (b : Bezier2 Vec3) (t : Float) : Vec3 :=
  let mt := 1.0 - t
  let d0 := b.p1.sub b.p0
  let d1 := b.p2.sub b.p1
  (d0.scale (2.0 * mt)).add (d1.scale (2.0 * t))

/-- Get the unit tangent at parameter t. -/
def tangentVec2 (b : Bezier2 Vec2) (t : Float) : Vec2 :=
  (b.derivativeVec2 t).normalize

/-- Get the unit tangent at parameter t for Vec3. -/
def tangentVec3 (b : Bezier2 Vec3) (t : Float) : Vec3 :=
  (b.derivativeVec3 t).normalize

/-- Split the curve at parameter t into two curves. -/
def splitVec2 (b : Bezier2 Vec2) (t : Float) : Bezier2 Vec2 × Bezier2 Vec2 :=
  -- De Casteljau algorithm
  let p01 := b.p0.lerp b.p1 t
  let p12 := b.p1.lerp b.p2 t
  let p012 := p01.lerp p12 t
  let left := Bezier2.mk b.p0 p01 p012
  let right := Bezier2.mk p012 p12 b.p2
  (left, right)

/-- Split the curve at parameter t into two curves for Vec3. -/
def splitVec3 (b : Bezier2 Vec3) (t : Float) : Bezier2 Vec3 × Bezier2 Vec3 :=
  let p01 := b.p0.lerp b.p1 t
  let p12 := b.p1.lerp b.p2 t
  let p012 := p01.lerp p12 t
  let left := Bezier2.mk b.p0 p01 p012
  let right := Bezier2.mk p012 p12 b.p2
  (left, right)

/-- Approximate the arc length using recursive subdivision. -/
partial def arcLengthVec2 (b : Bezier2 Vec2) (tolerance : Float := 0.001) : Float :=
  let chord := b.p0.distance b.p2
  let control := b.p0.distance b.p1 + b.p1.distance b.p2
  if control - chord < tolerance then
    (chord + control) / 2.0
  else
    let (left, right) := b.splitVec2 0.5
    left.arcLengthVec2 tolerance + right.arcLengthVec2 tolerance

/-- Approximate the arc length for Vec3. -/
partial def arcLengthVec3 (b : Bezier2 Vec3) (tolerance : Float := 0.001) : Float :=
  let chord := b.p0.distance b.p2
  let control := b.p0.distance b.p1 + b.p1.distance b.p2
  if control - chord < tolerance then
    (chord + control) / 2.0
  else
    let (left, right) := b.splitVec3 0.5
    left.arcLengthVec3 tolerance + right.arcLengthVec3 tolerance

end Bezier2

-- ============================================================================
-- Cubic Bezier Curve
-- ============================================================================

/-- A cubic Bezier curve defined by 4 control points.
    B(t) = (1-t)³P₀ + 3(1-t)²tP₁ + 3(1-t)t²P₂ + t³P₃ -/
structure Bezier3 (V : Type) where
  p0 : V  -- Start point
  p1 : V  -- First control point
  p2 : V  -- Second control point
  p3 : V  -- End point
  deriving Repr

namespace Bezier3

/-- Evaluate using Vec2. -/
def evalVec2 (b : Bezier3 Vec2) (t : Float) : Vec2 :=
  let t2 := t * t
  let t3 := t2 * t
  let mt := 1.0 - t
  let mt2 := mt * mt
  let mt3 := mt2 * mt
  let c0 := mt3
  let c1 := 3.0 * mt2 * t
  let c2 := 3.0 * mt * t2
  let c3 := t3
  Vec2.mk
    (c0 * b.p0.x + c1 * b.p1.x + c2 * b.p2.x + c3 * b.p3.x)
    (c0 * b.p0.y + c1 * b.p1.y + c2 * b.p2.y + c3 * b.p3.y)

/-- Evaluate using Vec3. -/
def evalVec3 (b : Bezier3 Vec3) (t : Float) : Vec3 :=
  let t2 := t * t
  let t3 := t2 * t
  let mt := 1.0 - t
  let mt2 := mt * mt
  let mt3 := mt2 * mt
  let c0 := mt3
  let c1 := 3.0 * mt2 * t
  let c2 := 3.0 * mt * t2
  let c3 := t3
  Vec3.mk
    (c0 * b.p0.x + c1 * b.p1.x + c2 * b.p2.x + c3 * b.p3.x)
    (c0 * b.p0.y + c1 * b.p1.y + c2 * b.p2.y + c3 * b.p3.y)
    (c0 * b.p0.z + c1 * b.p1.z + c2 * b.p2.z + c3 * b.p3.z)

/-- Compute the derivative (tangent) at parameter t.
    B'(t) = 3(1-t)²(P₁-P₀) + 6(1-t)t(P₂-P₁) + 3t²(P₃-P₂) -/
def derivativeVec2 (b : Bezier3 Vec2) (t : Float) : Vec2 :=
  let t2 := t * t
  let mt := 1.0 - t
  let mt2 := mt * mt
  let d0 := b.p1.sub b.p0
  let d1 := b.p2.sub b.p1
  let d2 := b.p3.sub b.p2
  let r0 := d0.scale (3.0 * mt2)
  let r1 := d1.scale (6.0 * mt * t)
  let r2 := d2.scale (3.0 * t2)
  r0.add (r1.add r2)

/-- Compute the derivative (tangent) at parameter t for Vec3. -/
def derivativeVec3 (b : Bezier3 Vec3) (t : Float) : Vec3 :=
  let t2 := t * t
  let mt := 1.0 - t
  let mt2 := mt * mt
  let d0 := b.p1.sub b.p0
  let d1 := b.p2.sub b.p1
  let d2 := b.p3.sub b.p2
  let r0 := d0.scale (3.0 * mt2)
  let r1 := d1.scale (6.0 * mt * t)
  let r2 := d2.scale (3.0 * t2)
  r0.add (r1.add r2)

/-- Get the unit tangent at parameter t. -/
def tangentVec2 (b : Bezier3 Vec2) (t : Float) : Vec2 :=
  (b.derivativeVec2 t).normalize

/-- Get the unit tangent at parameter t for Vec3. -/
def tangentVec3 (b : Bezier3 Vec3) (t : Float) : Vec3 :=
  (b.derivativeVec3 t).normalize

/-- Split the curve at parameter t using De Casteljau algorithm. -/
def splitVec2 (b : Bezier3 Vec2) (t : Float) : Bezier3 Vec2 × Bezier3 Vec2 :=
  let p01 := b.p0.lerp b.p1 t
  let p12 := b.p1.lerp b.p2 t
  let p23 := b.p2.lerp b.p3 t
  let p012 := p01.lerp p12 t
  let p123 := p12.lerp p23 t
  let p0123 := p012.lerp p123 t
  let left := Bezier3.mk b.p0 p01 p012 p0123
  let right := Bezier3.mk p0123 p123 p23 b.p3
  (left, right)

/-- Split the curve at parameter t for Vec3. -/
def splitVec3 (b : Bezier3 Vec3) (t : Float) : Bezier3 Vec3 × Bezier3 Vec3 :=
  let p01 := b.p0.lerp b.p1 t
  let p12 := b.p1.lerp b.p2 t
  let p23 := b.p2.lerp b.p3 t
  let p012 := p01.lerp p12 t
  let p123 := p12.lerp p23 t
  let p0123 := p012.lerp p123 t
  let left := Bezier3.mk b.p0 p01 p012 p0123
  let right := Bezier3.mk p0123 p123 p23 b.p3
  (left, right)

/-- Approximate the arc length using recursive subdivision. -/
partial def arcLengthVec2 (b : Bezier3 Vec2) (tolerance : Float := 0.001) : Float :=
  let chord := b.p0.distance b.p3
  let control := b.p0.distance b.p1 + b.p1.distance b.p2 + b.p2.distance b.p3
  if control - chord < tolerance then
    (chord + control) / 2.0
  else
    let (left, right) := b.splitVec2 0.5
    left.arcLengthVec2 tolerance + right.arcLengthVec2 tolerance

/-- Approximate the arc length for Vec3. -/
partial def arcLengthVec3 (b : Bezier3 Vec3) (tolerance : Float := 0.001) : Float :=
  let chord := b.p0.distance b.p3
  let control := b.p0.distance b.p1 + b.p1.distance b.p2 + b.p2.distance b.p3
  if control - chord < tolerance then
    (chord + control) / 2.0
  else
    let (left, right) := b.splitVec3 0.5
    left.arcLengthVec3 tolerance + right.arcLengthVec3 tolerance

/-- Compute the second derivative at parameter t. -/
def secondDerivativeVec2 (b : Bezier3 Vec2) (t : Float) : Vec2 :=
  let mt := 1.0 - t
  let d0 := (b.p2.sub (b.p1.scale 2.0)).add b.p0
  let d1 := (b.p3.sub (b.p2.scale 2.0)).add b.p1
  (d0.scale (6.0 * mt)).add (d1.scale (6.0 * t))

/-- Compute the second derivative for Vec3. -/
def secondDerivativeVec3 (b : Bezier3 Vec3) (t : Float) : Vec3 :=
  let mt := 1.0 - t
  let d0 := (b.p2.sub (b.p1.scale 2.0)).add b.p0
  let d1 := (b.p3.sub (b.p2.scale 2.0)).add b.p1
  (d0.scale (6.0 * mt)).add (d1.scale (6.0 * t))

/-- Compute the curvature at parameter t. -/
def curvatureVec2 (b : Bezier3 Vec2) (t : Float) : Float :=
  let d1 := b.derivativeVec2 t
  let d2 := b.secondDerivativeVec2 t
  let cross := d1.x * d2.y - d1.y * d2.x
  let lenCubed := Float.pow d1.length 3.0
  if lenCubed < Float.epsilon then 0.0
  else cross / lenCubed

end Bezier3

-- ============================================================================
-- Catmull-Rom Spline
-- ============================================================================

/-- Catmull-Rom spline segment defined by 4 control points.
    The curve passes through p1 and p2, with p0 and p3 influencing tangents. -/
structure CatmullRom (V : Type) where
  p0 : V  -- Previous point (for tangent calculation)
  p1 : V  -- Start point of this segment
  p2 : V  -- End point of this segment
  p3 : V  -- Next point (for tangent calculation)
  alpha : Float := 0.5  -- 0 = uniform, 0.5 = centripetal, 1 = chordal
  deriving Repr

namespace CatmullRom

/-- Evaluate the Catmull-Rom spline at parameter t ∈ [0, 1].
    Uses the standard Catmull-Rom formula with tension = 0.5. -/
def evalVec2 (s : CatmullRom Vec2) (t : Float) : Vec2 :=
  let t2 := t * t
  let t3 := t2 * t
  -- Catmull-Rom basis matrix coefficients
  let c0 := -0.5 * t3 + t2 - 0.5 * t
  let c1 := 1.5 * t3 - 2.5 * t2 + 1.0
  let c2 := -1.5 * t3 + 2.0 * t2 + 0.5 * t
  let c3 := 0.5 * t3 - 0.5 * t2
  Vec2.mk
    (c0 * s.p0.x + c1 * s.p1.x + c2 * s.p2.x + c3 * s.p3.x)
    (c0 * s.p0.y + c1 * s.p1.y + c2 * s.p2.y + c3 * s.p3.y)

/-- Evaluate for Vec3. -/
def evalVec3 (s : CatmullRom Vec3) (t : Float) : Vec3 :=
  let t2 := t * t
  let t3 := t2 * t
  let c0 := -0.5 * t3 + t2 - 0.5 * t
  let c1 := 1.5 * t3 - 2.5 * t2 + 1.0
  let c2 := -1.5 * t3 + 2.0 * t2 + 0.5 * t
  let c3 := 0.5 * t3 - 0.5 * t2
  Vec3.mk
    (c0 * s.p0.x + c1 * s.p1.x + c2 * s.p2.x + c3 * s.p3.x)
    (c0 * s.p0.y + c1 * s.p1.y + c2 * s.p2.y + c3 * s.p3.y)
    (c0 * s.p0.z + c1 * s.p1.z + c2 * s.p2.z + c3 * s.p3.z)

/-- Compute the derivative at parameter t. -/
def derivativeVec2 (s : CatmullRom Vec2) (t : Float) : Vec2 :=
  let t2 := t * t
  let c0 := -1.5 * t2 + 2.0 * t - 0.5
  let c1 := 4.5 * t2 - 5.0 * t
  let c2 := -4.5 * t2 + 4.0 * t + 0.5
  let c3 := 1.5 * t2 - t
  Vec2.mk
    (c0 * s.p0.x + c1 * s.p1.x + c2 * s.p2.x + c3 * s.p3.x)
    (c0 * s.p0.y + c1 * s.p1.y + c2 * s.p2.y + c3 * s.p3.y)

/-- Compute the derivative for Vec3. -/
def derivativeVec3 (s : CatmullRom Vec3) (t : Float) : Vec3 :=
  let t2 := t * t
  let c0 := -1.5 * t2 + 2.0 * t - 0.5
  let c1 := 4.5 * t2 - 5.0 * t
  let c2 := -4.5 * t2 + 4.0 * t + 0.5
  let c3 := 1.5 * t2 - t
  Vec3.mk
    (c0 * s.p0.x + c1 * s.p1.x + c2 * s.p2.x + c3 * s.p3.x)
    (c0 * s.p0.y + c1 * s.p1.y + c2 * s.p2.y + c3 * s.p3.y)
    (c0 * s.p0.z + c1 * s.p1.z + c2 * s.p2.z + c3 * s.p3.z)

/-- Get the unit tangent at parameter t. -/
def tangentVec2 (s : CatmullRom Vec2) (t : Float) : Vec2 :=
  (s.derivativeVec2 t).normalize

/-- Get the unit tangent for Vec3. -/
def tangentVec3 (s : CatmullRom Vec3) (t : Float) : Vec3 :=
  (s.derivativeVec3 t).normalize

end CatmullRom

-- ============================================================================
-- Spline Path (sequence of points with Catmull-Rom interpolation)
-- ============================================================================

/-- A path defined by a sequence of control points, interpolated with Catmull-Rom. -/
structure SplinePath2 where
  points : Array Vec2
  closed : Bool := false
  deriving Repr

/-- A path defined by a sequence of Vec3 control points. -/
structure SplinePath3 where
  points : Array Vec3
  closed : Bool := false
  deriving Repr

namespace SplinePath2

/-- Get the number of segments in the spline. -/
def segmentCount (sp : SplinePath2) : Nat :=
  if sp.points.size < 2 then 0
  else if sp.closed then sp.points.size
  else sp.points.size - 1

/-- Clamp index to valid range for open splines, or wrap for closed. -/
def wrapIndex (sp : SplinePath2) (i : Int) : Nat :=
  let n := sp.points.size
  if n == 0 then 0
  else if sp.closed then
    let m := i % (Int.ofNat n)
    if m < 0 then (m + Int.ofNat n).toNat else m.toNat
  else
    if i < 0 then 0
    else if i >= Int.ofNat n then n - 1
    else i.toNat

/-- Get the control points for a segment. -/
def getSegment (sp : SplinePath2) (segmentIndex : Nat) : Option (CatmullRom Vec2) :=
  let n := sp.points.size
  if n < 2 then none
  else if !sp.closed && segmentIndex >= n - 1 then none
  else if sp.closed && segmentIndex >= n then none
  else
    let i := Int.ofNat segmentIndex
    let p0 := sp.points.getD (sp.wrapIndex (i - 1)) Vec2.zero
    let p1 := sp.points.getD (sp.wrapIndex i) Vec2.zero
    let p2 := sp.points.getD (sp.wrapIndex (i + 1)) Vec2.zero
    let p3 := sp.points.getD (sp.wrapIndex (i + 2)) Vec2.zero
    some { p0 := p0, p1 := p1, p2 := p2, p3 := p3 }

/-- Evaluate the spline at a global parameter t ∈ [0, 1]. -/
def eval (sp : SplinePath2) (t : Float) : Vec2 :=
  let numSegs := sp.segmentCount
  if numSegs == 0 then
    sp.points.getD 0 Vec2.zero
  else
    let t' := Float.clamp t 0.0 1.0
    let segFloat := t' * numSegs.toFloat
    let segIndex := Float.min segFloat (numSegs.toFloat - 1.0) |> Float.floor |> Float.toUInt64 |> UInt64.toNat
    let localT := segFloat - segIndex.toFloat
    match sp.getSegment segIndex with
    | some seg => seg.evalVec2 localT
    | none => sp.points.getD 0 Vec2.zero

end SplinePath2

namespace SplinePath3

/-- Get the number of segments in the spline. -/
def segmentCount (sp : SplinePath3) : Nat :=
  if sp.points.size < 2 then 0
  else if sp.closed then sp.points.size
  else sp.points.size - 1

/-- Clamp index to valid range for open splines, or wrap for closed. -/
def wrapIndex (sp : SplinePath3) (i : Int) : Nat :=
  let n := sp.points.size
  if n == 0 then 0
  else if sp.closed then
    let m := i % (Int.ofNat n)
    if m < 0 then (m + Int.ofNat n).toNat else m.toNat
  else
    if i < 0 then 0
    else if i >= Int.ofNat n then n - 1
    else i.toNat

/-- Get the control points for a segment. -/
def getSegment (sp : SplinePath3) (segmentIndex : Nat) : Option (CatmullRom Vec3) :=
  let n := sp.points.size
  if n < 2 then none
  else if !sp.closed && segmentIndex >= n - 1 then none
  else if sp.closed && segmentIndex >= n then none
  else
    let i := Int.ofNat segmentIndex
    let p0 := sp.points.getD (sp.wrapIndex (i - 1)) Vec3.zero
    let p1 := sp.points.getD (sp.wrapIndex i) Vec3.zero
    let p2 := sp.points.getD (sp.wrapIndex (i + 1)) Vec3.zero
    let p3 := sp.points.getD (sp.wrapIndex (i + 2)) Vec3.zero
    some { p0 := p0, p1 := p1, p2 := p2, p3 := p3 }

/-- Evaluate the spline at a global parameter t ∈ [0, 1]. -/
def eval (sp : SplinePath3) (t : Float) : Vec3 :=
  let numSegs := sp.segmentCount
  if numSegs == 0 then
    sp.points.getD 0 Vec3.zero
  else
    let t' := Float.clamp t 0.0 1.0
    let segFloat := t' * numSegs.toFloat
    let segIndex := Float.min segFloat (numSegs.toFloat - 1.0) |> Float.floor |> Float.toUInt64 |> UInt64.toNat
    let localT := segFloat - segIndex.toFloat
    match sp.getSegment segIndex with
    | some seg => seg.evalVec3 localT
    | none => sp.points.getD 0 Vec3.zero

end SplinePath3

-- ============================================================================
-- Arc-Length Parameterization
-- ============================================================================

/-- A lookup table for arc-length parameterization.
    Maps arc-length s ∈ [0, totalLength] to parameter t ∈ [0, 1]. -/
structure ArcLengthTable where
  samples : Array (Float × Float)  -- (t, cumulativeLength)
  totalLength : Float
  deriving Repr

namespace ArcLengthTable

/-- Build an arc-length table from a curve evaluation function.
    sampleCount: number of samples (more = better accuracy) -/
def build (eval : Float → Vec2) (sampleCount : Nat := 100) : ArcLengthTable := Id.run do
  if sampleCount < 2 then
    return { samples := #[(0.0, 0.0), (1.0, 0.0)], totalLength := 0.0 }

  let dt := 1.0 / (sampleCount - 1).toFloat
  let mut samples : Array (Float × Float) := #[]
  let mut cumulativeLength := 0.0
  let mut prevPoint := eval 0.0

  for i in [:sampleCount] do
    let t := i.toFloat * dt
    let point := eval t
    if i > 0 then
      cumulativeLength := cumulativeLength + prevPoint.distance point
    samples := samples.push (t, cumulativeLength)
    prevPoint := point

  return { samples := samples, totalLength := cumulativeLength }

/-- Build for Vec3 curves. -/
def buildVec3 (eval : Float → Vec3) (sampleCount : Nat := 100) : ArcLengthTable := Id.run do
  if sampleCount < 2 then
    return { samples := #[(0.0, 0.0), (1.0, 0.0)], totalLength := 0.0 }

  let dt := 1.0 / (sampleCount - 1).toFloat
  let mut samples : Array (Float × Float) := #[]
  let mut cumulativeLength := 0.0
  let mut prevPoint := eval 0.0

  for i in [:sampleCount] do
    let t := i.toFloat * dt
    let point := eval t
    if i > 0 then
      cumulativeLength := cumulativeLength + prevPoint.distance point
    samples := samples.push (t, cumulativeLength)
    prevPoint := point

  return { samples := samples, totalLength := cumulativeLength }

/-- Convert arc-length s to parameter t using binary search. -/
def sToT (table : ArcLengthTable) (s : Float) : Float := Id.run do
  if table.samples.size < 2 || table.totalLength < Float.epsilon then
    return 0.0

  let s' := Float.clamp s 0.0 table.totalLength
  -- Binary search for the segment containing s'
  let mut lo := 0
  let mut hi := table.samples.size - 1

  while lo < hi do
    let mid := (lo + hi) / 2
    let (_, sampleS) := table.samples.getD mid (0.0, 0.0)
    if sampleS < s' then
      lo := mid + 1
    else
      hi := mid

  -- Linear interpolation within the segment
  if lo == 0 then
    let (t1, s1) := table.samples.getD 0 (0.0, 0.0)
    if s1 < Float.epsilon then return t1
    else return t1 * (s' / s1)
  else
    let (t0, s0) := table.samples.getD (lo - 1) (0.0, 0.0)
    let (t1, s1) := table.samples.getD lo (0.0, 0.0)
    let ds := s1 - s0
    if ds < Float.epsilon then return t0
    else
      let frac := (s' - s0) / ds
      return Float.lerp t0 t1 frac

/-- Convert normalized arc-length u ∈ [0, 1] to parameter t. -/
def uToT (table : ArcLengthTable) (u : Float) : Float :=
  table.sToT (u * table.totalLength)

/-- Evaluate the curve at arc-length s. -/
def evalAtLength (table : ArcLengthTable) (eval : Float → Vec2) (s : Float) : Vec2 :=
  eval (table.sToT s)

/-- Evaluate the curve at normalized arc-length u ∈ [0, 1]. -/
def evalAtU (table : ArcLengthTable) (eval : Float → Vec2) (u : Float) : Vec2 :=
  eval (table.uToT u)

end ArcLengthTable

-- ============================================================================
-- B-Spline (Basis Spline)
-- ============================================================================

/-- A B-spline curve defined by control points and knot vector.
    B-splines provide local control - moving one control point only affects
    a local portion of the curve.

    For a degree-k B-spline with n+1 control points, we need n+k+2 knots.
    Common cases:
    - Cubic B-spline (degree 3): smooth curves, C² continuity
    - Quadratic B-spline (degree 2): parabolic segments, C¹ continuity -/
structure BSpline (V : Type) where
  controlPoints : Array V
  knots : Array Float
  degree : Nat
  deriving Repr

namespace BSpline

/-- Create a uniform B-spline (equally spaced knots) of given degree.
    Automatically generates an open/clamped knot vector. -/
def uniform {V : Type} (points : Array V) (degree : Nat := 3) : BSpline V :=
  let n := points.size
  let k := degree
  if n < k + 1 then
    { controlPoints := points, knots := #[], degree := k }
  else
    -- Generate clamped uniform knot vector
    -- First k+1 knots are 0, last k+1 knots are 1, middle knots are uniform
    let numKnots := n + k + 1
    let numInternalKnots := numKnots - 2 * (k + 1)
    let knots := Id.run do
      let mut arr : Array Float := #[]
      -- First k+1 knots = 0
      for _ in [:k+1] do
        arr := arr.push 0.0
      -- Internal knots
      for i in [:numInternalKnots] do
        arr := arr.push ((i + 1).toFloat / (numInternalKnots + 1).toFloat)
      -- Last k+1 knots = 1
      for _ in [:k+1] do
        arr := arr.push 1.0
      return arr
    { controlPoints := points, knots := knots, degree := k }

/-- Cox-de Boor recursive basis function.
    N_{i,k}(t) where i is the knot index and k is the degree. -/
partial def basisFunction (knots : Array Float) (i : Nat) (k : Nat) (t : Float) : Float :=
  if k == 0 then
    -- Base case: degree 0
    let ti := knots.getD i 0.0
    let ti1 := knots.getD (i + 1) 0.0
    if t >= ti && t < ti1 then 1.0
    else if t == ti1 && ti1 == knots.getD (knots.size - 1) 1.0 then 1.0  -- Handle endpoint
    else 0.0
  else
    -- Recursive case
    let ti := knots.getD i 0.0
    let tik := knots.getD (i + k) 0.0
    let ti1 := knots.getD (i + 1) 0.0
    let tik1 := knots.getD (i + k + 1) 0.0

    let left :=
      if Float.abs (tik - ti) < Float.epsilon then 0.0
      else (t - ti) / (tik - ti) * basisFunction knots i (k - 1) t

    let right :=
      if Float.abs (tik1 - ti1) < Float.epsilon then 0.0
      else (tik1 - t) / (tik1 - ti1) * basisFunction knots (i + 1) (k - 1) t

    left + right

/-- Evaluate the B-spline at parameter t ∈ [0, 1] for Vec2. -/
def evalVec2 (b : BSpline Vec2) (t : Float) : Vec2 :=
  if b.controlPoints.isEmpty || b.knots.size < b.degree + 2 then Vec2.zero
  else Id.run do
    let t' := Float.clamp t 0.0 1.0
    let mut result := Vec2.zero
    for i in [:b.controlPoints.size] do
      let basis := basisFunction b.knots i b.degree t'
      let p := b.controlPoints.getD i Vec2.zero
      result := result.add (p.scale basis)
    return result

/-- Evaluate the B-spline at parameter t ∈ [0, 1] for Vec3. -/
def evalVec3 (b : BSpline Vec3) (t : Float) : Vec3 :=
  if b.controlPoints.isEmpty || b.knots.size < b.degree + 2 then Vec3.zero
  else Id.run do
    let t' := Float.clamp t 0.0 1.0
    let mut result := Vec3.zero
    for i in [:b.controlPoints.size] do
      let basis := basisFunction b.knots i b.degree t'
      let p := b.controlPoints.getD i Vec3.zero
      result := result.add (p.scale basis)
    return result

/-- Evaluate using de Boor's algorithm (more numerically stable).
    Returns the point on the curve at parameter t. -/
def deBoorVec2 (b : BSpline Vec2) (t : Float) : Vec2 :=
  if b.controlPoints.isEmpty || b.knots.size < b.degree + 2 then Vec2.zero
  else
    let t' := Float.clamp t 0.0 1.0
    let k := b.degree
    let n := b.controlPoints.size

    -- Find knot span (index of the knot interval containing t)
    let findSpan : Nat := Id.run do
      if t' >= b.knots.getD (n) 1.0 then return n - 1
      let mut low := k
      let mut high := n
      let mut mid := (low + high) / 2
      while t' < b.knots.getD mid 0.0 || t' >= b.knots.getD (mid + 1) 1.0 do
        if t' < b.knots.getD mid 0.0 then
          high := mid
        else
          low := mid
        mid := (low + high) / 2
        if low >= high - 1 then break
      return mid

    let span := findSpan

    -- de Boor's algorithm
    let points := Id.run do
      let mut d : Array Vec2 := #[]
      for j in [:k+1] do
        d := d.push (b.controlPoints.getD (span - k + j) Vec2.zero)

      for r in [1:k+1] do
        for j' in [:k-r+1] do
          let j := k - r - j'  -- Reverse iteration
          let i := span - k + j + r
          let ti := b.knots.getD i 0.0
          let tikr := b.knots.getD (i + k - r + 1) 0.0
          let alpha :=
            if Float.abs (tikr - ti) < Float.epsilon then 0.0
            else (t' - ti) / (tikr - ti)
          let dj := d.getD j Vec2.zero
          let dj1 := d.getD (j + 1) Vec2.zero
          d := d.set! j (dj.lerp dj1 alpha)
      return d

    points.getD 0 Vec2.zero

/-- Evaluate using de Boor's algorithm for Vec3. -/
def deBoorVec3 (b : BSpline Vec3) (t : Float) : Vec3 :=
  if b.controlPoints.isEmpty || b.knots.size < b.degree + 2 then Vec3.zero
  else
    let t' := Float.clamp t 0.0 1.0
    let k := b.degree
    let n := b.controlPoints.size

    let findSpan : Nat := Id.run do
      if t' >= b.knots.getD (n) 1.0 then return n - 1
      let mut low := k
      let mut high := n
      let mut mid := (low + high) / 2
      while t' < b.knots.getD mid 0.0 || t' >= b.knots.getD (mid + 1) 1.0 do
        if t' < b.knots.getD mid 0.0 then
          high := mid
        else
          low := mid
        mid := (low + high) / 2
        if low >= high - 1 then break
      return mid

    let span := findSpan

    let points := Id.run do
      let mut d : Array Vec3 := #[]
      for j in [:k+1] do
        d := d.push (b.controlPoints.getD (span - k + j) Vec3.zero)

      for r in [1:k+1] do
        for j' in [:k-r+1] do
          let j := k - r - j'
          let i := span - k + j + r
          let ti := b.knots.getD i 0.0
          let tikr := b.knots.getD (i + k - r + 1) 0.0
          let alpha :=
            if Float.abs (tikr - ti) < Float.epsilon then 0.0
            else (t' - ti) / (tikr - ti)
          let dj := d.getD j Vec3.zero
          let dj1 := d.getD (j + 1) Vec3.zero
          d := d.set! j (dj.lerp dj1 alpha)
      return d

    points.getD 0 Vec3.zero

/-- Sample the B-spline at n equally spaced parameters. -/
def sampleVec2 (b : BSpline Vec2) (numSamples : Nat) : Array Vec2 := Id.run do
  if numSamples < 2 then return #[b.evalVec2 0.0]
  let mut result : Array Vec2 := #[]
  for i in [:numSamples] do
    let t := i.toFloat / (numSamples - 1).toFloat
    result := result.push (b.evalVec2 t)
  return result

/-- Sample the B-spline for Vec3. -/
def sampleVec3 (b : BSpline Vec3) (numSamples : Nat) : Array Vec3 := Id.run do
  if numSamples < 2 then return #[b.evalVec3 0.0]
  let mut result : Array Vec3 := #[]
  for i in [:numSamples] do
    let t := i.toFloat / (numSamples - 1).toFloat
    result := result.push (b.evalVec3 t)
  return result

/-- Check if the B-spline is valid (has enough control points for degree). -/
def isValid {V : Type} (b : BSpline V) : Bool :=
  b.controlPoints.size >= b.degree + 1 &&
  b.knots.size == b.controlPoints.size + b.degree + 1

end BSpline

-- ============================================================================
-- Bezier Patch (Bicubic Surface)
-- ============================================================================

/-- A bicubic Bezier patch - a smooth surface defined by a 4×4 grid of control points.
    The surface S(u,v) is computed as the tensor product of cubic Bezier curves.

    Control points are stored row-major:
    p00 p01 p02 p03
    p10 p11 p12 p13
    p20 p21 p22 p23
    p30 p31 p32 p33 -/
structure BezierPatch where
  controlPoints : Array Vec3  -- 16 control points (4×4 grid)
  deriving Repr, Inhabited

namespace BezierPatch

/-- Create a Bezier patch from a 4×4 grid of control points. -/
def fromGrid (grid : Array (Array Vec3)) : BezierPatch :=
  let points := Id.run do
    let mut arr : Array Vec3 := #[]
    for row in grid do
      for p in row do
        arr := arr.push p
    return arr
  { controlPoints := points }

/-- Create a flat patch in the XY plane centered at origin. -/
def flat (width height : Float) : BezierPatch :=
  let hw := width / 2.0
  let hh := height / 2.0
  let points := Id.run do
    let mut arr : Array Vec3 := #[]
    for j in [:4] do
      for i in [:4] do
        let x := -hw + (i.toFloat / 3.0) * width
        let y := -hh + (j.toFloat / 3.0) * height
        arr := arr.push (Vec3.mk x y 0.0)
    return arr
  { controlPoints := points }

/-- Get control point at grid position (row, col). -/
def getPoint (p : BezierPatch) (row col : Nat) : Vec3 :=
  let idx := row * 4 + col
  p.controlPoints.getD idx Vec3.zero

/-- Set control point at grid position (row, col). -/
def setPoint (p : BezierPatch) (row col : Nat) (v : Vec3) : BezierPatch :=
  let idx := row * 4 + col
  if idx < p.controlPoints.size then
    { controlPoints := p.controlPoints.set! idx v }
  else p

/-- Cubic Bernstein polynomial: B_i(t) = C(3,i) * t^i * (1-t)^(3-i) -/
private def bernstein3 (i : Nat) (t : Float) : Float :=
  let mt := 1.0 - t
  match i with
  | 0 => mt * mt * mt           -- (1-t)³
  | 1 => 3.0 * mt * mt * t      -- 3(1-t)²t
  | 2 => 3.0 * mt * t * t       -- 3(1-t)t²
  | 3 => t * t * t              -- t³
  | _ => 0.0

/-- Evaluate the patch at parameters (u, v) ∈ [0,1] × [0,1].
    Uses the tensor product formula:
    S(u,v) = Σᵢ Σⱼ B_i(u) * B_j(v) * P_{i,j} -/
def eval (p : BezierPatch) (u v : Float) : Vec3 := Id.run do
  let u' := Float.clamp u 0.0 1.0
  let v' := Float.clamp v 0.0 1.0
  let mut result := Vec3.zero
  for j in [:4] do
    for i in [:4] do
      let basis := bernstein3 i u' * bernstein3 j v'
      let cp := p.getPoint j i
      result := result.add (cp.scale basis)
  return result

/-- Compute the partial derivative ∂S/∂u at (u, v). -/
def derivativeU (p : BezierPatch) (u v : Float) : Vec3 := Id.run do
  let u' := Float.clamp u 0.0 1.0
  let v' := Float.clamp v 0.0 1.0
  -- Derivative of Bernstein polynomials
  let dB (i : Nat) (t : Float) : Float :=
    let mt := 1.0 - t
    match i with
    | 0 => -3.0 * mt * mt
    | 1 => 3.0 * mt * mt - 6.0 * mt * t
    | 2 => 6.0 * mt * t - 3.0 * t * t
    | 3 => 3.0 * t * t
    | _ => 0.0
  let mut result := Vec3.zero
  for j in [:4] do
    for i in [:4] do
      let basis := dB i u' * bernstein3 j v'
      let cp := p.getPoint j i
      result := result.add (cp.scale basis)
  return result

/-- Compute the partial derivative ∂S/∂v at (u, v). -/
def derivativeV (p : BezierPatch) (u v : Float) : Vec3 := Id.run do
  let u' := Float.clamp u 0.0 1.0
  let v' := Float.clamp v 0.0 1.0
  let dB (i : Nat) (t : Float) : Float :=
    let mt := 1.0 - t
    match i with
    | 0 => -3.0 * mt * mt
    | 1 => 3.0 * mt * mt - 6.0 * mt * t
    | 2 => 6.0 * mt * t - 3.0 * t * t
    | 3 => 3.0 * t * t
    | _ => 0.0
  let mut result := Vec3.zero
  for j in [:4] do
    for i in [:4] do
      let basis := bernstein3 i u' * dB j v'
      let cp := p.getPoint j i
      result := result.add (cp.scale basis)
  return result

/-- Compute the surface normal at (u, v). -/
def normal (p : BezierPatch) (u v : Float) : Vec3 :=
  let du := p.derivativeU u v
  let dv := p.derivativeV u v
  (du.cross dv).normalize

/-- Sample the patch as a grid of points (for rendering).
    Returns (rows × cols) points. -/
def sample (p : BezierPatch) (rows cols : Nat) : Array (Array Vec3) := Id.run do
  let r := if rows < 2 then 2 else rows
  let c := if cols < 2 then 2 else cols
  let mut result : Array (Array Vec3) := #[]
  for j in [:r] do
    let v := j.toFloat / (r - 1).toFloat
    let mut row : Array Vec3 := #[]
    for i in [:c] do
      let u := i.toFloat / (c - 1).toFloat
      row := row.push (p.eval u v)
    result := result.push row
  return result

/-- Sample the patch as a flat array of points. -/
def sampleFlat (p : BezierPatch) (rows cols : Nat) : Array Vec3 := Id.run do
  let r := if rows < 2 then 2 else rows
  let c := if cols < 2 then 2 else cols
  let mut result : Array Vec3 := #[]
  for j in [:r] do
    let v := j.toFloat / (r - 1).toFloat
    for i in [:c] do
      let u := i.toFloat / (c - 1).toFloat
      result := result.push (p.eval u v)
  return result

/-- Sample with normals for rendering. Returns array of (position, normal). -/
def sampleWithNormals (p : BezierPatch) (rows cols : Nat) : Array (Vec3 × Vec3) := Id.run do
  let r := if rows < 2 then 2 else rows
  let c := if cols < 2 then 2 else cols
  let mut result : Array (Vec3 × Vec3) := #[]
  for j in [:r] do
    let v := j.toFloat / (r - 1).toFloat
    for i in [:c] do
      let u := i.toFloat / (c - 1).toFloat
      let pos := p.eval u v
      let norm := p.normal u v
      result := result.push (pos, norm)
  return result

/-- Extract a cubic Bezier curve along constant u (u-isocurve). -/
def isocurveU (p : BezierPatch) (u : Float) : Bezier3 Vec3 :=
  let u' := Float.clamp u 0.0 1.0
  let evalRow (j : Nat) : Vec3 := Id.run do
    let mut result := Vec3.zero
    for i in [:4] do
      result := result.add ((p.getPoint j i).scale (bernstein3 i u'))
    return result
  { p0 := evalRow 0, p1 := evalRow 1, p2 := evalRow 2, p3 := evalRow 3 }

/-- Extract a cubic Bezier curve along constant v (v-isocurve). -/
def isocurveV (p : BezierPatch) (v : Float) : Bezier3 Vec3 :=
  let v' := Float.clamp v 0.0 1.0
  let evalCol (i : Nat) : Vec3 := Id.run do
    let mut result := Vec3.zero
    for j in [:4] do
      result := result.add ((p.getPoint j i).scale (bernstein3 j v'))
    return result
  { p0 := evalCol 0, p1 := evalCol 1, p2 := evalCol 2, p3 := evalCol 3 }

/-- Check if the patch is valid (has exactly 16 control points). -/
def isValid (p : BezierPatch) : Bool := p.controlPoints.size == 16

/-- Translate all control points. -/
def translate (p : BezierPatch) (offset : Vec3) : BezierPatch :=
  { controlPoints := p.controlPoints.map (· + offset) }

/-- Scale all control points from origin. -/
def scale (p : BezierPatch) (factor : Float) : BezierPatch :=
  { controlPoints := p.controlPoints.map (·.scale factor) }

end BezierPatch

end Linalg
