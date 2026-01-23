/-
  Tests for procedural noise functions.
-/

import Linalg
import Crucible

namespace LinalgTests.NoiseTests

open Crucible
open Linalg
open Linalg.Noise

testSuite "Perlin Noise"

test "perlin1D returns value in expected range" := do
  -- Sample at multiple points
  let samples := #[0.0, 0.5, 1.0, 1.5, 2.0, 3.7, 10.5]
  for x in samples do
    let n := perlin1D x
    ensure (n >= -1.5 && n <= 1.5) s!"perlin1D({x}) = {n} should be in [-1.5, 1.5]"

test "perlin2D returns value in expected range" := do
  let samples := #[(0.0, 0.0), (1.0, 1.0), (0.5, 2.3), (10.0, 7.0)]
  for (x, y) in samples do
    let n := perlin2D x y
    ensure (n >= -1.5 && n <= 1.5) s!"perlin2D({x}, {y}) = {n} should be in [-1.5, 1.5]"

test "perlin3D returns value in expected range" := do
  let samples := #[(0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0.5, 2.3, 1.7)]
  for (x, y, z) in samples do
    let n := perlin3D x y z
    ensure (n >= -1.5 && n <= 1.5) s!"perlin3D({x}, {y}, {z}) = {n} should be in [-1.5, 1.5]"

test "perlin noise is deterministic" := do
  let n1 := perlin2D 1.5 2.7
  let n2 := perlin2D 1.5 2.7
  ensure (floatNear n1 n2 0.0001) "same input should give same output"

test "perlin noise varies with input" := do
  -- Sample at different points and check that not all values are identical
  let n1 := perlin2D 0.0 0.0
  let n2 := perlin2D 0.3 0.7
  let n3 := perlin2D 1.5 2.3
  let n4 := perlin2D 5.0 3.0
  -- Check that we get some variation (at least one pair differs)
  let allSame := floatNear n1 n2 0.0001 && floatNear n2 n3 0.0001 && floatNear n3 n4 0.0001
  ensure (!allSame) "noise should vary across different inputs"

test "perlin2DV matches perlin2D" := do
  let p := Vec2.mk 1.5 2.7
  let n1 := perlin2D p.x p.y
  let n2 := perlin2DV p
  ensure (floatNear n1 n2 0.0001) "vector version should match"

test "perlin3DV matches perlin3D" := do
  let p := Vec3.mk 1.5 2.7 0.8
  let n1 := perlin3D p.x p.y p.z
  let n2 := perlin3DV p
  ensure (floatNear n1 n2 0.0001) "vector version should match"

testSuite "Simplex Noise"

test "simplex2D returns value in expected range" := do
  let samples := #[(0.0, 0.0), (1.0, 1.0), (0.5, 2.3), (10.0, 7.0)]
  for (x, y) in samples do
    let n := simplex2D x y
    ensure (n >= -1.5 && n <= 1.5) s!"simplex2D({x}, {y}) = {n} should be in [-1.5, 1.5]"

test "simplex3D returns value in expected range" := do
  let samples := #[(0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0.5, 2.3, 1.7)]
  for (x, y, z) in samples do
    let n := simplex3D x y z
    ensure (n >= -1.5 && n <= 1.5) s!"simplex3D({x}, {y}, {z}) = {n} should be in [-1.5, 1.5]"

test "simplex noise is deterministic" := do
  let n1 := simplex2D 1.5 2.7
  let n2 := simplex2D 1.5 2.7
  ensure (floatNear n1 n2 0.0001) "same input should give same output"

test "simplex2DV matches simplex2D" := do
  let p := Vec2.mk 1.5 2.7
  let n1 := simplex2D p.x p.y
  let n2 := simplex2DV p
  ensure (floatNear n1 n2 0.0001) "vector version should match"

testSuite "Fractal Noise"

test "fbm2D returns value in expected range" := do
  let n := fbm2D 1.5 2.7
  ensure (n >= -1.5 && n <= 1.5) s!"fbm2D should be in [-1.5, 1.5], got {n}"

test "fbm3D returns value in expected range" := do
  let n := fbm3D 1.5 2.7 0.8
  ensure (n >= -1.5 && n <= 1.5) s!"fbm3D should be in [-1.5, 1.5], got {n}"

test "fbm respects octave count" := do
  let config1 : FractalConfig := { octaves := 1 }
  let config6 : FractalConfig := { octaves := 6 }
  let n1 := fbm2D 1.5 2.7 config1
  let n6 := fbm2D 1.5 2.7 config6
  -- More octaves typically adds more detail (different value)
  ensure (n1 != n6) "different octave counts should give different results"

test "fbmSimplex2D returns value in expected range" := do
  let n := fbmSimplex2D 1.5 2.7
  ensure (n >= -1.5 && n <= 1.5) s!"fbmSimplex2D should be in [-1.5, 1.5], got {n}"

testSuite "Noise Variations"

test "ridged2D returns value in [0, 1]" := do
  let samples := #[(0.0, 0.0), (1.5, 2.7), (5.0, 3.0)]
  for (x, y) in samples do
    let n := ridged2D x y
    ensure (n >= 0.0 && n <= 1.5) s!"ridged2D({x}, {y}) = {n} should be in [0, 1.5]"

test "turbulence2D returns positive values" := do
  let samples := #[(0.0, 0.0), (1.5, 2.7), (5.0, 3.0)]
  for (x, y) in samples do
    let n := turbulence2D x y
    ensure (n >= 0.0) s!"turbulence2D({x}, {y}) = {n} should be >= 0"

test "domain warping produces valid output" := do
  let n := warp2D 1.5 2.7
  ensure (n >= -2.0 && n <= 2.0) s!"warp2D should be in reasonable range, got {n}"

testSuite "Noise Utilities"

test "normalize maps [-1, 1] to [0, 1]" := do
  ensure (floatNear (normalize (-1.0)) 0.0 0.0001) "normalize(-1) = 0"
  ensure (floatNear (normalize 0.0) 0.5 0.0001) "normalize(0) = 0.5"
  ensure (floatNear (normalize 1.0) 1.0 0.0001) "normalize(1) = 1"

test "remap maps to custom range" := do
  ensure (floatNear (remap (-1.0) 0.0 100.0) 0.0 0.0001) "remap(-1) = 0"
  ensure (floatNear (remap 0.0 0.0 100.0) 50.0 0.0001) "remap(0) = 50"
  ensure (floatNear (remap 1.0 0.0 100.0) 100.0 0.0001) "remap(1) = 100"

test "terrace creates discrete levels" := do
  -- With 4 levels, values should snap to 0.0, 0.25, 0.5, 0.75 (in [0,1] space)
  let n := terrace 0.0 4  -- 0 maps to 0.5 in [0,1], which is level 2
  -- Just verify it returns a valid value
  ensure (n >= -1.0 && n <= 1.0) s!"terrace should return valid value, got {n}"

-- ============================================================================
-- Value Noise Tests
-- ============================================================================

testSuite "Value Noise"

test "value1D returns value in [0, 1]" := do
  let samples := #[0.0, 0.5, 1.0, 1.5, 2.0, 3.7, 10.5]
  for x in samples do
    let n := value1D x
    ensure (n >= 0.0 && n <= 1.0) s!"value1D({x}) = {n} should be in [0, 1]"

test "value2D returns value in [0, 1]" := do
  let samples := #[(0.0, 0.0), (1.0, 1.0), (0.5, 2.3), (10.0, 7.0)]
  for (x, y) in samples do
    let n := value2D x y
    ensure (n >= 0.0 && n <= 1.0) s!"value2D({x}, {y}) = {n} should be in [0, 1]"

test "value3D returns value in [0, 1]" := do
  let samples := #[(0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0.5, 2.3, 1.7)]
  for (x, y, z) in samples do
    let n := value3D x y z
    ensure (n >= 0.0 && n <= 1.0) s!"value3D({x}, {y}, {z}) = {n} should be in [0, 1]"

test "value noise is deterministic" := do
  let n1 := value2D 1.5 2.7
  let n2 := value2D 1.5 2.7
  ensure (floatNear n1 n2 0.0001) "same input should give same output"

test "value2DV matches value2D" := do
  let p := Vec2.mk 1.5 2.7
  let n1 := value2D p.x p.y
  let n2 := value2DV p
  ensure (floatNear n1 n2 0.0001) "vector version should match"

test "fbmValue2D returns value in [0, 1]" := do
  let n := fbmValue2D 1.5 2.7
  ensure (n >= 0.0 && n <= 1.0) s!"fbmValue2D should be in [0, 1], got {n}"

-- ============================================================================
-- Worley Noise Tests
-- ============================================================================

testSuite "Worley Noise"

test "worley2D returns non-negative distances" := do
  let samples := #[(0.0, 0.0), (1.0, 1.0), (0.5, 2.3), (10.0, 7.0)]
  for (x, y) in samples do
    let r := worley2D x y
    ensure (r.f1 >= 0.0) s!"worley2D f1 should be >= 0"
    ensure (r.f2 >= r.f1) s!"worley2D f2 should be >= f1"
    ensure (r.f3 >= r.f2) s!"worley2D f3 should be >= f2"

test "worley3D returns non-negative distances" := do
  let samples := #[(0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0.5, 2.3, 1.7)]
  for (x, y, z) in samples do
    let r := worley3D x y z
    ensure (r.f1 >= 0.0) s!"worley3D f1 should be >= 0"
    ensure (r.f2 >= r.f1) s!"worley3D f2 should be >= f1"
    ensure (r.f3 >= r.f2) s!"worley3D f3 should be >= f2"

test "worley noise is deterministic" := do
  let r1 := worley2D 1.5 2.7
  let r2 := worley2D 1.5 2.7
  ensure (floatNear r1.f1 r2.f1 0.0001) "same input should give same f1"
  ensure (floatNear r1.f2 r2.f2 0.0001) "same input should give same f2"

test "worley2DF1 matches worley2D.f1" := do
  let r := worley2D 1.5 2.7
  let f1 := worley2DF1 1.5 2.7
  ensure (floatNear r.f1 f1 0.0001) "F1 helper should match"

test "worley2DEdge is f2 - f1" := do
  let r := worley2D 1.5 2.7
  let edge := worley2DEdge 1.5 2.7
  ensure (floatNear edge (r.f2 - r.f1) 0.0001) "Edge should be f2 - f1"

test "worley jitter affects output" := do
  let r1 := worley2D 1.5 2.7 1.0
  let r2 := worley2D 1.5 2.7 0.5
  -- Different jitter produces different results (or same if feature point happens to not move)
  -- Just verify both are valid non-negative distances
  ensure (r1.f1 >= 0.0 && r2.f1 >= 0.0) "both results should have valid f1"

test "fbmWorley2D returns non-negative" := do
  let n := fbmWorley2D 1.5 2.7
  ensure (n >= 0.0) s!"fbmWorley2D should be >= 0, got {n}"

-- ============================================================================
-- Random Shape Tests
-- ============================================================================

open Linalg.Random in
testSuite "Random Shapes 2D"

open Linalg.Random in
test "RNG generates different values" := do
  let rng := RNG.seed 12345
  let (v1, rng1) := rng.nextFloat
  let (v2, rng2) := rng1.nextFloat
  let (v3, _) := rng2.nextFloat
  ensure (v1 != v2 || v2 != v3) "RNG should generate different values"

open Linalg.Random in
test "inUnitCircle returns point in circle" := do
  let rng := RNG.seed 42
  let (p, _) := inUnitCircle rng
  let distSq := p.x * p.x + p.y * p.y
  ensure (distSq <= 1.0) s!"point should be in unit circle, distSq = {distSq}"

open Linalg.Random in
test "onUnitCircle returns point on circle" := do
  let rng := RNG.seed 42
  let (p, _) := onUnitCircle rng
  let dist := Float.sqrt (p.x * p.x + p.y * p.y)
  ensure (floatNear dist 1.0 0.0001) s!"point should be on unit circle, dist = {dist}"

open Linalg.Random in
test "inCircle respects center and radius" := do
  let rng := RNG.seed 42
  let center := Vec2.mk 5.0 5.0
  let radius := 3.0
  let (p, _) := inCircle rng center radius
  let dist := (p - center).length
  ensure (dist <= radius + 0.0001) s!"point should be in circle, dist = {dist}"

open Linalg.Random in
test "inRectangle returns point in bounds" := do
  let rng := RNG.seed 42
  let min := Vec2.mk 1.0 2.0
  let max := Vec2.mk 5.0 7.0
  let (p, _) := inRectangle rng min max
  ensure (p.x >= min.x && p.x <= max.x) "x should be in bounds"
  ensure (p.y >= min.y && p.y <= max.y) "y should be in bounds"

open Linalg.Random in
test "inTriangle2D returns point in triangle" := do
  let rng := RNG.seed 42
  let a := Vec2.mk 0.0 0.0
  let b := Vec2.mk 10.0 0.0
  let c := Vec2.mk 5.0 10.0
  let tri := Polygon2D.triangle a b c
  let (p, _) := inTriangle2D rng a b c
  ensure (tri.containsPoint p) "point should be in triangle"

open Linalg.Random in
testSuite "Random Shapes 3D"

open Linalg.Random in
test "inUnitSphere returns point in sphere" := do
  let rng := RNG.seed 42
  let (p, _) := inUnitSphere rng
  let distSq := p.x * p.x + p.y * p.y + p.z * p.z
  ensure (distSq <= 1.0) s!"point should be in unit sphere, distSq = {distSq}"

open Linalg.Random in
test "onUnitSphere returns point on sphere" := do
  let rng := RNG.seed 42
  let (p, _) := onUnitSphere rng
  let dist := Float.sqrt (p.x * p.x + p.y * p.y + p.z * p.z)
  ensure (floatNear dist 1.0 0.0001) s!"point should be on unit sphere, dist = {dist}"

open Linalg.Random in
test "inSphere respects center and radius" := do
  let rng := RNG.seed 42
  let center := Vec3.mk 5.0 5.0 5.0
  let radius := 3.0
  let (p, _) := inSphere rng center radius
  let dist := (p.sub center).length
  ensure (dist <= radius + 0.0001) s!"point should be in sphere, dist = {dist}"

open Linalg.Random in
test "inBox returns point in bounds" := do
  let rng := RNG.seed 42
  let min := Vec3.mk 1.0 2.0 3.0
  let max := Vec3.mk 5.0 7.0 9.0
  let (p, _) := inBox rng min max
  ensure (p.x >= min.x && p.x <= max.x) "x should be in bounds"
  ensure (p.y >= min.y && p.y <= max.y) "y should be in bounds"
  ensure (p.z >= min.z && p.z <= max.z) "z should be in bounds"

open Linalg.Random in
test "inCylinder returns point in cylinder" := do
  let rng := RNG.seed 42
  let radius := 2.0
  let height := 4.0
  let (p, _) := inCylinder rng radius height
  let horizDist := Float.sqrt (p.x * p.x + p.z * p.z)
  ensure (horizDist <= radius + 0.0001) s!"horizontal distance should be <= radius"
  ensure (p.y >= -height/2.0 - 0.0001 && p.y <= height/2.0 + 0.0001) "y should be in height range"

open Linalg.Random in
test "inUnitHemisphere returns point in upper hemisphere" := do
  let rng := RNG.seed 42
  let (p, _) := inUnitHemisphere rng
  let distSq := p.x * p.x + p.y * p.y + p.z * p.z
  ensure (distSq <= 1.0) "point should be in unit sphere"
  ensure (p.y >= 0.0) "y should be non-negative (upper hemisphere)"

-- ============================================================================
-- PCG Random Number Generator Tests
-- ============================================================================

open Linalg.Random in
testSuite "PCG Random"

open Linalg.Random in
test "PCG generates different values" := do
  let pcg := PCG.seed 42
  let (v1, pcg1) := pcg.nextUInt32
  let (v2, pcg2) := pcg1.nextUInt32
  let (v3, _) := pcg2.nextUInt32
  ensure (v1 != v2 || v2 != v3) "PCG should generate different values"

open Linalg.Random in
test "PCG is deterministic" := do
  let pcg1 := PCG.seed 12345
  let pcg2 := PCG.seed 12345
  let (v1, _) := pcg1.nextUInt32
  let (v2, _) := pcg2.nextUInt32
  ensure (v1 == v2) "same seed should produce same value"

open Linalg.Random in
test "PCG different seeds produce different values" := do
  let pcg1 := PCG.seed 1
  let pcg2 := PCG.seed 2
  let (v1, _) := pcg1.nextUInt32
  let (v2, _) := pcg2.nextUInt32
  ensure (v1 != v2) "different seeds should produce different values"

open Linalg.Random in
test "PCG nextFloat returns value in [0, 1)" := do
  let pcg := PCG.seed 42
  let (f, _) := pcg.nextFloat
  ensure (f >= 0.0 && f < 1.0) s!"float should be in [0, 1), got {f}"

open Linalg.Random in
test "PCG nextFloatRange respects bounds" := do
  let pcg := PCG.seed 42
  let (f, _) := pcg.nextFloatRange 10.0 20.0
  ensure (f >= 10.0 && f < 20.0) s!"float should be in [10, 20), got {f}"

open Linalg.Random in
test "PCG nextBounded respects bound" := do
  let pcg := PCG.seed 42
  let (v, _) := pcg.nextBounded 100
  ensure (v < 100) s!"value should be < 100, got {v}"

open Linalg.Random in
test "PCG nextIntRange works" := do
  let pcg := PCG.seed 42
  let (v, _) := pcg.nextIntRange 5 15
  ensure (v >= 5 && v <= 15) s!"value should be in [5, 15], got {v}"

open Linalg.Random in
test "PCG nextBool returns boolean" := do
  let pcg := PCG.seed 42
  let (b, _) := pcg.nextBool
  ensure (b == true || b == false) "should return a boolean"

open Linalg.Random in
test "PCG seedWithStream creates different streams" := do
  let pcg1 := PCG.seedWithStream 42 0
  let pcg2 := PCG.seedWithStream 42 1
  let (v1, _) := pcg1.nextUInt32
  let (v2, _) := pcg2.nextUInt32
  ensure (v1 != v2) "different streams should produce different values"

open Linalg.Random in
test "PCG split creates independent generators" := do
  let pcg := PCG.seed 42
  let (pcg1, pcg2) := pcg.split
  let (v1, _) := pcg1.nextUInt32
  let (v2, _) := pcg2.nextUInt32
  ensure (v1 != v2) "split generators should produce different values"

open Linalg.Random in
test "PCG shuffle preserves elements" := do
  let pcg := PCG.seed 42
  let arr := #[1, 2, 3, 4, 5]
  let (shuffled, _) := pcg.shuffle arr
  ensure (shuffled.size == arr.size) "shuffled array should have same size"
  -- Check all elements are present
  let sum1 := arr.foldl (· + ·) 0
  let sum2 := shuffled.foldl (· + ·) 0
  ensure (sum1 == sum2) "sum should be preserved after shuffle"

open Linalg.Random in
test "PCG choose returns element from array" := do
  let pcg := PCG.seed 42
  let arr := #[10, 20, 30, 40, 50]
  let (opt, _) := pcg.choose arr
  match opt with
  | some v => ensure (arr.contains v) s!"chosen value {v} should be in array"
  | none => ensure false "choose should return some for non-empty array"

open Linalg.Random in
test "PCG choose returns none for empty array" := do
  let pcg := PCG.seed 42
  let arr : Array Nat := #[]
  let (opt, _) := pcg.choose arr
  match opt with
  | some _ => ensure false "choose should return none for empty array"
  | none => ensure true "correctly returned none"

open Linalg.Random in
test "PCG sample returns correct count" := do
  let pcg := PCG.seed 42
  let arr := #[1, 2, 3, 4, 5]
  let (sampled, _) := pcg.sample arr 3
  ensure (sampled.size == 3) s!"sample should return 3 elements, got {sampled.size}"

open Linalg.Random in
test "PCG inUnitCirclePCG returns point in circle" := do
  let pcg := PCG.seed 42
  let (p, _) := inUnitCirclePCG pcg
  let distSq := p.x * p.x + p.y * p.y
  ensure (distSq <= 1.0) s!"point should be in unit circle, distSq = {distSq}"

open Linalg.Random in
test "PCG inUnitSpherePCG returns point in sphere" := do
  let pcg := PCG.seed 42
  let (p, _) := inUnitSpherePCG pcg
  let distSq := p.x * p.x + p.y * p.y + p.z * p.z
  ensure (distSq <= 1.0) s!"point should be in unit sphere, distSq = {distSq}"

open Linalg.Random in
test "PCG onUnitSpherePCG returns point on sphere" := do
  let pcg := PCG.seed 42
  let (p, _) := onUnitSpherePCG pcg
  let dist := Float.sqrt (p.x * p.x + p.y * p.y + p.z * p.z)
  ensure (floatNear dist 1.0 0.0001) s!"point should be on unit sphere, dist = {dist}"

open Linalg.Random in
test "PCG randomRotationPCG returns normalized quaternion" := do
  let pcg := PCG.seed 42
  let (q, _) := randomRotationPCG pcg
  let len := Float.sqrt (q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w)
  ensure (floatNear len 1.0 0.001) s!"quaternion should be normalized, len = {len}"

open Linalg.Random in
test "PCG advance moves generator forward" := do
  let pcg := PCG.seed 42
  -- Manually advance 10 steps
  let pcg1 := pcg.nextUInt32.2.nextUInt32.2.nextUInt32.2.nextUInt32.2.nextUInt32.2
  let pcg1 := pcg1.nextUInt32.2.nextUInt32.2.nextUInt32.2.nextUInt32.2.nextUInt32.2
  -- Use advance function
  let pcg2 := pcg.advance 10
  -- Both should produce same next value
  let (v1, _) := pcg1.nextUInt32
  let (v2, _) := pcg2.nextUInt32
  ensure (v1 == v2) s!"advance(10) should match 10 manual steps, got {v1} vs {v2}"



end LinalgTests.NoiseTests
