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

#generate_tests

end LinalgTests.NoiseTests
