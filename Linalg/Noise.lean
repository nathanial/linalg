/-
  Procedural noise functions for terrain, textures, and effects.

  Includes:
  - Perlin noise (classic gradient noise) in 1D, 2D, 3D
  - Simplex noise (improved gradient noise) in 2D, 3D
  - Fractal/FBM noise (layered octaves)
  - Utility functions (ridged, turbulence, domain warping)
-/

import Linalg.Core
import Linalg.Vec2
import Linalg.Vec3

namespace Linalg

namespace Noise

-- ============================================================================
-- Permutation Table
-- ============================================================================

/-- Standard permutation table (Ken Perlin's original). -/
private def permutation : Array UInt8 := #[
  151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225,
  140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148,
  247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32,
  57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175,
  74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122,
  60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54,
  65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169,
  200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64,
  52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212,
  207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213,
  119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
  129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104,
  218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241,
  81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157,
  184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93,
  222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
]

/-- Get permutation value with wrapping. -/
@[inline]
private def perm (i : Int) : UInt8 :=
  let idx := (i % 256 + 256) % 256
  permutation.getD idx.toNat 0

/-- Get permutation as Float in [0, 1]. -/
@[inline]
private def permFloat (i : Int) : Float :=
  (perm i).toFloat / 255.0

-- ============================================================================
-- Gradient Functions
-- ============================================================================

/-- 1D gradient. -/
@[inline]
private def grad1 (hash : UInt8) (x : Float) : Float :=
  if hash &&& 1 == 0 then x else -x

/-- 2D gradient vectors. -/
@[inline]
private def grad2 (hash : UInt8) (x y : Float) : Float :=
  let h := hash &&& 7
  let u := if h < 4 then x else y
  let v := if h < 4 then y else x
  (if h &&& 1 == 0 then u else -u) + (if h &&& 2 == 0 then v else -v)

/-- 3D gradient vectors. -/
@[inline]
private def grad3 (hash : UInt8) (x y z : Float) : Float :=
  let h := hash &&& 15
  let u := if h < 8 then x else y
  let v := if h < 4 then y else if h == 12 || h == 14 then x else z
  (if h &&& 1 == 0 then u else -u) + (if h &&& 2 == 0 then v else -v)

-- ============================================================================
-- Interpolation
-- ============================================================================

/-- Fade function (6t^5 - 15t^4 + 10t^3) for smooth interpolation. -/
@[inline]
private def fade (t : Float) : Float :=
  t * t * t * (t * (t * 6.0 - 15.0) + 10.0)

-- ============================================================================
-- Perlin Noise
-- ============================================================================

/-- 1D Perlin noise. Returns value in approximately [-1, 1]. -/
def perlin1D (x : Float) : Float :=
  let xi := Float.floor x |> Float.toUInt64 |> UInt64.toNat |> Int.ofNat
  let xf := x - Float.floor x

  let u := fade xf

  let a := perm xi
  let b := perm (xi + 1)

  let ga := grad1 a xf
  let gb := grad1 b (xf - 1.0)

  Float.lerp ga gb u

/-- 2D Perlin noise. Returns value in approximately [-1, 1]. -/
def perlin2D (x y : Float) : Float :=
  let xi := Float.floor x |> Float.toUInt64 |> UInt64.toNat |> Int.ofNat
  let yi := Float.floor y |> Float.toUInt64 |> UInt64.toNat |> Int.ofNat
  let xf := x - Float.floor x
  let yf := y - Float.floor y

  let u := fade xf
  let v := fade yf

  let aa := perm (xi + (perm yi).toNat)
  let ab := perm (xi + (perm (yi + 1)).toNat)
  let ba := perm (xi + 1 + (perm yi).toNat)
  let bb := perm (xi + 1 + (perm (yi + 1)).toNat)

  let x1 := Float.lerp (grad2 aa xf yf) (grad2 ba (xf - 1.0) yf) u
  let x2 := Float.lerp (grad2 ab xf (yf - 1.0)) (grad2 bb (xf - 1.0) (yf - 1.0)) u

  Float.lerp x1 x2 v

/-- 2D Perlin noise from Vec2. -/
def perlin2DV (p : Vec2) : Float := perlin2D p.x p.y

/-- 3D Perlin noise. Returns value in approximately [-1, 1]. -/
def perlin3D (x y z : Float) : Float :=
  let xi := Float.floor x |> Float.toUInt64 |> UInt64.toNat |> Int.ofNat
  let yi := Float.floor y |> Float.toUInt64 |> UInt64.toNat |> Int.ofNat
  let zi := Float.floor z |> Float.toUInt64 |> UInt64.toNat |> Int.ofNat
  let xf := x - Float.floor x
  let yf := y - Float.floor y
  let zf := z - Float.floor z

  let u := fade xf
  let v := fade yf
  let w := fade zf

  let aaa := perm (xi + (perm (yi + (perm zi).toNat)).toNat)
  let aab := perm (xi + (perm (yi + (perm (zi + 1)).toNat)).toNat)
  let aba := perm (xi + (perm (yi + 1 + (perm zi).toNat)).toNat)
  let abb := perm (xi + (perm (yi + 1 + (perm (zi + 1)).toNat)).toNat)
  let baa := perm (xi + 1 + (perm (yi + (perm zi).toNat)).toNat)
  let bab := perm (xi + 1 + (perm (yi + (perm (zi + 1)).toNat)).toNat)
  let bba := perm (xi + 1 + (perm (yi + 1 + (perm zi).toNat)).toNat)
  let bbb := perm (xi + 1 + (perm (yi + 1 + (perm (zi + 1)).toNat)).toNat)

  let x1 := Float.lerp (grad3 aaa xf yf zf) (grad3 baa (xf - 1.0) yf zf) u
  let x2 := Float.lerp (grad3 aba xf (yf - 1.0) zf) (grad3 bba (xf - 1.0) (yf - 1.0) zf) u
  let y1 := Float.lerp x1 x2 v

  let x3 := Float.lerp (grad3 aab xf yf (zf - 1.0)) (grad3 bab (xf - 1.0) yf (zf - 1.0)) u
  let x4 := Float.lerp (grad3 abb xf (yf - 1.0) (zf - 1.0)) (grad3 bbb (xf - 1.0) (yf - 1.0) (zf - 1.0)) u
  let y2 := Float.lerp x3 x4 v

  Float.lerp y1 y2 w

/-- 3D Perlin noise from Vec3. -/
def perlin3DV (p : Vec3) : Float := perlin3D p.x p.y p.z

-- ============================================================================
-- Simplex Noise
-- ============================================================================

/-- Skewing factor for 2D simplex noise. -/
private def F2 : Float := 0.5 * (Float.sqrt 3.0 - 1.0)

/-- Unskewing factor for 2D simplex noise. -/
private def G2 : Float := (3.0 - Float.sqrt 3.0) / 6.0

/-- 2D Simplex noise. Returns value in approximately [-1, 1]. -/
def simplex2D (x y : Float) : Float := Id.run do
  -- Skew input space to determine simplex cell
  let s := (x + y) * F2
  let i := Float.floor (x + s)
  let j := Float.floor (y + s)

  -- Unskew back to (x, y) space
  let t := (i + j) * G2
  let x0 := x - (i - t)
  let y0 := y - (j - t)

  -- Determine which simplex we're in
  let (i1, j1) := if x0 > y0 then (1.0, 0.0) else (0.0, 1.0)

  -- Offsets for corners
  let x1 := x0 - i1 + G2
  let y1 := y0 - j1 + G2
  let x2 := x0 - 1.0 + 2.0 * G2
  let y2 := y0 - 1.0 + 2.0 * G2

  -- Hash coordinates
  let ii := (Float.toUInt64 i).toNat % 256
  let jj := (Float.toUInt64 j).toNat % 256

  let gi0 := perm (ii + (perm jj).toNat)
  let gi1 := perm (ii + i1.toUInt64.toNat + (perm (jj + j1.toUInt64.toNat)).toNat)
  let gi2 := perm (ii + 1 + (perm (jj + 1)).toNat)

  -- Calculate contribution from each corner
  let mut n := 0.0

  let t0 := 0.5 - x0 * x0 - y0 * y0
  if t0 >= 0.0 then
    let t0' := t0 * t0
    n := n + t0' * t0' * grad2 gi0 x0 y0

  let t1 := 0.5 - x1 * x1 - y1 * y1
  if t1 >= 0.0 then
    let t1' := t1 * t1
    n := n + t1' * t1' * grad2 gi1 x1 y1

  let t2 := 0.5 - x2 * x2 - y2 * y2
  if t2 >= 0.0 then
    let t2' := t2 * t2
    n := n + t2' * t2' * grad2 gi2 x2 y2

  -- Scale to [-1, 1]
  return 70.0 * n

/-- 2D Simplex noise from Vec2. -/
def simplex2DV (p : Vec2) : Float := simplex2D p.x p.y

/-- Skewing factor for 3D simplex noise. -/
private def F3 : Float := 1.0 / 3.0

/-- Unskewing factor for 3D simplex noise. -/
private def G3 : Float := 1.0 / 6.0

/-- 3D Simplex noise. Returns value in approximately [-1, 1]. -/
def simplex3D (x y z : Float) : Float := Id.run do
  -- Skew input space
  let s := (x + y + z) * F3
  let i := Float.floor (x + s)
  let j := Float.floor (y + s)
  let k := Float.floor (z + s)

  -- Unskew back
  let t := (i + j + k) * G3
  let x0 := x - (i - t)
  let y0 := y - (j - t)
  let z0 := z - (k - t)

  -- Determine simplex
  let (i1, j1, k1, i2, j2, k2) :=
    if x0 >= y0 then
      if y0 >= z0 then (1.0, 0.0, 0.0, 1.0, 1.0, 0.0)
      else if x0 >= z0 then (1.0, 0.0, 0.0, 1.0, 0.0, 1.0)
      else (0.0, 0.0, 1.0, 1.0, 0.0, 1.0)
    else
      if y0 < z0 then (0.0, 0.0, 1.0, 0.0, 1.0, 1.0)
      else if x0 < z0 then (0.0, 1.0, 0.0, 0.0, 1.0, 1.0)
      else (0.0, 1.0, 0.0, 1.0, 1.0, 0.0)

  -- Offsets for corners
  let x1 := x0 - i1 + G3
  let y1 := y0 - j1 + G3
  let z1 := z0 - k1 + G3
  let x2 := x0 - i2 + 2.0 * G3
  let y2 := y0 - j2 + 2.0 * G3
  let z2 := z0 - k2 + 2.0 * G3
  let x3 := x0 - 1.0 + 3.0 * G3
  let y3 := y0 - 1.0 + 3.0 * G3
  let z3 := z0 - 1.0 + 3.0 * G3

  -- Hash coordinates
  let ii := (Float.toUInt64 i).toNat % 256
  let jj := (Float.toUInt64 j).toNat % 256
  let kk := (Float.toUInt64 k).toNat % 256

  let gi0 := perm (ii + (perm (jj + (perm kk).toNat)).toNat)
  let gi1 := perm (ii + i1.toUInt64.toNat + (perm (jj + j1.toUInt64.toNat + (perm (kk + k1.toUInt64.toNat)).toNat)).toNat)
  let gi2 := perm (ii + i2.toUInt64.toNat + (perm (jj + j2.toUInt64.toNat + (perm (kk + k2.toUInt64.toNat)).toNat)).toNat)
  let gi3 := perm (ii + 1 + (perm (jj + 1 + (perm (kk + 1)).toNat)).toNat)

  -- Calculate contributions
  let mut n := 0.0

  let t0 := 0.6 - x0 * x0 - y0 * y0 - z0 * z0
  if t0 >= 0.0 then
    let t0' := t0 * t0
    n := n + t0' * t0' * grad3 gi0 x0 y0 z0

  let t1 := 0.6 - x1 * x1 - y1 * y1 - z1 * z1
  if t1 >= 0.0 then
    let t1' := t1 * t1
    n := n + t1' * t1' * grad3 gi1 x1 y1 z1

  let t2 := 0.6 - x2 * x2 - y2 * y2 - z2 * z2
  if t2 >= 0.0 then
    let t2' := t2 * t2
    n := n + t2' * t2' * grad3 gi2 x2 y2 z2

  let t3 := 0.6 - x3 * x3 - y3 * y3 - z3 * z3
  if t3 >= 0.0 then
    let t3' := t3 * t3
    n := n + t3' * t3' * grad3 gi3 x3 y3 z3

  -- Scale to [-1, 1]
  return 32.0 * n

/-- 3D Simplex noise from Vec3. -/
def simplex3DV (p : Vec3) : Float := simplex3D p.x p.y p.z

-- ============================================================================
-- Fractal Brownian Motion (FBM)
-- ============================================================================

/-- Configuration for fractal noise. -/
structure FractalConfig where
  octaves : Nat := 6
  lacunarity : Float := 2.0    -- Frequency multiplier per octave
  persistence : Float := 0.5   -- Amplitude multiplier per octave (also called gain)
  deriving Repr, Inhabited

/-- Default fractal configuration. -/
def FractalConfig.default : FractalConfig := {}

/-- 2D Fractal Brownian Motion using Perlin noise. -/
def fbm2D (x y : Float) (config : FractalConfig := {}) : Float := Id.run do
  let mut total := 0.0
  let mut frequency := 1.0
  let mut amplitude := 1.0
  let mut maxValue := 0.0

  for _ in [:config.octaves] do
    total := total + perlin2D (x * frequency) (y * frequency) * amplitude
    maxValue := maxValue + amplitude
    amplitude := amplitude * config.persistence
    frequency := frequency * config.lacunarity

  return total / maxValue

/-- 2D FBM from Vec2. -/
def fbm2DV (p : Vec2) (config : FractalConfig := {}) : Float :=
  fbm2D p.x p.y config

/-- 3D Fractal Brownian Motion using Perlin noise. -/
def fbm3D (x y z : Float) (config : FractalConfig := {}) : Float := Id.run do
  let mut total := 0.0
  let mut frequency := 1.0
  let mut amplitude := 1.0
  let mut maxValue := 0.0

  for _ in [:config.octaves] do
    total := total + perlin3D (x * frequency) (y * frequency) (z * frequency) * amplitude
    maxValue := maxValue + amplitude
    amplitude := amplitude * config.persistence
    frequency := frequency * config.lacunarity

  return total / maxValue

/-- 3D FBM from Vec3. -/
def fbm3DV (p : Vec3) (config : FractalConfig := {}) : Float :=
  fbm3D p.x p.y p.z config

/-- 2D FBM using Simplex noise (often smoother). -/
def fbmSimplex2D (x y : Float) (config : FractalConfig := {}) : Float := Id.run do
  let mut total := 0.0
  let mut frequency := 1.0
  let mut amplitude := 1.0
  let mut maxValue := 0.0

  for _ in [:config.octaves] do
    total := total + simplex2D (x * frequency) (y * frequency) * amplitude
    maxValue := maxValue + amplitude
    amplitude := amplitude * config.persistence
    frequency := frequency * config.lacunarity

  return total / maxValue

/-- 3D FBM using Simplex noise. -/
def fbmSimplex3D (x y z : Float) (config : FractalConfig := {}) : Float := Id.run do
  let mut total := 0.0
  let mut frequency := 1.0
  let mut amplitude := 1.0
  let mut maxValue := 0.0

  for _ in [:config.octaves] do
    total := total + simplex3D (x * frequency) (y * frequency) (z * frequency) * amplitude
    maxValue := maxValue + amplitude
    amplitude := amplitude * config.persistence
    frequency := frequency * config.lacunarity

  return total / maxValue

-- ============================================================================
-- Noise Variations
-- ============================================================================

/-- Ridged multifractal noise (2D). Good for mountains and ridges. -/
def ridged2D (x y : Float) (config : FractalConfig := {}) : Float := Id.run do
  let mut total := 0.0
  let mut frequency := 1.0
  let mut amplitude := 1.0
  let mut maxValue := 0.0

  for _ in [:config.octaves] do
    let n := perlin2D (x * frequency) (y * frequency)
    -- Transform: ridge = 1 - |noise|
    let ridge := 1.0 - Float.abs' n
    -- Square for sharper ridges
    let ridge := ridge * ridge
    total := total + ridge * amplitude
    maxValue := maxValue + amplitude
    amplitude := amplitude * config.persistence
    frequency := frequency * config.lacunarity

  return total / maxValue

/-- Ridged multifractal noise (3D). -/
def ridged3D (x y z : Float) (config : FractalConfig := {}) : Float := Id.run do
  let mut total := 0.0
  let mut frequency := 1.0
  let mut amplitude := 1.0
  let mut maxValue := 0.0

  for _ in [:config.octaves] do
    let n := perlin3D (x * frequency) (y * frequency) (z * frequency)
    let ridge := 1.0 - Float.abs' n
    let ridge := ridge * ridge
    total := total + ridge * amplitude
    maxValue := maxValue + amplitude
    amplitude := amplitude * config.persistence
    frequency := frequency * config.lacunarity

  return total / maxValue

/-- Turbulence noise (2D). Sum of absolute values - good for fire, clouds. -/
def turbulence2D (x y : Float) (config : FractalConfig := {}) : Float := Id.run do
  let mut total := 0.0
  let mut frequency := 1.0
  let mut amplitude := 1.0
  let mut maxValue := 0.0

  for _ in [:config.octaves] do
    total := total + Float.abs' (perlin2D (x * frequency) (y * frequency)) * amplitude
    maxValue := maxValue + amplitude
    amplitude := amplitude * config.persistence
    frequency := frequency * config.lacunarity

  return total / maxValue

/-- Turbulence noise (3D). -/
def turbulence3D (x y z : Float) (config : FractalConfig := {}) : Float := Id.run do
  let mut total := 0.0
  let mut frequency := 1.0
  let mut amplitude := 1.0
  let mut maxValue := 0.0

  for _ in [:config.octaves] do
    total := total + Float.abs' (perlin3D (x * frequency) (y * frequency) (z * frequency)) * amplitude
    maxValue := maxValue + amplitude
    amplitude := amplitude * config.persistence
    frequency := frequency * config.lacunarity

  return total / maxValue

-- ============================================================================
-- Domain Warping
-- ============================================================================

/-- Domain warping: distort coordinates using noise before sampling.
    Creates organic, flowing patterns. -/
def warp2D (x y : Float) (strength : Float := 4.0) : Float :=
  let qx := fbm2D x y
  let qy := fbm2D (x + 5.2) (y + 1.3)
  fbm2D (x + strength * qx) (y + strength * qy)

/-- Advanced domain warping with two layers of distortion. -/
def warp2DAdvanced (x y : Float) (strength1 strength2 : Float := 4.0) : Float :=
  let qx := fbm2D x y
  let qy := fbm2D (x + 5.2) (y + 1.3)
  let rx := fbm2D (x + strength1 * qx + 1.7) (y + strength1 * qy + 9.2)
  let ry := fbm2D (x + strength1 * qx + 8.3) (y + strength1 * qy + 2.8)
  fbm2D (x + strength2 * rx) (y + strength2 * ry)

-- ============================================================================
-- Utility Functions
-- ============================================================================

/-- Remap noise from [-1, 1] to [0, 1]. -/
@[inline]
def normalize (n : Float) : Float := (n + 1.0) * 0.5

/-- Remap noise from [-1, 1] to [lo, hi]. -/
@[inline]
def remap (n : Float) (lo hi : Float) : Float :=
  lo + (n + 1.0) * 0.5 * (hi - lo)

/-- Apply a power curve to noise (in [0, 1] range).
    power < 1: flattens valleys, power > 1: sharpens peaks. -/
def redistribute (n : Float) (power : Float) : Float :=
  let n01 := normalize n
  Float.pow n01 power * 2.0 - 1.0

/-- Terrace/posterize noise into discrete levels. -/
def terrace (n : Float) (levels : Nat) : Float :=
  if levels == 0 then n
  else
    let n01 := normalize n
    let step := 1.0 / levels.toFloat
    let level := Float.floor (n01 / step) * step
    level * 2.0 - 1.0

end Noise

end Linalg
