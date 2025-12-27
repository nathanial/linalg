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
import Linalg.Quat

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

-- ============================================================================
-- Value Noise
-- ============================================================================

/-- 1D Value noise. Returns value in [0, 1].
    Simpler than Perlin - just interpolates random values at integer points. -/
def value1D (x : Float) : Float :=
  let xi := Float.floor x
  let xf := x - xi
  let t := fade xf  -- Use smooth interpolation

  let i := xi.toUInt64.toNat
  let v0 := permFloat i
  let v1 := permFloat (i + 1)

  Float.lerp v0 v1 t

/-- 2D Value noise. Returns value in [0, 1]. -/
def value2D (x y : Float) : Float :=
  let xi := Float.floor x
  let yi := Float.floor y
  let xf := x - xi
  let yf := y - yi
  let u := fade xf
  let v := fade yf

  let i := xi.toUInt64.toNat
  let j := yi.toUInt64.toNat

  let v00 := permFloat (i + (perm j).toNat)
  let v10 := permFloat (i + 1 + (perm j).toNat)
  let v01 := permFloat (i + (perm (j + 1)).toNat)
  let v11 := permFloat (i + 1 + (perm (j + 1)).toNat)

  let x1 := Float.lerp v00 v10 u
  let x2 := Float.lerp v01 v11 u
  Float.lerp x1 x2 v

/-- 2D Value noise from Vec2. -/
def value2DV (p : Vec2) : Float := value2D p.x p.y

/-- 3D Value noise. Returns value in [0, 1]. -/
def value3D (x y z : Float) : Float :=
  let xi := Float.floor x
  let yi := Float.floor y
  let zi := Float.floor z
  let xf := x - xi
  let yf := y - yi
  let zf := z - zi
  let u := fade xf
  let v := fade yf
  let w := fade zf

  let i := xi.toUInt64.toNat
  let j := yi.toUInt64.toNat
  let k := zi.toUInt64.toNat

  let v000 := permFloat (i + (perm (j + (perm k).toNat)).toNat)
  let v100 := permFloat (i + 1 + (perm (j + (perm k).toNat)).toNat)
  let v010 := permFloat (i + (perm (j + 1 + (perm k).toNat)).toNat)
  let v110 := permFloat (i + 1 + (perm (j + 1 + (perm k).toNat)).toNat)
  let v001 := permFloat (i + (perm (j + (perm (k + 1)).toNat)).toNat)
  let v101 := permFloat (i + 1 + (perm (j + (perm (k + 1)).toNat)).toNat)
  let v011 := permFloat (i + (perm (j + 1 + (perm (k + 1)).toNat)).toNat)
  let v111 := permFloat (i + 1 + (perm (j + 1 + (perm (k + 1)).toNat)).toNat)

  let x1 := Float.lerp v000 v100 u
  let x2 := Float.lerp v010 v110 u
  let y1 := Float.lerp x1 x2 v

  let x3 := Float.lerp v001 v101 u
  let x4 := Float.lerp v011 v111 u
  let y2 := Float.lerp x3 x4 v

  Float.lerp y1 y2 w

/-- 3D Value noise from Vec3. -/
def value3DV (p : Vec3) : Float := value3D p.x p.y p.z

/-- 2D Value FBM (fractal brownian motion with value noise). -/
def fbmValue2D (x y : Float) (config : FractalConfig := {}) : Float := Id.run do
  let mut total := 0.0
  let mut frequency := 1.0
  let mut amplitude := 1.0
  let mut maxValue := 0.0

  for _ in [:config.octaves] do
    total := total + value2D (x * frequency) (y * frequency) * amplitude
    maxValue := maxValue + amplitude
    amplitude := amplitude * config.persistence
    frequency := frequency * config.lacunarity

  return total / maxValue

/-- 3D Value FBM. -/
def fbmValue3D (x y z : Float) (config : FractalConfig := {}) : Float := Id.run do
  let mut total := 0.0
  let mut frequency := 1.0
  let mut amplitude := 1.0
  let mut maxValue := 0.0

  for _ in [:config.octaves] do
    total := total + value3D (x * frequency) (y * frequency) (z * frequency) * amplitude
    maxValue := maxValue + amplitude
    amplitude := amplitude * config.persistence
    frequency := frequency * config.lacunarity

  return total / maxValue

-- ============================================================================
-- Worley (Cellular) Noise
-- ============================================================================

/-- Hash function to get pseudo-random offset for Worley feature points. -/
@[inline]
private def worleyHash (i j : Int) (n : Nat) : Float :=
  let h := perm (i + (perm j).toNat * 17 + n * 31)
  h.toFloat / 255.0

/-- Hash function for 3D Worley. -/
@[inline]
private def worleyHash3 (i j k : Int) (n : Nat) : Float :=
  let h := perm (i + (perm (j + (perm k).toNat * 17)).toNat + n * 31)
  h.toFloat / 255.0

/-- Result of Worley noise calculation. -/
structure WorleyResult where
  /-- Distance to closest feature point. -/
  f1 : Float
  /-- Distance to second closest feature point. -/
  f2 : Float
  /-- Distance to third closest feature point. -/
  f3 : Float
deriving Repr, Inhabited

/-- Convert Int to Float. -/
@[inline]
private def intToFloat (i : Int) : Float :=
  if i >= 0 then i.toNat.toFloat
  else -((-i).toNat.toFloat)

/-- 2D Worley (Cellular) noise.
    Returns distances to the three closest feature points.
    Common uses:
    - f1: cell pattern (distance to nearest)
    - f2 - f1: cell edges
    - f2: inverse cells -/
def worley2D (x y : Float) (jitter : Float := 1.0) : WorleyResult := Id.run do
  let xi := Float.floor x |> Float.toUInt64 |> UInt64.toNat |> Int.ofNat
  let yi := Float.floor y |> Float.toUInt64 |> UInt64.toNat |> Int.ofNat

  let mut f1 := Float.infinity
  let mut f2 := Float.infinity
  let mut f3 := Float.infinity

  -- Check 3x3 neighborhood of cells
  for di in [-1, 0, 1] do
    for dj in [-1, 0, 1] do
      let ci := xi + di
      let cj := yi + dj

      -- Get pseudo-random feature point position within cell
      let px := intToFloat ci + worleyHash ci cj 0 * jitter
      let py := intToFloat cj + worleyHash ci cj 1 * jitter

      -- Distance to feature point
      let dx := px - x
      let dy := py - y
      let dist := Float.sqrt (dx * dx + dy * dy)

      -- Update closest distances
      if dist < f1 then
        f3 := f2
        f2 := f1
        f1 := dist
      else if dist < f2 then
        f3 := f2
        f2 := dist
      else if dist < f3 then
        f3 := dist

  return ⟨f1, f2, f3⟩

/-- 2D Worley noise from Vec2. -/
def worley2DV (p : Vec2) (jitter : Float := 1.0) : WorleyResult :=
  worley2D p.x p.y jitter

/-- 2D Worley noise returning just the distance to closest point (F1).
    Returns value typically in [0, ~1.5]. -/
def worley2DF1 (x y : Float) (jitter : Float := 1.0) : Float :=
  (worley2D x y jitter).f1

/-- 2D Worley noise returning cell edge pattern (F2 - F1). -/
def worley2DEdge (x y : Float) (jitter : Float := 1.0) : Float :=
  let r := worley2D x y jitter
  r.f2 - r.f1

/-- 3D Worley (Cellular) noise.
    Returns distances to the three closest feature points. -/
def worley3D (x y z : Float) (jitter : Float := 1.0) : WorleyResult := Id.run do
  let xi := Float.floor x |> Float.toUInt64 |> UInt64.toNat |> Int.ofNat
  let yi := Float.floor y |> Float.toUInt64 |> UInt64.toNat |> Int.ofNat
  let zi := Float.floor z |> Float.toUInt64 |> UInt64.toNat |> Int.ofNat

  let mut f1 := Float.infinity
  let mut f2 := Float.infinity
  let mut f3 := Float.infinity

  -- Check 3x3x3 neighborhood of cells
  for di in [-1, 0, 1] do
    for dj in [-1, 0, 1] do
      for dk in [-1, 0, 1] do
        let ci := xi + di
        let cj := yi + dj
        let ck := zi + dk

        -- Get pseudo-random feature point position within cell
        let px := intToFloat ci + worleyHash3 ci cj ck 0 * jitter
        let py := intToFloat cj + worleyHash3 ci cj ck 1 * jitter
        let pz := intToFloat ck + worleyHash3 ci cj ck 2 * jitter

        -- Distance to feature point
        let dx := px - x
        let dy := py - y
        let dz := pz - z
        let dist := Float.sqrt (dx * dx + dy * dy + dz * dz)

        -- Update closest distances
        if dist < f1 then
          f3 := f2
          f2 := f1
          f1 := dist
        else if dist < f2 then
          f3 := f2
          f2 := dist
        else if dist < f3 then
          f3 := dist

  return ⟨f1, f2, f3⟩

/-- 3D Worley noise from Vec3. -/
def worley3DV (p : Vec3) (jitter : Float := 1.0) : WorleyResult :=
  worley3D p.x p.y p.z jitter

/-- 3D Worley noise returning just F1 (distance to closest). -/
def worley3DF1 (x y z : Float) (jitter : Float := 1.0) : Float :=
  (worley3D x y z jitter).f1

/-- 3D Worley noise returning cell edge pattern (F2 - F1). -/
def worley3DEdge (x y z : Float) (jitter : Float := 1.0) : Float :=
  let r := worley3D x y z jitter
  r.f2 - r.f1

/-- 2D Worley FBM (fractal cellular noise). -/
def fbmWorley2D (x y : Float) (config : FractalConfig := {}) (jitter : Float := 1.0) : Float := Id.run do
  let mut total := 0.0
  let mut frequency := 1.0
  let mut amplitude := 1.0
  let mut maxValue := 0.0

  for _ in [:config.octaves] do
    total := total + worley2DF1 (x * frequency) (y * frequency) jitter * amplitude
    maxValue := maxValue + amplitude
    amplitude := amplitude * config.persistence
    frequency := frequency * config.lacunarity

  return total / maxValue

/-- 3D Worley FBM. -/
def fbmWorley3D (x y z : Float) (config : FractalConfig := {}) (jitter : Float := 1.0) : Float := Id.run do
  let mut total := 0.0
  let mut frequency := 1.0
  let mut amplitude := 1.0
  let mut maxValue := 0.0

  for _ in [:config.octaves] do
    total := total + worley3DF1 (x * frequency) (y * frequency) (z * frequency) jitter * amplitude
    maxValue := maxValue + amplitude
    amplitude := amplitude * config.persistence
    frequency := frequency * config.lacunarity

  return total / maxValue

end Noise

-- ============================================================================
-- Random Point Generation in Shapes
-- ============================================================================

namespace Random

/-- Simple random number generator state using the permutation table.
    For game use cases where quality is less critical than speed. -/
structure RNG where
  state : UInt64
deriving Repr, Inhabited

namespace RNG

/-- Create RNG with seed. -/
def seed (s : UInt64) : RNG := ⟨s⟩

/-- Generate next random UInt64 (xorshift64). -/
def nextUInt64 (rng : RNG) : UInt64 × RNG :=
  let x := rng.state
  let x := x ^^^ (x >>> 12)
  let x := x ^^^ (x <<< 25)
  let x := x ^^^ (x >>> 27)
  (x * 0x2545F4914F6CDD1D, ⟨x⟩)

/-- Generate random Float in [0, 1). -/
def nextFloat (rng : RNG) : Float × RNG :=
  let (u, rng') := rng.nextUInt64
  ((u.toFloat / UInt64.toFloat (UInt64.ofNat 0xFFFFFFFFFFFFFFFF)), rng')

/-- Generate random Float in [lo, hi). -/
def nextFloatRange (rng : RNG) (lo hi : Float) : Float × RNG :=
  let (f, rng') := rng.nextFloat
  (lo + f * (hi - lo), rng')

/-- Generate two random Floats in [0, 1). -/
def nextFloat2 (rng : RNG) : (Float × Float) × RNG :=
  let (f1, rng1) := rng.nextFloat
  let (f2, rng2) := rng1.nextFloat
  ((f1, f2), rng2)

/-- Generate three random Floats in [0, 1). -/
def nextFloat3 (rng : RNG) : (Float × Float × Float) × RNG :=
  let (f1, rng1) := rng.nextFloat
  let (f2, rng2) := rng1.nextFloat
  let (f3, rng3) := rng2.nextFloat
  ((f1, f2, f3), rng3)

end RNG

-- ============================================================================
-- PCG Random Number Generator
-- ============================================================================

/-- PCG (Permuted Congruential Generator) - A high-quality, fast PRNG.
    Uses PCG-XSH-RR variant: 64-bit state, 32-bit output.
    Better statistical quality than xorshift while remaining fast.

    Reference: https://www.pcg-random.org/ -/
structure PCG where
  /-- Internal 64-bit state -/
  state : UInt64
  /-- Stream selector (increment, must be odd) -/
  inc : UInt64
deriving Repr, Inhabited

namespace PCG

/-- The LCG multiplier used by PCG. -/
private def multiplier : UInt64 := 6364136223846793005

/-- Create a PCG generator from a seed.
    Uses a default stream. -/
def seed (s : UInt64) : PCG :=
  -- Initialize with default increment (must be odd)
  let inc : UInt64 := 1442695040888963407
  -- Initialize state
  let pcg : PCG := ⟨0, inc⟩
  -- Advance once and add seed
  let state1 := pcg.state * multiplier + inc
  let state2 := (state1 + s) * multiplier + inc
  ⟨state2, inc⟩

/-- Create a PCG generator with specific seed and stream.
    Different streams produce independent sequences. -/
def seedWithStream (s : UInt64) (stream : UInt64) : PCG :=
  -- Increment must be odd, so we use 2*stream + 1
  let inc := stream * 2 + 1
  let pcg : PCG := ⟨0, inc⟩
  let state1 := pcg.state * multiplier + inc
  let state2 := (state1 + s) * multiplier + inc
  ⟨state2, inc⟩

/-- Generate next random UInt32 using PCG-XSH-RR output function. -/
def nextUInt32 (pcg : PCG) : UInt32 × PCG :=
  let oldState := pcg.state
  -- Advance LCG
  let newState := oldState * multiplier + pcg.inc
  -- Calculate output using XSH-RR (xorshift high, random rotate)
  -- xorshifted = ((oldState >> 18) ^ oldState) >> 27
  let xorshifted := ((oldState >>> 18) ^^^ oldState) >>> 27
  -- rot = oldState >> 59 (top 5 bits)
  let rot := oldState >>> 59
  -- Rotate right: (xorshifted >> rot) | (xorshifted << (32 - rot))
  let rot32 := rot.toUInt32 &&& 31  -- 0-31
  let xorshifted32 := xorshifted.toUInt32
  let negRot32 := (32 - rot32) &&& 31  -- Handle rot=0 case
  let result := (xorshifted32 >>> rot32) ||| (xorshifted32 <<< negRot32)
  (result, ⟨newState, pcg.inc⟩)

/-- Generate next random UInt64 by combining two UInt32 outputs. -/
def nextUInt64 (pcg : PCG) : UInt64 × PCG :=
  let (lo, pcg1) := pcg.nextUInt32
  let (hi, pcg2) := pcg1.nextUInt32
  let result := (hi.toUInt64 <<< 32) ||| lo.toUInt64
  (result, pcg2)

/-- Generate random Float in [0, 1).
    Uses 32 bits of randomness for the mantissa. -/
def nextFloat (pcg : PCG) : Float × PCG :=
  let (u, pcg') := pcg.nextUInt32
  -- Divide by 2^32 to get [0, 1)
  let f := u.toNat.toFloat / 4294967296.0
  (f, pcg')

/-- Generate random Float in [0, 1) using full 64-bit precision. -/
def nextFloatFull (pcg : PCG) : Float × PCG :=
  let (u, pcg') := pcg.nextUInt64
  (u.toFloat / UInt64.toFloat (UInt64.ofNat 0xFFFFFFFFFFFFFFFF), pcg')

/-- Generate random Float in [lo, hi). -/
def nextFloatRange (pcg : PCG) (lo hi : Float) : Float × PCG :=
  let (f, pcg') := pcg.nextFloat
  (lo + f * (hi - lo), pcg')

/-- Generate two random Floats in [0, 1). -/
def nextFloat2 (pcg : PCG) : (Float × Float) × PCG :=
  let (f1, pcg1) := pcg.nextFloat
  let (f2, pcg2) := pcg1.nextFloat
  ((f1, f2), pcg2)

/-- Generate three random Floats in [0, 1). -/
def nextFloat3 (pcg : PCG) : (Float × Float × Float) × PCG :=
  let (f1, pcg1) := pcg.nextFloat
  let (f2, pcg2) := pcg1.nextFloat
  let (f3, pcg3) := pcg2.nextFloat
  ((f1, f2, f3), pcg3)

/-- Helper for bounded random generation with fuel. -/
private def nextBoundedGo (threshold bound : UInt32) (pcg : PCG) : Nat → UInt32 × PCG
  | 0 => (0, pcg)
  | fuel + 1 =>
    let (r, pcg') := pcg.nextUInt32
    if r >= threshold then (r % bound, pcg')
    else nextBoundedGo threshold bound pcg' fuel

/-- Generate random UInt32 in [0, bound).
    Uses rejection sampling to avoid modulo bias. -/
def nextBounded (pcg : PCG) (bound : UInt32) : UInt32 × PCG :=
  if bound == 0 then (0, pcg)
  else
    -- Calculate threshold for rejection to avoid bias
    -- threshold = -bound % bound = (2^32 - bound) % bound
    let threshold := (0 - bound) % bound
    -- Rejection sampling loop
    nextBoundedGo threshold bound pcg 100

/-- Generate random Int in [lo, hi] (inclusive). -/
def nextIntRange (pcg : PCG) (lo hi : Int) : Int × PCG :=
  if hi < lo then (lo, pcg)
  else
    let range := (hi - lo + 1).toNat
    let (r, pcg') := pcg.nextBounded range.toUInt32
    (lo + r.toNat, pcg')

/-- Generate a random boolean. -/
def nextBool (pcg : PCG) : Bool × PCG :=
  let (r, pcg') := pcg.nextUInt32
  (r &&& 1 == 1, pcg')

/-- Generate a random boolean with given probability of true. -/
def nextBoolWeighted (pcg : PCG) (probTrue : Float) : Bool × PCG :=
  let (f, pcg') := pcg.nextFloat
  (f < probTrue, pcg')

/-- Advance the generator by n steps (useful for parallel sequences).
    This is useful for creating independent subsequences for parallel processing. -/
def advance (pcg : PCG) (n : UInt64) : PCG :=
  if n == 0 then pcg
  else advanceLoop pcg n.toNat
where
  advanceLoop (pcg : PCG) : Nat → PCG
    | 0 => pcg
    | k + 1 => advanceLoop (pcg.nextUInt32.2) k

/-- Split the generator to create an independent stream.
    Useful for parallel random number generation. -/
def split (pcg : PCG) : PCG × PCG :=
  let (seed1, pcg1) := pcg.nextUInt64
  let (stream, pcg2) := pcg1.nextUInt64
  let newPcg := seedWithStream seed1 stream
  (newPcg, pcg2)

/-- Shuffle an array using Fisher-Yates algorithm. -/
def shuffle {α : Type} [Inhabited α] (pcg : PCG) (arr : Array α) : Array α × PCG := Id.run do
  if arr.size <= 1 then return (arr, pcg)
  let mut result := arr
  let mut rng := pcg
  for i in [0:arr.size - 1] do
    let remaining := arr.size - i
    let (j, rng') := rng.nextBounded remaining.toUInt32
    rng := rng'
    let swapIdx := i + j.toNat
    let tmp := result[i]!
    result := result.set! i result[swapIdx]!
    result := result.set! swapIdx tmp
  return (result, rng)

/-- Pick a random element from an array. Returns none if empty. -/
def choose {α : Type} [Inhabited α] (pcg : PCG) (arr : Array α) : Option α × PCG :=
  if arr.isEmpty then (none, pcg)
  else
    let (idx, pcg') := pcg.nextBounded arr.size.toUInt32
    (some arr[idx.toNat]!, pcg')

/-- Pick n random elements from an array (with replacement). -/
def sample {α : Type} [Inhabited α] (pcg : PCG) (arr : Array α) (n : Nat) : Array α × PCG := Id.run do
  if arr.isEmpty then return (#[], pcg)
  let mut result := #[]
  let mut rng := pcg
  for _ in [:n] do
    let (idx, rng') := rng.nextBounded arr.size.toUInt32
    rng := rng'
    result := result.push arr[idx.toNat]!
  return (result, rng)

/-- Pick n random elements from an array (without replacement).
    If n > arr.size, returns shuffled copy of arr. -/
def sampleWithoutReplacement {α : Type} [Inhabited α] (pcg : PCG) (arr : Array α) (n : Nat) : Array α × PCG :=
  let (shuffled, pcg') := shuffle pcg arr
  (shuffled.toList.take n |>.toArray, pcg')

end PCG

-- ============================================================================
-- Random Points in 2D Shapes (PCG versions)
-- ============================================================================

/-- Random point in unit circle using PCG. -/
def inUnitCirclePCG (pcg : PCG) : Vec2 × PCG := Id.run do
  let mut rng := pcg
  for _ in [:100] do
    let ((x, y), rng') := rng.nextFloat2
    rng := rng'
    let px := x * 2.0 - 1.0
    let py := y * 2.0 - 1.0
    if px * px + py * py <= 1.0 then
      return (Vec2.mk px py, rng)
  return (Vec2.zero, rng)

/-- Random point in circle using PCG. -/
def inCirclePCG (pcg : PCG) (center : Vec2) (radius : Float) : Vec2 × PCG :=
  let (p, pcg') := inUnitCirclePCG pcg
  (center + p.scale radius, pcg')

/-- Random point on unit circle circumference using PCG. -/
def onUnitCirclePCG (pcg : PCG) : Vec2 × PCG :=
  let (f, pcg') := pcg.nextFloat
  let angle := f * 2.0 * Float.pi
  (Vec2.mk (Float.cos angle) (Float.sin angle), pcg')

/-- Random point in unit sphere using PCG. -/
def inUnitSpherePCG (pcg : PCG) : Vec3 × PCG := Id.run do
  let mut rng := pcg
  for _ in [:100] do
    let ((x, y, z), rng') := rng.nextFloat3
    rng := rng'
    let px := x * 2.0 - 1.0
    let py := y * 2.0 - 1.0
    let pz := z * 2.0 - 1.0
    if px * px + py * py + pz * pz <= 1.0 then
      return (Vec3.mk px py pz, rng)
  return (Vec3.zero, rng)

/-- Random point on unit sphere surface using PCG. -/
def onUnitSpherePCG (pcg : PCG) : Vec3 × PCG :=
  let ((u, v), pcg') := pcg.nextFloat2
  let theta := u * 2.0 * Float.pi
  let phi := Float.acos (2.0 * v - 1.0)
  let sinPhi := Float.sin phi
  (Vec3.mk (sinPhi * Float.cos theta) (sinPhi * Float.sin theta) (Float.cos phi), pcg')

/-- Random point in axis-aligned box using PCG. -/
def inBoxPCG (pcg : PCG) (min max : Vec3) : Vec3 × PCG :=
  let ((fx, fy, fz), pcg') := pcg.nextFloat3
  let x := min.x + fx * (max.x - min.x)
  let y := min.y + fy * (max.y - min.y)
  let z := min.z + fz * (max.z - min.z)
  (Vec3.mk x y z, pcg')

/-- Generate a random unit vector (direction) using PCG. -/
def randomDirectionPCG (pcg : PCG) : Vec3 × PCG := onUnitSpherePCG pcg

/-- Generate a random rotation (as quaternion) using PCG.
    Uses uniform distribution over SO(3). -/
def randomRotationPCG (pcg : PCG) : Quat × PCG :=
  let ((u1, u2, u3), pcg') := pcg.nextFloat3
  -- Use Shoemake's method for uniform quaternion generation
  let sqrt1MinusU1 := Float.sqrt (1.0 - u1)
  let sqrtU1 := Float.sqrt u1
  let theta1 := 2.0 * Float.pi * u2
  let theta2 := 2.0 * Float.pi * u3
  let q := Quat.mk
    (sqrt1MinusU1 * Float.sin theta1)
    (sqrt1MinusU1 * Float.cos theta1)
    (sqrtU1 * Float.sin theta2)
    (sqrtU1 * Float.cos theta2)
  (q, pcg')

/-- Generate a random color (RGB, each channel in [0, 1]) using PCG. -/
def randomColorPCG (pcg : PCG) : (Float × Float × Float) × PCG :=
  pcg.nextFloat3

/-- Generate a random color with specified saturation/value ranges (HSV) using PCG.
    Returns RGB tuple. -/
def randomColorHSVPCG (pcg : PCG) (satMin satMax valMin valMax : Float) : (Float × Float × Float) × PCG :=
  let ((h, s, v), pcg') := pcg.nextFloat3
  let hue := h * 360.0
  let sat := satMin + s * (satMax - satMin)
  let val := valMin + v * (valMax - valMin)
  -- HSV to RGB conversion
  let c := val * sat
  -- Float modulo: a mod b = a - b * floor(a / b)
  let hueDiv60 := hue / 60.0
  let hueDiv60Mod2 := hueDiv60 - 2.0 * Float.floor (hueDiv60 / 2.0)
  let x := c * (1.0 - Float.abs (hueDiv60Mod2 - 1.0))
  let m := val - c
  let (r', g', b') :=
    if hue < 60.0 then (c, x, 0.0)
    else if hue < 120.0 then (x, c, 0.0)
    else if hue < 180.0 then (0.0, c, x)
    else if hue < 240.0 then (0.0, x, c)
    else if hue < 300.0 then (x, 0.0, c)
    else (c, 0.0, x)
  ((r' + m, g' + m, b' + m), pcg')

-- ============================================================================
-- Random Points in 2D Shapes
-- ============================================================================

/-- Random point in unit circle (radius 1, centered at origin).
    Uses rejection sampling for uniform distribution. -/
def inUnitCircle (rng : RNG) : Vec2 × RNG := Id.run do
  let mut rng := rng
  for _ in [:100] do  -- Limit iterations
    let ((x, y), rng') := rng.nextFloat2
    rng := rng'
    let px := x * 2.0 - 1.0
    let py := y * 2.0 - 1.0
    if px * px + py * py <= 1.0 then
      return (Vec2.mk px py, rng)
  -- Fallback (shouldn't happen often)
  return (Vec2.zero, rng)

/-- Random point in circle with given center and radius. -/
def inCircle (rng : RNG) (center : Vec2) (radius : Float) : Vec2 × RNG :=
  let (p, rng') := inUnitCircle rng
  (center + p.scale radius, rng')

/-- Random point on unit circle circumference. -/
def onUnitCircle (rng : RNG) : Vec2 × RNG :=
  let (f, rng') := rng.nextFloat
  let angle := f * 2.0 * Float.pi
  (Vec2.mk (Float.cos angle) (Float.sin angle), rng')

/-- Random point on circle circumference. -/
def onCircle (rng : RNG) (center : Vec2) (radius : Float) : Vec2 × RNG :=
  let (p, rng') := onUnitCircle rng
  (center + p.scale radius, rng')

/-- Random point in axis-aligned rectangle. -/
def inRectangle (rng : RNG) (min max : Vec2) : Vec2 × RNG :=
  let ((fx, fy), rng') := rng.nextFloat2
  let x := min.x + fx * (max.x - min.x)
  let y := min.y + fy * (max.y - min.y)
  (Vec2.mk x y, rng')

/-- Random point in triangle using barycentric coordinates. -/
def inTriangle2D (rng : RNG) (a b c : Vec2) : Vec2 × RNG :=
  let ((u, v), rng') := rng.nextFloat2
  -- Ensure point is in triangle (not parallelogram)
  let (u, v) := if u + v > 1.0 then (1.0 - u, 1.0 - v) else (u, v)
  let w := 1.0 - u - v
  let p := a.scale w + b.scale u + c.scale v
  (p, rng')

/-- Random point in ring (annulus) between inner and outer radius. -/
def inRing (rng : RNG) (center : Vec2) (innerRadius outerRadius : Float) : Vec2 × RNG :=
  let ((fr, fa), rng') := rng.nextFloat2
  -- Use sqrt for uniform distribution in area
  let r := Float.sqrt (innerRadius * innerRadius + fr * (outerRadius * outerRadius - innerRadius * innerRadius))
  let angle := fa * 2.0 * Float.pi
  let p := center + Vec2.mk (r * Float.cos angle) (r * Float.sin angle)
  (p, rng')

-- ============================================================================
-- Random Points in 3D Shapes
-- ============================================================================

/-- Random point in unit sphere (radius 1, centered at origin).
    Uses rejection sampling for uniform distribution. -/
def inUnitSphere (rng : RNG) : Vec3 × RNG := Id.run do
  let mut rng := rng
  for _ in [:100] do
    let ((x, y, z), rng') := rng.nextFloat3
    rng := rng'
    let px := x * 2.0 - 1.0
    let py := y * 2.0 - 1.0
    let pz := z * 2.0 - 1.0
    if px * px + py * py + pz * pz <= 1.0 then
      return (Vec3.mk px py pz, rng)
  return (Vec3.zero, rng)

/-- Random point in sphere with given center and radius. -/
def inSphere (rng : RNG) (center : Vec3) (radius : Float) : Vec3 × RNG :=
  let (p, rng') := inUnitSphere rng
  (center.add (p.scale radius), rng')

/-- Random point on unit sphere surface (uniform distribution). -/
def onUnitSphere (rng : RNG) : Vec3 × RNG :=
  let ((u, v), rng') := rng.nextFloat2
  let theta := u * 2.0 * Float.pi
  let phi := Float.acos (2.0 * v - 1.0)
  let sinPhi := Float.sin phi
  (Vec3.mk (sinPhi * Float.cos theta) (sinPhi * Float.sin theta) (Float.cos phi), rng')

/-- Random point on sphere surface. -/
def onSphere (rng : RNG) (center : Vec3) (radius : Float) : Vec3 × RNG :=
  let (p, rng') := onUnitSphere rng
  (center.add (p.scale radius), rng')

/-- Random point in axis-aligned box. -/
def inBox (rng : RNG) (min max : Vec3) : Vec3 × RNG :=
  let ((fx, fy, fz), rng') := rng.nextFloat3
  let x := min.x + fx * (max.x - min.x)
  let y := min.y + fy * (max.y - min.y)
  let z := min.z + fz * (max.z - min.z)
  (Vec3.mk x y z, rng')

/-- Random point in hemisphere (upper half of unit sphere, +Y up). -/
def inUnitHemisphere (rng : RNG) : Vec3 × RNG :=
  let (p, rng') := inUnitSphere rng
  let p := if p.y < 0.0 then Vec3.mk p.x (-p.y) p.z else p
  (p, rng')

/-- Random point on hemisphere surface (cosine-weighted, good for lighting). -/
def onHemisphereCosine (rng : RNG) (normal : Vec3) : Vec3 × RNG :=
  let ((u, v), rng') := rng.nextFloat2
  let r := Float.sqrt u
  let theta := 2.0 * Float.pi * v
  let x := r * Float.cos theta
  let y := r * Float.sin theta
  let z := Float.sqrt (1.0 - u)

  -- Build tangent space (simplified - assumes normal isn't (0,1,0))
  let tangent := if Float.abs' normal.y < 0.999
    then Vec3.unitY.cross normal |>.normalize
    else Vec3.unitX.cross normal |>.normalize
  let bitangent := normal.cross tangent

  -- Transform to world space
  let p := tangent.scale x |>.add (bitangent.scale y) |>.add (normal.scale z)
  (p, rng')

/-- Random point in cone (apex at origin, pointing in +Y direction).
    angle: half-angle of cone in radians
    height: height of cone -/
def inCone (rng : RNG) (angle height : Float) : Vec3 × RNG :=
  let ((fh, fr, fa), rng') := rng.nextFloat3
  -- Random height (cube root for uniform volume distribution)
  let h := Float.pow fh (1.0 / 3.0) * height
  -- Radius at this height
  let maxR := h * Float.tan angle
  -- Random radius (sqrt for uniform area distribution)
  let r := Float.sqrt fr * maxR
  -- Random angle
  let theta := fa * 2.0 * Float.pi
  (Vec3.mk (r * Float.cos theta) h (r * Float.sin theta), rng')

/-- Random point in cylinder (centered at origin, axis along Y).
    radius: cylinder radius
    height: full height of cylinder -/
def inCylinder (rng : RNG) (radius height : Float) : Vec3 × RNG :=
  let ((fr, fa, fh), rng') := rng.nextFloat3
  -- Random radius (sqrt for uniform area distribution)
  let r := Float.sqrt fr * radius
  let theta := fa * 2.0 * Float.pi
  let y := (fh - 0.5) * height
  (Vec3.mk (r * Float.cos theta) y (r * Float.sin theta), rng')

/-- Random point on cylinder surface (side only, not caps). -/
def onCylinderSide (rng : RNG) (radius height : Float) : Vec3 × RNG :=
  let ((fa, fh), rng') := rng.nextFloat2
  let theta := fa * 2.0 * Float.pi
  let y := (fh - 0.5) * height
  (Vec3.mk (radius * Float.cos theta) y (radius * Float.sin theta), rng')

/-- Random point in triangle (3D). -/
def inTriangle3D (rng : RNG) (a b c : Vec3) : Vec3 × RNG :=
  let ((u, v), rng') := rng.nextFloat2
  let (u, v) := if u + v > 1.0 then (1.0 - u, 1.0 - v) else (u, v)
  let w := 1.0 - u - v
  let p := a.scale w |>.add (b.scale u) |>.add (c.scale v)
  (p, rng')

/-- Random direction in a cone around the given direction.
    halfAngle: maximum angle from the direction in radians -/
def inConeDirection (rng : RNG) (direction : Vec3) (halfAngle : Float) : Vec3 × RNG :=
  let ((u, v), rng') := rng.nextFloat2
  -- Random angle within cone
  let theta := 2.0 * Float.pi * u
  let cosAngle := 1.0 - v * (1.0 - Float.cos halfAngle)
  let sinAngle := Float.sqrt (1.0 - cosAngle * cosAngle)

  -- Local direction
  let localX := sinAngle * Float.cos theta
  let localY := sinAngle * Float.sin theta
  let localZ := cosAngle

  -- Build tangent space
  let dir := direction.normalize
  let tangent := if Float.abs' dir.y < 0.999
    then Vec3.unitY.cross dir |>.normalize
    else Vec3.unitX.cross dir |>.normalize
  let bitangent := dir.cross tangent

  -- Transform to world space
  let p := tangent.scale localX |>.add (bitangent.scale localY) |>.add (dir.scale localZ)
  (p.normalize, rng')

end Random

end Linalg
