/-
  Sampling utilities: low-discrepancy sequences and blue-noise sampling.
-/

import Init.Data.Nat.Bitwise
import Linalg.Core
import Linalg.Vec2
import Linalg.Vec3
import Linalg.Noise

namespace Linalg
namespace Sampling

/-- Radical inverse in the given base. -/
def radicalInverse (base : Nat) (index : Nat) : Float := Id.run do
  if base <= 1 then
    return 0.0
  let baseF := Float.ofNat base
  let mut i := index
  let mut f := 1.0 / baseF
  let mut result := 0.0
  while i > 0 do
    let digit := i % base
    result := result + Float.ofNat digit * f
    i := i / base
    f := f / baseF
  return result

/-- Halton sequence element (index starts at 0, skips the zero sample). -/
def halton (index : Nat) (base : Nat) : Float :=
  radicalInverse base (index + 1)

/-- Halton sequence in 2D. -/
def halton2 (index : Nat) (baseX : Nat := 2) (baseY : Nat := 3) : Vec2 :=
  Vec2.mk (halton index baseX) (halton index baseY)

/-- Halton sequence in 3D. -/
def halton3 (index : Nat) (baseX : Nat := 2) (baseY : Nat := 3) (baseZ : Nat := 5) : Vec3 :=
  Vec3.mk (halton index baseX) (halton index baseY) (halton index baseZ)

/-- Generate a Halton sequence in 2D. -/
def halton2Sequence (count : Nat) (baseX : Nat := 2) (baseY : Nat := 3) (start : Nat := 0) : Array Vec2 := Id.run do
  let mut result := #[]
  for i in [:count] do
    result := result.push (halton2 (start + i) baseX baseY)
  return result

/-- Hammersley point set in 2D. -/
def hammersley2 (index : Nat) (count : Nat) (base : Nat := 2) : Vec2 :=
  let u := if count == 0 then 0.0 else Float.ofNat index / Float.ofNat count
  let v := radicalInverse base index
  Vec2.mk u v

/-- Generate a Hammersley point set in 2D. -/
def hammersley2Sequence (count : Nat) (base : Nat := 2) : Array Vec2 := Id.run do
  let mut result := #[]
  for i in [:count] do
    result := result.push (hammersley2 i count base)
  return result

private def sobolBits : Nat := 30

private def sobolDirections (s : Nat) (a : Nat) (mInit : Array Nat) (bits : Nat := sobolBits) : Array Nat := Id.run do
  let mut m : Array Nat := #[]
  for i in [:bits] do
    if i < s then
      m := m.push (mInit[i]!)
    else
      let mut value := Nat.xor (m[i - s]!) ((m[i - s]!) >>> s)
      for k in [1:s] do
        let bitIndex := s - 1 - k
        if Nat.testBit a bitIndex then
          value := Nat.xor value ((m[i - k]!) >>> k)
      m := m.push value
  let mut v : Array Nat := #[]
  for i in [:bits] do
    v := v.push (m[i]! <<< (bits - i - 1))
  return v

private def sobolDir1 (bits : Nat := sobolBits) : Array Nat := Id.run do
  let mut v : Array Nat := #[]
  for i in [:bits] do
    v := v.push (1 <<< (bits - i - 1))
  return v

private def sobolDir2 (bits : Nat := sobolBits) : Array Nat :=
  sobolDirections 3 1 #[1, 3, 5] bits

private def sobolComponent (index : Nat) (dirs : Array Nat) (bits : Nat := sobolBits) : Float := Id.run do
  let mut x := 0
  let mut n := index
  let mut bit := 0
  while n > 0 && bit < dirs.size do
    if n % 2 == 1 then
      x := Nat.xor x (dirs[bit]!)
    n := n / 2
    bit := bit + 1
  let denom := Float.ofNat (1 <<< bits)
  return Float.ofNat x / denom

/-- Sobol sequence in 2D. -/
def sobol2 (index : Nat) : Vec2 :=
  let x := sobolComponent index (sobolDir1)
  let y := sobolComponent index (sobolDir2)
  Vec2.mk x y

/-- Generate a Sobol sequence in 2D. -/
def sobol2Sequence (count : Nat) (start : Nat := 0) : Array Vec2 := Id.run do
  let mut result := #[]
  for i in [:count] do
    result := result.push (sobol2 (start + i))
  return result

/-- Configuration for Poisson-disk (blue-noise) sampling. -/
structure PoissonConfig where
  radius : Float
  k : Nat := 30
  maxPoints : Nat := 0
deriving Repr, Inhabited

private def gridIndex (x y width : Nat) : Nat := y * width + x

private def cellCoord (min : Vec2) (cellInv : Float) (p : Vec2) : Nat × Nat :=
  let fx := (p.x - min.x) * cellInv
  let fy := (p.y - min.y) * cellInv
  let cx := fx |> Float.floor |> Float.toUInt64 |> UInt64.toNat
  let cy := fy |> Float.floor |> Float.toUInt64 |> UInt64.toNat
  (cx, cy)

/-- Poisson-disk sampling in a 2D axis-aligned box. -/
def poissonDisk2D (pcg : Random.PCG) (min max : Vec2) (config : PoissonConfig) : Array Vec2 × Random.PCG := Id.run do
  if config.radius <= 0.0 then
    return (#[], pcg)
  let size := max - min
  if size.x <= 0.0 || size.y <= 0.0 then
    return (#[], pcg)
  let cellSize := config.radius / Float.sqrt 2.0
  let cellInv := 1.0 / cellSize
  let gridW := Float.ceil (size.x / cellSize) |> Float.toUInt64 |> UInt64.toNat
  let gridH := Float.ceil (size.y / cellSize) |> Float.toUInt64 |> UInt64.toNat
  if gridW == 0 || gridH == 0 then
    return (#[], pcg)
  let mut grid : Array (Option Nat) := Array.replicate (gridW * gridH) none
  let ((fx, fy), pcg1) := pcg.nextFloat2
  let first := Vec2.mk (min.x + fx * size.x) (min.y + fy * size.y)
  let mut points : Array Vec2 := #[first]
  let mut active : Array Nat := #[0]
  let (cx0, cy0) := cellCoord min cellInv first
  grid := grid.set! (gridIndex cx0 cy0 gridW) (some 0)
  let mut rng := pcg1
  let radius2 := config.radius * config.radius
  let twoPi := Float.twoPi
  while !active.isEmpty do
    if config.maxPoints > 0 && points.size >= config.maxPoints then
      break
    let (pick, rng') := rng.nextBounded active.size.toUInt32
    rng := rng'
    let activeIdx := pick.toNat
    let baseIndex := active[activeIdx]!
    let basePoint := points[baseIndex]!
    let mut found := false
    let mut tries := 0
    while tries < config.k && !found do
      let (fa, rng1) := rng.nextFloat
      let (fr, rng2) := rng1.nextFloat
      rng := rng2
      let angle := fa * twoPi
      let radius := config.radius * (1.0 + fr)
      let offset := Vec2.mk (Float.cos angle * radius) (Float.sin angle * radius)
      let candidate := basePoint + offset
      if candidate.x >= min.x && candidate.x <= max.x && candidate.y >= min.y && candidate.y <= max.y then
        let (cx, cy) := cellCoord min cellInv candidate
        if cx < gridW && cy < gridH then
          let xStart := if cx >= 2 then cx - 2 else 0
          let yStart := if cy >= 2 then cy - 2 else 0
          let xEnd := Nat.min (cx + 3) gridW
          let yEnd := Nat.min (cy + 3) gridH
          let mut ok := true
          for y in [yStart:yEnd] do
            for x in [xStart:xEnd] do
              match grid[gridIndex x y gridW]! with
              | none => ()
              | some idx =>
                  let p := points[idx]!
                  if Vec2.distanceSquared p candidate < radius2 then
                    ok := false
          if ok then
            let newIndex := points.size
            points := points.push candidate
            active := active.push newIndex
            grid := grid.set! (gridIndex cx cy gridW) (some newIndex)
            found := true
      tries := tries + 1
    if !found then
      let lastIdx := active.size - 1
      if activeIdx < lastIdx then
        let last := active[lastIdx]!
        active := (active.set! activeIdx last).pop
      else
        active := active.pop
  return (points, rng)

/-- Alias for Poisson-disk sampling to emphasize blue-noise distribution. -/
def blueNoise2D (pcg : Random.PCG) (min max : Vec2) (config : PoissonConfig) : Array Vec2 × Random.PCG :=
  poissonDisk2D pcg min max config

end Sampling
end Linalg
