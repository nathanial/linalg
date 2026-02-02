/-
  Tests for sampling utilities.
-/

import Linalg
import Crucible

namespace LinalgTests.SamplingTests

open Crucible
open Linalg
open Linalg.Sampling

testSuite "Sampling Sequences"

test "Halton2 first sample" := do
  let p := halton2 0 2 3
  ensure (floatNear p.x 0.5 0.0001) "halton x"
  ensure (floatNear p.y (1.0 / 3.0) 0.0001) "halton y"

test "Hammersley2 first samples" := do
  let p0 := hammersley2 0 4 2
  let p1 := hammersley2 1 4 2
  ensure (floatNear p0.x 0.0 0.0001) "hammersley x0"
  ensure (floatNear p0.y 0.0 0.0001) "hammersley y0"
  ensure (floatNear p1.x 0.25 0.0001) "hammersley x1"
  ensure (floatNear p1.y 0.5 0.0001) "hammersley y1"

test "Sobol2 first samples" := do
  let p0 := sobol2 0
  let p1 := sobol2 1
  let p2 := sobol2 2
  let p3 := sobol2 3
  ensure (floatNear p0.x 0.0 0.0001) "sobol x0"
  ensure (floatNear p0.y 0.0 0.0001) "sobol y0"
  ensure (floatNear p1.x 0.5 0.0001) "sobol x1"
  ensure (floatNear p1.y 0.5 0.0001) "sobol y1"
  ensure (floatNear p2.x 0.25 0.0001) "sobol x2"
  ensure (floatNear p2.y 0.75 0.0001) "sobol y2"
  ensure (floatNear p3.x 0.75 0.0001) "sobol x3"
  ensure (floatNear p3.y 0.25 0.0001) "sobol y3"

testSuite "Sampling Blue Noise"

test "Poisson disk separation" := do
  let pcg := Random.PCG.seed 12345
  let config : PoissonConfig := { radius := 0.15, k := 20, maxPoints := 64 }
  let (points, _) := poissonDisk2D pcg (Vec2.mk 0.0 0.0) (Vec2.mk 1.0 1.0) config
  ensure (points.size > 0) "should generate points"
  let r2 := config.radius * config.radius - 0.000001
  let mut ok := true
  for i in [:points.size] do
    let p := points[i]!
    if p.x < 0.0 || p.x > 1.0 || p.y < 0.0 || p.y > 1.0 then
      ok := false
    for j in [i+1:points.size] do
      let d2 := Vec2.distanceSquared p points[j]!
      if d2 < r2 then
        ok := false
  ensure ok "points should be separated and within bounds"

end LinalgTests.SamplingTests
