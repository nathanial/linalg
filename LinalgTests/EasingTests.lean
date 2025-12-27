/-
  Tests for Easing functions.
-/

import Linalg
import Crucible

namespace LinalgTests.EasingTests

open Crucible
open Linalg

testSuite "Easing Boundaries"

test "linear at 0 returns 0" := do
  ensure (floatNear (Easing.linear 0.0) 0.0 0.0001) "should be 0"

test "linear at 1 returns 1" := do
  ensure (floatNear (Easing.linear 1.0) 1.0 0.0001) "should be 1"

test "quadIn at 0 returns 0" := do
  ensure (floatNear (Easing.quadIn 0.0) 0.0 0.0001) "should be 0"

test "quadIn at 1 returns 1" := do
  ensure (floatNear (Easing.quadIn 1.0) 1.0 0.0001) "should be 1"

test "quadOut at 0 returns 0" := do
  ensure (floatNear (Easing.quadOut 0.0) 0.0 0.0001) "should be 0"

test "quadOut at 1 returns 1" := do
  ensure (floatNear (Easing.quadOut 1.0) 1.0 0.0001) "should be 1"

test "sineIn at 0 returns 0" := do
  ensure (floatNear (Easing.sineIn 0.0) 0.0 0.0001) "should be 0"

test "sineIn at 1 returns 1" := do
  ensure (floatNear (Easing.sineIn 1.0) 1.0 0.0001) "should be 1"

testSuite "Easing Midpoints"

test "linear at 0.5 returns 0.5" := do
  ensure (floatNear (Easing.linear 0.5) 0.5 0.0001) "should be 0.5"

test "quadInOut at 0.5 returns 0.5" := do
  ensure (floatNear (Easing.quadInOut 0.5) 0.5 0.0001) "should be 0.5"

test "cubicInOut at 0.5 returns 0.5" := do
  ensure (floatNear (Easing.cubicInOut 0.5) 0.5 0.0001) "should be 0.5"

test "sineInOut at 0.5 returns 0.5" := do
  ensure (floatNear (Easing.sineInOut 0.5) 0.5 0.0001) "should be 0.5"

testSuite "Easing Properties"

test "quadIn is slower than linear at start" := do
  let t := 0.3
  ensure (Easing.quadIn t < Easing.linear t) "quadIn should be slower at start"

test "quadOut is faster than linear at start" := do
  let t := 0.3
  ensure (Easing.quadOut t > Easing.linear t) "quadOut should be faster at start"

test "smoothStep at 0.5 returns 0.5" := do
  ensure (floatNear (Easing.smoothStep 0.5) 0.5 0.0001) "should be 0.5"

test "smootherStep at 0.5 returns 0.5" := do
  ensure (floatNear (Easing.smootherStep 0.5) 0.5 0.0001) "should be 0.5"

testSuite "Easing Apply"

test "apply interpolates between values" := do
  let result := Easing.apply Easing.linear 10.0 20.0 0.5
  ensure (floatNear result 15.0 0.0001) "should interpolate to 15"

test "apply with quadIn at midpoint" := do
  let result := Easing.apply Easing.quadIn 0.0 100.0 0.5
  ensure (floatNear result 25.0 0.0001) "quadIn at 0.5 should give 25"

#generate_tests

end LinalgTests.EasingTests
