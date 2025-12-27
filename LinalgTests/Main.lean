/-
  Linalg Test Suite
  Main entry point for running all tests.
-/

import LinalgTests.Vec2Tests
import LinalgTests.Vec3Tests
import LinalgTests.Vec4Tests
import LinalgTests.Mat2Tests
import LinalgTests.Mat3Tests
import LinalgTests.Mat4Tests
import LinalgTests.QuatTests
import LinalgTests.GeometryTests
import LinalgTests.TriangleTests
import LinalgTests.EasingTests
import LinalgTests.TransformTests
import LinalgTests.FrustumTests
import LinalgTests.CurveTests
import LinalgTests.OBBTests
import LinalgTests.CapsuleTests
import Crucible

open Crucible

def main : IO UInt32 := do
  IO.println "Linalg Game Math Library Tests"
  IO.println "==============================="
  IO.println ""

  let result ← runAllSuites

  IO.println ""
  IO.println "==============================="

  if result != 0 then
    IO.println "Some tests failed!"
    return 1
  else
    IO.println "All tests passed!"
    return 0
