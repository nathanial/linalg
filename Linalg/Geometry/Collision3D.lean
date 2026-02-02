/-
  3D Convex Collision Detection

  Implements:
  - GJK (Gilbert-Johnson-Keerthi) for general convex shape intersection
  - EPA (Expanding Polytope Algorithm) for penetration depth / MTV
  - MPR (Minkowski Portal Refinement) interface (delegates to GJK/EPA)
-/

import Linalg.Core
import Linalg.Vec3
import Linalg.Geometry.Sphere
import Linalg.Geometry.AABB
import Linalg.Geometry.OBB
import Linalg.Geometry.Capsule
import Linalg.Geometry.Triangle

namespace Linalg

-- ============================================================================
-- Collision Result Types
-- ============================================================================

/-- Result of a 3D collision test with penetration information. -/
structure CollisionResult3D where
  /-- Whether the shapes are colliding. -/
  colliding : Bool
  /-- Minimum Translation Vector direction (unit).
      Points from shape A toward shape B. -/
  mtv : Vec3
  /-- Penetration depth (magnitude of MTV). -/
  depth : Float
deriving Repr, BEq, Inhabited

namespace CollisionResult3D

/-- No collision result. -/
def none : CollisionResult3D := ⟨false, Vec3.zero, 0.0⟩

/-- Create a collision result with MTV. -/
def collision (mtv : Vec3) (depth : Float) : CollisionResult3D :=
  ⟨true, mtv, depth⟩

end CollisionResult3D

-- ============================================================================
-- Support Function Typeclass (for GJK/EPA/MPR)
-- ============================================================================

/-- Typeclass for shapes that can provide a support point in 3D. -/
class Support3D (α : Type) where
  /-- Get the support point in the given direction. -/
  support : α → Vec3 → Vec3

-- ============================================================================
-- GJK (Gilbert-Johnson-Keerthi Algorithm)
-- ============================================================================

namespace Collision3D

namespace GJK

/-- Support point in Minkowski difference with original points tracked. -/
structure MinkowskiPoint where
  /-- Point in Minkowski difference space (a - b). -/
  point : Vec3
  /-- Support point from shape A. -/
  supportA : Vec3
  /-- Support point from shape B. -/
  supportB : Vec3
deriving Repr, Inhabited

/-- Get support point in Minkowski difference A - B. -/
def minkowskiSupport {α β : Type} [Support3D α] [Support3D β] (a : α) (b : β) (dir : Vec3) : MinkowskiPoint :=
  let sa := Support3D.support a dir
  let sb := Support3D.support b (dir.scale (-1.0))
  ⟨sa - sb, sa, sb⟩

/-- Pick a perpendicular direction to a vector (non-zero if input is non-zero). -/
private def perpendicular (v : Vec3) : Vec3 :=
  let ax := Float.abs' v.x
  let ay := Float.abs' v.y
  let az := Float.abs' v.z
  let axis :=
    if ax < ay then
      if ax < az then Vec3.unitX else Vec3.unitZ
    else
      if ay < az then Vec3.unitY else Vec3.unitZ
  v.cross axis

/-- Process a 2-point simplex (line). -/
private def processLine (simplex : Array MinkowskiPoint) : Bool × Array MinkowskiPoint × Vec3 :=
  let a := simplex[1]!.point  -- Most recent point
  let b := simplex[0]!.point
  let ab := b - a
  let ao := Vec3.zero - a
  if ab.dot ao > 0.0 then
    let dir := (ab.cross ao).cross ab
    let dir := if dir.lengthSquared > Float.epsilon then dir else perpendicular ab
    (false, simplex, dir)
  else
    (false, #[simplex[1]!], ao)

/-- Process a 3-point simplex (triangle). -/
private def processTriangle (simplex : Array MinkowskiPoint) : Bool × Array MinkowskiPoint × Vec3 :=
  let a := simplex[2]!.point  -- Most recent point
  let b := simplex[1]!.point
  let c := simplex[0]!.point
  let ab := b - a
  let ac := c - a
  let ao := Vec3.zero - a
  let abc := ab.cross ac

  let abPerp := ab.cross abc
  if abPerp.dot ao > 0.0 then
    if ab.dot ao > 0.0 then
      let dir := (ab.cross ao).cross ab
      (false, #[simplex[1]!, simplex[2]!], dir)
    else
      (false, #[simplex[2]!], ao)
  else
    let acPerp := abc.cross ac
    if acPerp.dot ao > 0.0 then
      if ac.dot ao > 0.0 then
        let dir := (ac.cross ao).cross ac
        (false, #[simplex[0]!, simplex[2]!], dir)
      else
        (false, #[simplex[2]!], ao)
    else
      let dir := if abc.dot ao > 0.0 then abc else abc.scale (-1.0)
      (false, simplex, dir)

/-- Compute outward normal of face (a,b,c) relative to opposite point d. -/
private def faceNormalOut (a b c d : Vec3) : Vec3 :=
  let n := (b - a).cross (c - a)
  if n.dot (d - a) > 0.0 then n.scale (-1.0) else n

/-- Process a 4-point simplex (tetrahedron). -/
private def processTetrahedron (simplex : Array MinkowskiPoint) : Bool × Array MinkowskiPoint × Vec3 :=
  let a := simplex[3]!.point  -- Most recent point
  let b := simplex[2]!.point
  let c := simplex[1]!.point
  let d := simplex[0]!.point
  let ao := Vec3.zero - a

  let nABC := faceNormalOut a b c d
  if nABC.dot ao > 0.0 then
    (false, #[simplex[1]!, simplex[2]!, simplex[3]!], nABC)
  else
    let nACD := faceNormalOut a c d b
    if nACD.dot ao > 0.0 then
      (false, #[simplex[0]!, simplex[1]!, simplex[3]!], nACD)
    else
      let nADB := faceNormalOut a d b c
      if nADB.dot ao > 0.0 then
        (false, #[simplex[0]!, simplex[2]!, simplex[3]!], nADB)
      else
        let nBCD := faceNormalOut b c d a
        let bo := Vec3.zero - b
        if nBCD.dot bo > 0.0 then
          (false, #[simplex[0]!, simplex[1]!, simplex[2]!], nBCD)
        else
          (true, simplex, Vec3.zero)

/-- Check if simplex contains the origin and update search direction.
    Returns (containsOrigin, newSimplex, newDirection). -/
def processSimplex (simplex : Array MinkowskiPoint) : Bool × Array MinkowskiPoint × Vec3 :=
  match simplex.size with
  | 2 => processLine simplex
  | 3 => processTriangle simplex
  | 4 => processTetrahedron simplex
  | _ => (false, simplex, Vec3.zero)

/-- Maximum GJK iterations to prevent infinite loops. -/
def maxIterations : Nat := 32

/-- GJK loop for intersection test. -/
private partial def gjkLoop {α β : Type} [Support3D α] [Support3D β]
    (a : α) (b : β) (simplex : Array MinkowskiPoint) (dir : Vec3) (iter : Nat) : Bool :=
  if iter >= maxIterations then false
  else if dir.lengthSquared < Float.epsilon then true
  else
    let support := minkowskiSupport a b dir
    if support.point.dot dir < 0.0 then
      false
    else
      let newSimplex := simplex.push support
      let (containsOrigin, simplex', newDir) := processSimplex newSimplex
      if containsOrigin then true
      else gjkLoop a b simplex' newDir (iter + 1)

/-- Run GJK algorithm to test if two convex shapes intersect. -/
def intersects {α β : Type} [Support3D α] [Support3D β] (a : α) (b : β) : Bool :=
  let initialDir := Vec3.unitX
  let support1 := minkowskiSupport a b initialDir
  if support1.point.dot initialDir < 0.0 then
    false
  else
    let dir := Vec3.zero - support1.point
    gjkLoop a b #[support1] dir 0

/-- GJK loop that returns the final simplex. -/
private partial def gjkLoopWithSimplex {α β : Type} [Support3D α] [Support3D β]
    (a : α) (b : β) (simplex : Array MinkowskiPoint) (dir : Vec3) (iter : Nat) : Option (Array MinkowskiPoint) :=
  if iter >= maxIterations then none
  else if dir.lengthSquared < Float.epsilon then some simplex
  else
    let support := minkowskiSupport a b dir
    if support.point.dot dir < 0.0 then none
    else
      let newSimplex := simplex.push support
      let (containsOrigin, simplex', newDir) := processSimplex newSimplex
      if containsOrigin then some simplex'
      else gjkLoopWithSimplex a b simplex' newDir (iter + 1)

/-- Direction from A toward B using support points in a simplex. -/
private def simplexDirAB (simplex : Array MinkowskiPoint) : Vec3 :=
  let (sumA, sumB) := simplex.foldl (init := (Vec3.zero, Vec3.zero)) (fun (accA, accB) p =>
    (accA + p.supportA, accB + p.supportB)
  )
  sumB - sumA

/-- Ensure normal points from A toward B. -/
private def orientNormal (normal : Vec3) (dirAB : Vec3) : Vec3 :=
  if normal.dot dirAB >= 0.0 then normal else normal.scale (-1.0)

/-- EPA face representation. -/
private structure Face where
  i1 : Nat
  i2 : Nat
  i3 : Nat
  normal : Vec3
  distance : Float
deriving Repr, Inhabited

/-- Build a face with outward normal (away from origin). -/
private def makeFace (verts : Array MinkowskiPoint) (i1 i2 i3 : Nat) : Face :=
  let a := verts[i1]!.point
  let b := verts[i2]!.point
  let c := verts[i3]!.point
  let n := (b - a).cross (c - a)
  let nLenSq := n.lengthSquared
  if nLenSq < Float.epsilon then
    let normal := Vec3.unitX
    let dist := Float.abs' (normal.dot a)
    { i1 := i1, i2 := i2, i3 := i3, normal := normal, distance := dist }
  else
    let normal := n.normalize
    let dist := normal.dot a
    if dist < 0.0 then
      { i1 := i1, i2 := i3, i3 := i2, normal := normal.scale (-1.0), distance := -dist }
    else
      { i1 := i1, i2 := i2, i3 := i3, normal := normal, distance := dist }

/-- Remove an element at index from array. -/
private def eraseAt (arr : Array (Nat × Nat)) (idx : Nat) : Array (Nat × Nat) :=
  arr.foldl (init := (#[], 0)) (fun (result, i) elem =>
    if i == idx then (result, i + 1)
    else (result.push elem, i + 1)
  ) |>.1

/-- Find index of reverse edge if present. -/
private def findReverseEdgeIdx (edges : Array (Nat × Nat)) (edge : Nat × Nat) : Option Nat :=
  edges.foldl (init := (none, 0)) (fun (found, i) elem =>
    match found with
    | some _ => (found, i + 1)
    | none =>
      let (a, b) := elem
      if a == edge.2 && b == edge.1 then (some i, i + 1) else (none, i + 1)
  ) |>.1

/-- Add edge to boundary list or remove reverse edge if it exists. -/
private def addOrRemoveEdge (edges : Array (Nat × Nat)) (edge : Nat × Nat) : Array (Nat × Nat) :=
  match findReverseEdgeIdx edges edge with
  | some idx => eraseAt edges idx
  | none => edges.push edge

/-- Direction from A toward B for a face. -/
private def faceDirAB (verts : Array MinkowskiPoint) (face : Face) : Vec3 :=
  let a1 := verts[face.i1]!.supportA
  let a2 := verts[face.i2]!.supportA
  let a3 := verts[face.i3]!.supportA
  let b1 := verts[face.i1]!.supportB
  let b2 := verts[face.i2]!.supportB
  let b3 := verts[face.i3]!.supportB
  (b1 + b2 + b3) - (a1 + a2 + a3)

/-- EPA iterations to get penetration depth and MTV. -/
private def epa {α β : Type} [Support3D α] [Support3D β]
    (a : α) (b : β) (simplex : Array MinkowskiPoint) : CollisionResult3D := Id.run do
  if simplex.size < 4 then
    let dirAB := simplexDirAB simplex
    let normal :=
      if simplex.size == 3 then
        let p0 := simplex[0]!.point
        let p1 := simplex[1]!.point
        let p2 := simplex[2]!.point
        let n := (p1 - p0).cross (p2 - p0)
        if n.lengthSquared > Float.epsilon then n.normalize else Vec3.unitX
      else if simplex.size == 2 then
        let d := simplex[1]!.point - simplex[0]!.point
        let n := perpendicular d
        if n.lengthSquared > Float.epsilon then n.normalize else Vec3.unitX
      else
        Vec3.unitX
    let oriented := orientNormal normal dirAB
    return CollisionResult3D.collision oriented Float.epsilon

  let mut vertices := simplex
  let mut faces : Array Face := #[
    makeFace vertices 0 1 2,
    makeFace vertices 0 3 1,
    makeFace vertices 0 2 3,
    makeFace vertices 1 3 2
  ]

  let maxIterations := 64
  let tolerance := 1.0e-4

  for _ in [:maxIterations] do
    if faces.isEmpty then
      return CollisionResult3D.collision Vec3.unitX Float.epsilon

    let mut bestIdx := 0
    let mut bestDist := faces[0]!.distance
    for i in [1:faces.size] do
      let dist := faces[i]!.distance
      if dist < bestDist then
        bestDist := dist
        bestIdx := i

    let bestFace := faces[bestIdx]!
    let support := minkowskiSupport a b bestFace.normal
    let supportDist := bestFace.normal.dot support.point

    if supportDist - bestFace.distance <= tolerance then
      let dirAB := faceDirAB vertices bestFace
      let normal := orientNormal bestFace.normal dirAB
      return CollisionResult3D.collision normal bestFace.distance

    let newIndex := vertices.size
    vertices := vertices.push support

    let mut newFaces : Array Face := #[]
    let mut edges : Array (Nat × Nat) := #[]
    for face in faces do
      let aPoint := vertices[face.i1]!.point
      if face.normal.dot (support.point - aPoint) > 0.0 then
        edges := addOrRemoveEdge edges (face.i1, face.i2)
        edges := addOrRemoveEdge edges (face.i2, face.i3)
        edges := addOrRemoveEdge edges (face.i3, face.i1)
      else
        newFaces := newFaces.push face

    for edge in edges do
      let (i1, i2) := edge
      newFaces := newFaces.push (makeFace vertices i1 i2 newIndex)

    faces := newFaces

  let fallbackFace := faces.getD 0 (makeFace vertices 0 1 2)
  let dirAB := faceDirAB vertices fallbackFace
  let normal := orientNormal fallbackFace.normal dirAB
  CollisionResult3D.collision normal fallbackFace.distance

/-- GJK with collision result (includes EPA for MTV). -/
def collision {α β : Type} [Support3D α] [Support3D β] (a : α) (b : β) : CollisionResult3D :=
  let initialDir := Vec3.unitX
  let support1 := minkowskiSupport a b initialDir
  if support1.point.dot initialDir < 0.0 then
    CollisionResult3D.none
  else
    let dir := Vec3.zero - support1.point
    match gjkLoopWithSimplex a b #[support1] dir 0 with
    | none => CollisionResult3D.none
    | some simplex => epa a b simplex

end GJK

-- ============================================================================
-- MPR (Minkowski Portal Refinement)
-- ============================================================================

namespace MPR

/-- MPR intersection test (delegates to GJK for robustness). -/
def intersects {α β : Type} [Support3D α] [Support3D β] (a : α) (b : β) : Bool :=
  GJK.intersects a b

/-- MPR collision result (delegates to GJK/EPA). -/
def collision {α β : Type} [Support3D α] [Support3D β] (a : α) (b : β) : CollisionResult3D :=
  GJK.collision a b

end MPR

end Collision3D

-- ============================================================================
-- Support3D Instances
-- ============================================================================

instance : Support3D Sphere where
  support sphere dir :=
    let normalizedDir := if dir.lengthSquared > Float.epsilon
                         then dir.normalize
                         else Vec3.unitX
    sphere.center + normalizedDir.scale sphere.radius

instance : Support3D AABB where
  support aabb dir :=
    Vec3.mk
      (if dir.x >= 0.0 then aabb.max.x else aabb.min.x)
      (if dir.y >= 0.0 then aabb.max.y else aabb.min.y)
      (if dir.z >= 0.0 then aabb.max.z else aabb.min.z)

instance : Support3D OBB where
  support obb dir :=
    let localDir := obb.orientation.inverse.rotateVec3 dir
    let localSupport := Vec3.mk
      (if localDir.x >= 0.0 then obb.halfExtents.x else -obb.halfExtents.x)
      (if localDir.y >= 0.0 then obb.halfExtents.y else -obb.halfExtents.y)
      (if localDir.z >= 0.0 then obb.halfExtents.z else -obb.halfExtents.z)
    obb.localToWorld localSupport

instance : Support3D Capsule where
  support capsule dir :=
    let base := if dir.dot capsule.a >= dir.dot capsule.b then capsule.a else capsule.b
    let normalizedDir := if dir.lengthSquared > Float.epsilon
                         then dir.normalize
                         else Vec3.unitX
    base + normalizedDir.scale capsule.radius

instance : Support3D Triangle where
  support tri dir :=
    let d0 := dir.dot tri.v0
    let d1 := dir.dot tri.v1
    let d2 := dir.dot tri.v2
    if d0 >= d1 && d0 >= d2 then tri.v0
    else if d1 >= d2 then tri.v1
    else tri.v2

-- ============================================================================
-- Convenience Functions
-- ============================================================================

/-- Test if two convex shapes collide using GJK/EPA (3D). -/
def collideGJK3D {α β : Type} [Support3D α] [Support3D β] (a : α) (b : β) : CollisionResult3D :=
  Collision3D.GJK.collision a b

/-- Quick intersection test using GJK (3D). -/
def intersectsGJK3D {α β : Type} [Support3D α] [Support3D β] (a : α) (b : β) : Bool :=
  Collision3D.GJK.intersects a b

/-- Test if two convex shapes collide using MPR (3D). -/
def collideMPR3D {α β : Type} [Support3D α] [Support3D β] (a : α) (b : β) : CollisionResult3D :=
  Collision3D.MPR.collision a b

/-- Quick intersection test using MPR (3D). -/
def intersectsMPR3D {α β : Type} [Support3D α] [Support3D β] (a : α) (b : β) : Bool :=
  Collision3D.MPR.intersects a b

end Linalg
