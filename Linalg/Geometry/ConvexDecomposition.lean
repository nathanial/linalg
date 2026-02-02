/-
  Convex decomposition for triangle meshes (3D).

  This is a fast, conservative decomposition: triangles are clustered by
  centroid splits, and each cluster produces a convex hull proxy.
-/

import Linalg.Geometry.Mesh
import Linalg.Geometry.ConvexHull3D
import Linalg.Geometry.AABB

namespace Linalg

/-- A convex piece produced by mesh decomposition. -/
structure ConvexPiece where
  hull : ConvexHull3D
  bounds : AABB
  triangleIndices : Array Nat
deriving Repr, Inhabited

namespace ConvexPiece

/-- Number of triangles referenced by this piece. -/
@[inline]
def triangleCount (piece : ConvexPiece) : Nat :=
  piece.triangleIndices.size

/-- Collect triangles for this piece from the source mesh. -/
def triangles (piece : ConvexPiece) (mesh : Mesh) : Array Triangle :=
  piece.triangleIndices.foldl (fun acc idx =>
    match Mesh.triangle? mesh idx with
    | none => acc
    | some t => acc.push t
  ) #[]

end ConvexPiece

/-- Configuration for convex decomposition. -/
structure ConvexDecompositionConfig where
  /-- Max triangles per convex piece (0 = no limit). -/
  maxTrianglesPerPart : Nat := 64
  /-- Max recursion depth. -/
  maxDepth : Nat := 10
  /-- Minimum centroid extent to allow splitting. -/
  minSplitExtent : Float := 0.001
deriving Repr, BEq, Inhabited

namespace ConvexDecomposition

private structure TriInfo where
  index : Nat
  bounds : AABB
  centroid : Vec3
deriving Repr, Inhabited

private def triInfo (mesh : Mesh) (triIdx : Nat) : Option TriInfo :=
  match Mesh.triangle? mesh triIdx with
  | none => none
  | some tri =>
    let (minPt, maxPt) := tri.boundingBox
    let bounds := AABB.fromMinMax minPt maxPt
    some { index := triIdx, bounds := bounds, centroid := tri.centroid }

private def triInfos (mesh : Mesh) : Array TriInfo := Id.run do
  let count := mesh.triangleCount
  let mut infos : Array TriInfo := #[]
  for i in [:count] do
    match triInfo mesh i with
    | none => ()
    | some info => infos := infos.push info
  return infos

private def axisValue (v : Vec3) (axis : Fin 3) : Float :=
  match axis with
  | ⟨0, _⟩ => v.x
  | ⟨1, _⟩ => v.y
  | ⟨2, _⟩ => v.z

private def axisExtent (b : AABB) (axis : Fin 3) : Float :=
  match axis with
  | ⟨0, _⟩ => b.max.x - b.min.x
  | ⟨1, _⟩ => b.max.y - b.min.y
  | ⟨2, _⟩ => b.max.z - b.min.z

private def longestAxis (b : AABB) : Fin 3 :=
  let ex := b.max.x - b.min.x
  let ey := b.max.y - b.min.y
  let ez := b.max.z - b.min.z
  if ex >= ey && ex >= ez then ⟨0, by decide⟩
  else if ey >= ez then ⟨1, by decide⟩
  else ⟨2, by decide⟩

private def centroidBounds (prims : Array TriInfo) : AABB :=
  prims.foldl (fun acc p => AABB.expand acc p.centroid)
    (AABB.fromPoint prims[0]!.centroid)

private def mergeBounds (prims : Array TriInfo) : AABB :=
  prims.foldl (fun acc p => AABB.merge acc p.bounds) prims[0]!.bounds

private def collectVertexIndices (mesh : Mesh) (triIndices : Array Nat) : Array Nat := Id.run do
  let mut used := Array.replicate mesh.vertexCount false
  for triIdx in triIndices do
    let base := triIdx * 3
    if base + 2 < mesh.indices.size then
      let i0 := mesh.indices[base]!
      let i1 := mesh.indices[base + 1]!
      let i2 := mesh.indices[base + 2]!
      if i0 < used.size then used := used.set! i0 true
      if i1 < used.size then used := used.set! i1 true
      if i2 < used.size then used := used.set! i2 true
  let mut result : Array Nat := #[]
  for i in [:used.size] do
    if used[i]! then
      result := result.push i
  return result

private def collectPoints (mesh : Mesh) (triIndices : Array Nat) : Array Vec3 :=
  let indices := collectVertexIndices mesh triIndices
  indices.map (fun i => mesh.vertices[i]!)

private def makePiece (mesh : Mesh) (prims : Array TriInfo) : ConvexPiece :=
  let triangleIndices := prims.map (·.index)
  let points := collectPoints mesh triangleIndices
  let bounds :=
    if points.isEmpty then
      mergeBounds prims
    else
      points.foldl (fun acc p => AABB.expand acc p) (AABB.fromPoint points[0]!)
  let hull := ConvexHull3D.quickHull points
  { hull := hull, bounds := bounds, triangleIndices := triangleIndices }

private partial def buildParts (mesh : Mesh) (prims : Array TriInfo) (depth : Nat)
    (config : ConvexDecompositionConfig) : Array ConvexPiece :=
  if prims.isEmpty then
    #[]
  else
    let limit := if config.maxTrianglesPerPart == 0 then prims.size else config.maxTrianglesPerPart
    if depth >= config.maxDepth || prims.size <= limit then
      #[makePiece mesh prims]
    else
      let cBounds := centroidBounds prims
      let axis := longestAxis cBounds
      let extent := axisExtent cBounds axis
      if extent < config.minSplitExtent then
        #[makePiece mesh prims]
      else
        let sorted := prims.qsort fun a b => axisValue a.centroid axis < axisValue b.centroid axis
        let mid := sorted.size / 2
        let left := sorted.extract 0 mid
        let right := sorted.extract mid sorted.size
        if left.isEmpty || right.isEmpty then
          #[makePiece mesh prims]
        else
          buildParts mesh left (depth + 1) config ++ buildParts mesh right (depth + 1) config

/--
  Decompose a triangle mesh into convex pieces using a fast centroid split.
  Each piece is represented as a convex hull proxy over its triangle vertices.
-/
def decompose (mesh : Mesh) (config : ConvexDecompositionConfig := {}) : Array ConvexPiece :=
  if !mesh.isValid then
    #[]
  else
    let prims := triInfos mesh
    if prims.isEmpty then #[] else buildParts mesh prims 0 config

/-- Decompose and return only the convex hulls. -/
def decomposeHulls (mesh : Mesh) (config : ConvexDecompositionConfig := {}) : Array ConvexHull3D :=
  (decompose mesh config).map (·.hull)

end ConvexDecomposition

namespace Mesh

/-- Convenience wrapper for convex decomposition. -/
def convexDecompose (mesh : Mesh) (config : ConvexDecompositionConfig := {}) : Array ConvexPiece :=
  ConvexDecomposition.decompose mesh config

end Mesh

end Linalg
