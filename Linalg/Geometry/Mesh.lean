/-
  Triangle mesh core types and utilities.
-/

import Linalg.Core
import Linalg.Vec2
import Linalg.Vec3
import Linalg.Geometry.Triangle
import Linalg.Geometry.AABB
import Linalg.Geometry.Ray
import Linalg.Geometry.Intersection
import Std.Data.HashMap

namespace Linalg

/-- Triangle mesh with vertex positions and triangle indices. -/
structure Mesh where
  vertices : Array Vec3
  indices : Array Nat
deriving Repr, Inhabited

/-- Ray hit on a mesh triangle. -/
structure MeshHit where
  triangle : Nat
  t : Float
  point : Vec3
  normal : Vec3
  u : Float
  v : Float
  w : Float
deriving Repr, Inhabited

/-- Adjacency info for a triangle (neighbors across edges 0,1,2). -/
structure TriangleAdjacency where
  e0 : Option Nat
  e1 : Option Nat
  e2 : Option Nat
deriving Repr, Inhabited

/-- Adjacency info for a mesh. -/
structure MeshAdjacency where
  faces : Array TriangleAdjacency
deriving Repr, Inhabited

/-- Half-edge structure for mesh connectivity. -/
structure HalfEdge where
  origin : Nat
  dest : Nat
  face : Nat
  next : Nat
  prev : Nat
  twin : Option Nat
deriving Repr, Inhabited

/-- Half-edge mesh representation. -/
structure HalfEdgeMesh where
  halfEdges : Array HalfEdge
  faceEdge : Array Nat
  vertexEdge : Array (Option Nat)
deriving Repr, Inhabited

namespace TriangleAdjacency

/-- Empty adjacency (no neighbors). -/
def empty : TriangleAdjacency := { e0 := none, e1 := none, e2 := none }

/-- Set neighbor across a given edge index. -/
def setEdge (adj : TriangleAdjacency) (edge : Nat) (value : Option Nat) : TriangleAdjacency :=
  match edge with
  | 0 => { adj with e0 := value }
  | 1 => { adj with e1 := value }
  | _ => { adj with e2 := value }

/-- Get neighbor across a given edge index. -/
def getEdge (adj : TriangleAdjacency) (edge : Nat) : Option Nat :=
  match edge with
  | 0 => adj.e0
  | 1 => adj.e1
  | _ => adj.e2

end TriangleAdjacency

namespace MeshAdjacency

/-- Create an empty adjacency structure with given face count. -/
def empty (faceCount : Nat) : MeshAdjacency :=
  { faces := Array.replicate faceCount TriangleAdjacency.empty }

/-- Get neighbor across an edge for a face. -/
def neighbor (adj : MeshAdjacency) (face edge : Nat) : Option Nat :=
  if face < adj.faces.size then
    TriangleAdjacency.getEdge adj.faces[face]! edge
  else none

/-- Get adjacency for a face. -/
def face (adj : MeshAdjacency) (faceIdx : Nat) : Option TriangleAdjacency :=
  if faceIdx < adj.faces.size then
    some adj.faces[faceIdx]!
  else none

end MeshAdjacency

namespace Mesh

/-- Create an empty mesh. -/
def empty : Mesh := { vertices := #[], indices := #[] }

/-- Create a mesh from vertices and indices. -/
def fromVerticesIndices (vertices : Array Vec3) (indices : Array Nat) : Mesh :=
  { vertices, indices }

/-- Number of vertices. -/
def vertexCount (m : Mesh) : Nat := m.vertices.size

/-- Number of indices. -/
def indexCount (m : Mesh) : Nat := m.indices.size

/-- Number of triangles (index count / 3). -/
def triangleCount (m : Mesh) : Nat := m.indices.size / 3

/-- Check if index count is a multiple of 3. -/
def hasValidIndexCount (m : Mesh) : Bool := m.indices.size % 3 == 0

/-- Check if all indices are within vertex range. -/
def indicesInRange (m : Mesh) : Bool :=
  m.indices.foldl (fun acc idx => acc && idx < m.vertices.size) true

/-- Check if mesh indices are valid. -/
def isValid (m : Mesh) : Bool := m.hasValidIndexCount && m.indicesInRange

/-- Get vertex at index (returns zero if out of range). -/
def vertex (m : Mesh) (i : Nat) : Vec3 :=
  if i < m.vertices.size then m.vertices[i]! else Vec3.zero

/-- Get triangle indices for a triangle index. -/
def triangleIndices? (m : Mesh) (triIdx : Nat) : Option (Nat × Nat × Nat) :=
  let base := triIdx * 3
  if base + 2 < m.indices.size then
    let i0 := m.indices[base]!
    let i1 := m.indices[base + 1]!
    let i2 := m.indices[base + 2]!
    if i0 < m.vertices.size && i1 < m.vertices.size && i2 < m.vertices.size then
      some (i0, i1, i2)
    else none
  else none

/-- Get triangle at index, if valid. -/
def triangle? (m : Mesh) (triIdx : Nat) : Option Triangle :=
  match m.triangleIndices? triIdx with
  | none => none
  | some (i0, i1, i2) =>
    some (Triangle.mk' m.vertices[i0]! m.vertices[i1]! m.vertices[i2]!)

/-- Get all triangles (returns none if mesh is invalid). -/
def triangles (m : Mesh) : Option (Array Triangle) := Id.run do
  if !m.isValid then
    return none
  let count := m.triangleCount
  let mut tris : Array Triangle := Array.mkEmpty count
  for i in [:count] do
    let base := i * 3
    let i0 := m.indices[base]!
    let i1 := m.indices[base + 1]!
    let i2 := m.indices[base + 2]!
    tris := tris.push (Triangle.mk' m.vertices[i0]! m.vertices[i1]! m.vertices[i2]!)
  return some tris

/-- Compute mesh bounds from all vertices. -/
def bounds (m : Mesh) : AABB :=
  if m.vertices.isEmpty then
    AABB.fromMinMax Vec3.zero Vec3.zero
  else
    m.vertices.foldl (fun b v => AABB.expand b v) (AABB.fromPoint m.vertices[0]!)

/-- Compute bounds of a triangle in the mesh. -/
def triangleBounds (m : Mesh) (triIdx : Nat) : Option AABB :=
  match m.triangle? triIdx with
  | none => none
  | some tri =>
    let (minPt, maxPt) := tri.boundingBox
    some (AABB.fromMinMax minPt maxPt)

/-- Compute triangle normals (unnormalized). -/
def triangleNormals (m : Mesh) : Array Vec3 := Id.run do
  let count := m.triangleCount
  let mut normals : Array Vec3 := Array.mkEmpty count
  if !m.isValid then
    for _ in [:count] do
      normals := normals.push Vec3.zero
    return normals
  for i in [:count] do
    let base := i * 3
    let i0 := m.indices[base]!
    let i1 := m.indices[base + 1]!
    let i2 := m.indices[base + 2]!
    let v0 := m.vertices[i0]!
    let v1 := m.vertices[i1]!
    let v2 := m.vertices[i2]!
    normals := normals.push ((v1.sub v0).cross (v2.sub v0))
  return normals

/-- Compute triangle unit normals. -/
def triangleUnitNormals (m : Mesh) : Array Vec3 :=
  m.triangleNormals.map (fun n => n.normalize)

/-- Compute per-vertex normals (area-weighted). -/
def vertexNormals (m : Mesh) : Array Vec3 := Id.run do
  let mut normals := Array.replicate m.vertexCount Vec3.zero
  if !m.isValid then
    return normals
  let triCount := m.triangleCount
  for i in [:triCount] do
    let base := i * 3
    let i0 := m.indices[base]!
    let i1 := m.indices[base + 1]!
    let i2 := m.indices[base + 2]!
    let v0 := m.vertices[i0]!
    let v1 := m.vertices[i1]!
    let v2 := m.vertices[i2]!
    let n := (v1.sub v0).cross (v2.sub v0)
    normals := normals.set! i0 (normals[i0]!.add n)
    normals := normals.set! i1 (normals[i1]!.add n)
    normals := normals.set! i2 (normals[i2]!.add n)
  return normals.map (fun n => n.normalize)

/-- Compute per-vertex tangents from positions and UVs. -/
def vertexTangents (m : Mesh) (uvs : Array Vec2) : Array Vec3 := Id.run do
  let mut tangents := Array.replicate m.vertexCount Vec3.zero
  if !m.isValid || uvs.size != m.vertexCount then
    return tangents
  let triCount := m.triangleCount
  for i in [:triCount] do
    let base := i * 3
    let i0 := m.indices[base]!
    let i1 := m.indices[base + 1]!
    let i2 := m.indices[base + 2]!
    let v0 := m.vertices[i0]!
    let v1 := m.vertices[i1]!
    let v2 := m.vertices[i2]!
    let uv0 := uvs[i0]!
    let uv1 := uvs[i1]!
    let uv2 := uvs[i2]!
    let edge1 := v1.sub v0
    let edge2 := v2.sub v0
    let duv1 := uv1.sub uv0
    let duv2 := uv2.sub uv0
    let det := duv1.x * duv2.y - duv1.y * duv2.x
    if Float.abs' det > Float.epsilon then
      let f := 1.0 / det
      let tangent := (edge1.scale duv2.y).sub (edge2.scale duv1.y) |>.scale f
      tangents := tangents.set! i0 (tangents[i0]!.add tangent)
      tangents := tangents.set! i1 (tangents[i1]!.add tangent)
      tangents := tangents.set! i2 (tangents[i2]!.add tangent)
  return tangents.map (fun t => t.normalize)

private structure UndirectedEdgeKey where
  a : Nat
  b : Nat
deriving Repr, BEq, Inhabited, Hashable

private def undirectedKey (i j : Nat) : UndirectedEdgeKey :=
  if i <= j then { a := i, b := j } else { a := j, b := i }

private structure DirectedEdgeKey where
  origin : Nat
  dest : Nat
deriving Repr, BEq, Inhabited, Hashable

private def linkTwin (halfEdges : Array HalfEdge)
    (edgeMap : Std.HashMap DirectedEdgeKey Nat) (origin dest idx : Nat) :
    Array HalfEdge × Std.HashMap DirectedEdgeKey Nat :=
  let key : DirectedEdgeKey := { origin := origin, dest := dest }
  let revKey : DirectedEdgeKey := { origin := dest, dest := origin }
  match edgeMap.get? revKey with
  | some twinIdx =>
    let he := halfEdges[idx]!
    let halfEdges := halfEdges.set! idx { he with twin := some twinIdx }
    let twin := halfEdges[twinIdx]!
    let halfEdges := halfEdges.set! twinIdx { twin with twin := some idx }
    (halfEdges, edgeMap.insert key idx)
  | none =>
    (halfEdges, edgeMap.insert key idx)

/-- Build triangle adjacency for the mesh. -/
def buildAdjacency (m : Mesh) : Option MeshAdjacency := Id.run do
  if !m.isValid then
    return none
  let faceCount := m.triangleCount
  let mut faces := Array.replicate faceCount TriangleAdjacency.empty
  let mut edgeMap : Std.HashMap UndirectedEdgeKey (Nat × Nat) := {}
  for f in [:faceCount] do
    let base := f * 3
    let i0 := m.indices[base]!
    let i1 := m.indices[base + 1]!
    let i2 := m.indices[base + 2]!
    let edges : Array (Nat × Nat × Nat) := #[(i0, i1, 0), (i1, i2, 1), (i2, i0, 2)]
    for e in edges do
      let (a, b, edgeIdx) := e
      let key := undirectedKey a b
      match edgeMap.get? key with
      | none =>
        edgeMap := edgeMap.insert key (f, edgeIdx)
      | some (otherFace, otherEdge) =>
        let adjOther := faces[otherFace]!
        faces := faces.set! otherFace (TriangleAdjacency.setEdge adjOther otherEdge (some f))
        let adjCurr := faces[f]!
        faces := faces.set! f (TriangleAdjacency.setEdge adjCurr edgeIdx (some otherFace))
  return some { faces := faces }

/-- Build half-edge connectivity for the mesh. -/
def buildHalfEdge (m : Mesh) : Option HalfEdgeMesh := Id.run do
  if !m.isValid then
    return none
  let faceCount := m.triangleCount
  let mut halfEdges : Array HalfEdge := #[]
  let mut faceEdge : Array Nat := Array.mkEmpty faceCount
  let mut vertexEdge : Array (Option Nat) := Array.replicate m.vertexCount none
  let mut edgeMap : Std.HashMap DirectedEdgeKey Nat := {}
  for f in [:faceCount] do
    let base := f * 3
    let i0 := m.indices[base]!
    let i1 := m.indices[base + 1]!
    let i2 := m.indices[base + 2]!
    let startIdx := halfEdges.size
    let he0 : HalfEdge := ⟨i0, i1, f, startIdx + 1, startIdx + 2, none⟩
    let he1 : HalfEdge := ⟨i1, i2, f, startIdx + 2, startIdx, none⟩
    let he2 : HalfEdge := ⟨i2, i0, f, startIdx, startIdx + 1, none⟩
    halfEdges := halfEdges.push he0
    halfEdges := halfEdges.push he1
    halfEdges := halfEdges.push he2
    faceEdge := faceEdge.push startIdx
    match vertexEdge[i0]! with
    | none => vertexEdge := vertexEdge.set! i0 (some startIdx)
    | some _ => ()
    match vertexEdge[i1]! with
    | none => vertexEdge := vertexEdge.set! i1 (some (startIdx + 1))
    | some _ => ()
    match vertexEdge[i2]! with
    | none => vertexEdge := vertexEdge.set! i2 (some (startIdx + 2))
    | some _ => ()
    let (halfEdges', edgeMap') := linkTwin halfEdges edgeMap i0 i1 startIdx
    halfEdges := halfEdges'
    edgeMap := edgeMap'
    let (halfEdges', edgeMap') := linkTwin halfEdges edgeMap i1 i2 (startIdx + 1)
    halfEdges := halfEdges'
    edgeMap := edgeMap'
    let (halfEdges', edgeMap') := linkTwin halfEdges edgeMap i2 i0 (startIdx + 2)
    halfEdges := halfEdges'
    edgeMap := edgeMap'
  return some { halfEdges := halfEdges, faceEdge := faceEdge, vertexEdge := vertexEdge }

/-- Ray cast against all mesh triangles, returning closest hit. -/
def rayCast (m : Mesh) (ray : Ray) (cullBackface : Bool := false) : Option MeshHit := Id.run do
  if !m.isValid then
    return none
  match Intersection.rayAABB ray m.bounds with
  | none => return none
  | some _ => ()
  let mut best : Option MeshHit := none
  let triCount := m.triangleCount
  for i in [:triCount] do
    let base := i * 3
    let i0 := m.indices[base]!
    let i1 := m.indices[base + 1]!
    let i2 := m.indices[base + 2]!
    let tri := Triangle.mk' m.vertices[i0]! m.vertices[i1]! m.vertices[i2]!
    match Intersection.rayTriangleBarycentric ray tri cullBackface with
    | none => ()
    | some hit =>
      let maxT := match best with
        | some b => b.t
        | none => Float.infinity
      if hit.t < maxT then
        best := some ⟨i, hit.t, hit.point, hit.normal, hit.u, hit.v, hit.w⟩
  return best

/-- Ray cast against all mesh triangles, returning all hits. -/
def rayCastAll (m : Mesh) (ray : Ray) (cullBackface : Bool := false) : Array MeshHit := Id.run do
  let mut hits : Array MeshHit := #[]
  if !m.isValid then
    return hits
  match Intersection.rayAABB ray m.bounds with
  | none => return hits
  | some _ => ()
  let triCount := m.triangleCount
  for i in [:triCount] do
    let base := i * 3
    let i0 := m.indices[base]!
    let i1 := m.indices[base + 1]!
    let i2 := m.indices[base + 2]!
    let tri := Triangle.mk' m.vertices[i0]! m.vertices[i1]! m.vertices[i2]!
    match Intersection.rayTriangleBarycentric ray tri cullBackface with
    | none => ()
    | some hit =>
      hits := hits.push ⟨i, hit.t, hit.point, hit.normal, hit.u, hit.v, hit.w⟩
  return hits

/-- Check if any triangle is hit by the ray within maxT. -/
def rayAny (m : Mesh) (ray : Ray) (maxT : Float := Float.infinity)
    (cullBackface : Bool := false) : Bool := Id.run do
  if !m.isValid then
    return false
  match Intersection.rayAABB ray m.bounds with
  | none => return false
  | some _ => ()
  let triCount := m.triangleCount
  let mut hit := false
  for i in [:triCount] do
    if !hit then
      let base := i * 3
      let i0 := m.indices[base]!
      let i1 := m.indices[base + 1]!
      let i2 := m.indices[base + 2]!
      let tri := Triangle.mk' m.vertices[i0]! m.vertices[i1]! m.vertices[i2]!
      match Intersection.rayTriangle ray tri cullBackface with
      | none => ()
      | some info =>
        if info.t <= maxT then
          hit := true
  return hit

end Mesh

end Linalg
