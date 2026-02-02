/-
  BVH helpers for triangle meshes.
-/

import Linalg.Geometry.Mesh
import Linalg.Spatial.BVH

namespace Linalg.Spatial

open Linalg

namespace MeshBVH

/-- Build a BVH from mesh triangles. -/
def build (mesh : Mesh) (config : BVHConfig := {}) : Option BVH :=
  match Mesh.triangles mesh with
  | none => none
  | some tris =>
    let instBounded : Bounded3D Triangle := {
      bounds := fun tri =>
        let (minPt, maxPt) := tri.boundingBox
        AABB.fromMinMax minPt maxPt
    }
    let instInhabited : Inhabited Triangle := {
      default := Triangle.mk' Vec3.zero Vec3.zero Vec3.zero
    }
    some (@BVH.build Triangle instBounded instInhabited tris config)

/-- Ray cast against a mesh using a BVH. -/
def rayCast (mesh : Mesh) (bvh : BVH) (ray : Ray) (cullBackface : Bool := false) : Option MeshHit :=
  let hitTest := fun idx =>
    match Mesh.triangle? mesh idx with
    | none => none
    | some tri => Intersection.rayTriangle ray tri cullBackface
  match BVH.rayCast bvh ray hitTest with
  | none => none
  | some hit =>
    match Mesh.triangle? mesh hit.index with
    | none => none
    | some tri =>
      let bc := tri.barycentric hit.point
      some ⟨hit.index, hit.t, hit.point, hit.normal, bc.u, bc.v, bc.w⟩

/-- Ray cast against a mesh using a BVH, returning all hits. -/
def rayCastAll (mesh : Mesh) (bvh : BVH) (ray : Ray) (cullBackface : Bool := false) : Array MeshHit :=
  let hitTest := fun idx =>
    match Mesh.triangle? mesh idx with
    | none => none
    | some tri => Intersection.rayTriangle ray tri cullBackface
  let hits := BVH.rayCastAll bvh ray hitTest
  hits.foldl (fun acc h =>
    match Mesh.triangle? mesh h.index with
    | none => acc
    | some tri =>
      let bc := tri.barycentric h.point
      acc.push ⟨h.index, h.t, h.point, h.normal, bc.u, bc.v, bc.w⟩
  ) #[]

/-- Check if any triangle is hit by the ray within maxT using BVH. -/
def rayAny (mesh : Mesh) (bvh : BVH) (ray : Ray) (maxT : Float := Float.infinity)
    (cullBackface : Bool := false) : Bool :=
  let hitTest := fun idx =>
    match Mesh.triangle? mesh idx with
    | none => false
    | some tri =>
      match Intersection.rayTriangle ray tri cullBackface with
      | none => false
      | some hit => hit.t <= maxT
  BVH.rayAny bvh ray maxT hitTest

end MeshBVH

end Linalg.Spatial
