/-
  Bounding Volume Hierarchy for ray tracing acceleration.
  Uses Surface Area Heuristic (SAH) for optimal tree construction.
-/

import Linalg.Spatial.Common
import Linalg.Geometry.Intersection

namespace Linalg.Spatial

open Linalg

/-- BVH node. -/
inductive BVHNode where
  /-- Branch node with two children. -/
  | branch (bounds : AABB) (left right : BVHNode)
  /-- Leaf node containing item indices. -/
  | leaf (bounds : AABB) (indices : Array Nat)
deriving Repr, Inhabited

namespace BVHNode

/-- Get the bounds of a node. -/
def bounds : BVHNode → AABB
  | branch b _ _ => b
  | leaf b _ => b

/-- Check if a node is a leaf. -/
def isLeaf : BVHNode → Bool
  | branch _ _ _ => false
  | leaf _ _ => true

end BVHNode

/-- Bounding Volume Hierarchy. -/
structure BVH where
  /-- Root node. -/
  root : BVHNode
  /-- Total number of items. -/
  itemCount : Nat
  /-- Configuration used for building. -/
  config : BVHConfig
deriving Repr, Inhabited

namespace BVH

/-- Bucket for SAH evaluation. -/
private structure SAHBucket where
  count : Nat
  bounds : AABB
deriving Repr

instance : Inhabited SAHBucket where
  default := { count := 0, bounds := AABB.fromMinMax Vec3.zero Vec3.zero }

/-- Create an empty bucket. -/
private def SAHBucket.empty : SAHBucket :=
  { count := 0, bounds := AABB.fromMinMax Vec3.zero Vec3.zero }

/-- Merge a primitive into a bucket. -/
private def SAHBucket.add (b : SAHBucket) (primBounds : AABB) : SAHBucket :=
  if b.count == 0 then
    { count := 1, bounds := primBounds }
  else
    { count := b.count + 1, bounds := AABB.merge b.bounds primBounds }

/-- Primitive info for building. -/
private structure PrimInfo where
  index : Nat
  bounds : AABB
  centroid : Vec3
deriving Repr

instance : Inhabited PrimInfo where
  default := { index := 0, bounds := default, centroid := Vec3.zero }

private def getAxisValue (v : Vec3) (axis : Fin 3) : Float :=
  match axis with
  | ⟨0, _⟩ => v.x
  | ⟨1, _⟩ => v.y
  | ⟨2, _⟩ => v.z

private def getAxisExtent (b : AABB) (axis : Fin 3) : Float :=
  match axis with
  | ⟨0, _⟩ => b.max.x - b.min.x
  | ⟨1, _⟩ => b.max.y - b.min.y
  | ⟨2, _⟩ => b.max.z - b.min.z

private def partitionPrims (prims : Array PrimInfo) (axis : Fin 3) (splitBucket : Nat)
    (centroidBounds : AABB) (config : BVHConfig) : Array PrimInfo × Array PrimInfo :=
  let nBuckets := config.sahBuckets
  let extent := getAxisExtent centroidBounds axis
  prims.foldl (fun (left, right) prim =>
    let offset := (getAxisValue prim.centroid axis - getAxisValue centroidBounds.min axis) / extent
    let b := min ((offset * nBuckets.toFloat).toUInt32.toNat) (nBuckets - 1)
    if b < splitBucket then (left.push prim, right)
    else (left, right.push prim)
  ) (#[], #[])

private def findBestSplit (prims : Array PrimInfo) (centroidBounds : AABB) (nodeBounds : AABB)
    (config : BVHConfig) : Fin 3 × Nat × Float :=
  let nBuckets := config.sahBuckets
  let axes := #[⟨0, by omega⟩, ⟨1, by omega⟩, ⟨2, by omega⟩]
  -- Find best split across all axes
  axes.foldl (fun (bestAxis, bestSplit, bestCost) axis =>
    let extent := getAxisExtent centroidBounds axis
    if extent < 0.0001 then
      (bestAxis, bestSplit, bestCost)
    else
      -- Initialize buckets
      let buckets := (List.range nBuckets).foldl (fun arr _ =>
        arr.push SAHBucket.empty
      ) #[]
      -- Assign primitives to buckets
      let buckets := prims.foldl (fun bkts prim =>
        let offset := (getAxisValue prim.centroid axis - getAxisValue centroidBounds.min axis) / extent
        let b := min ((offset * nBuckets.toFloat).toUInt32.toNat) (nBuckets - 1)
        bkts.set! b (bkts[b]!.add prim.bounds)
      ) buckets
      -- Evaluate SAH cost for each split position
      let (_, splitCost) := (List.range (nBuckets - 1)).foldl (fun (bestIdx, bestC) i =>
        -- Compute cost for splitting after bucket i
        let (leftCount, leftBounds) := (List.range (i + 1)).foldl (fun (c, b) j =>
          let bkt := buckets[j]!
          if bkt.count == 0 then (c, b)
          else if c == 0 then (bkt.count, bkt.bounds)
          else (c + bkt.count, AABB.merge b bkt.bounds)
        ) (0, AABB.fromMinMax Vec3.zero Vec3.zero)
        let (rightCount, rightBounds) := (List.range (nBuckets - i - 1)).foldl (fun (c, b) j =>
          let bkt := buckets[i + 1 + j]!
          if bkt.count == 0 then (c, b)
          else if c == 0 then (bkt.count, bkt.bounds)
          else (c + bkt.count, AABB.merge b bkt.bounds)
        ) (0, AABB.fromMinMax Vec3.zero Vec3.zero)
        if leftCount == 0 || rightCount == 0 then (bestIdx, bestC)
        else
          let cost := config.traversalCost +
            (leftBounds.surfaceArea / nodeBounds.surfaceArea) * leftCount.toFloat * config.intersectionCost +
            (rightBounds.surfaceArea / nodeBounds.surfaceArea) * rightCount.toFloat * config.intersectionCost
          if cost < bestC then (i + 1, cost) else (bestIdx, bestC)
      ) (0, Float.infinity)
      if splitCost < bestCost then (axis, splitCost.toUInt32.toNat, splitCost)
      else (bestAxis, bestSplit, bestCost)
  ) (⟨0, by omega⟩, 0, Float.infinity)

private partial def buildBVHNode (prims : Array PrimInfo) (config : BVHConfig) : BVHNode :=
  if prims.size <= config.maxLeafItems then
    -- Create leaf
    let bounds := prims.foldl (fun acc p => AABB.merge acc p.bounds)
      prims[0]!.bounds
    BVHNode.leaf bounds (prims.map (·.index))
  else
    -- Compute overall bounds and centroid bounds
    let bounds := prims.foldl (fun acc p => AABB.merge acc p.bounds)
      prims[0]!.bounds
    let centroidBounds := prims.foldl (fun acc p => AABB.expand acc p.centroid)
      (AABB.fromPoint prims[0]!.centroid)

    -- Find best split using SAH
    let (bestAxis, bestSplit, bestCost) := findBestSplit prims centroidBounds bounds config

    -- Compare with leaf cost
    let leafCost := prims.size.toFloat * config.intersectionCost
    if bestCost >= leafCost || bestSplit == 0 || bestSplit >= prims.size then
      -- Leaf is cheaper
      BVHNode.leaf bounds (prims.map (·.index))
    else
      -- Partition primitives
      let (leftPrims, rightPrims) := partitionPrims prims bestAxis bestSplit centroidBounds config
      if leftPrims.isEmpty || rightPrims.isEmpty then
        -- Couldn't partition, make leaf
        BVHNode.leaf bounds (prims.map (·.index))
      else
        -- Recursively build children
        let left := buildBVHNode leftPrims config
        let right := buildBVHNode rightPrims config
        BVHNode.branch bounds left right

/-- Build a BVH from an array of items using SAH. -/
def build {α : Type} [Bounded3D α] [Inhabited α] (items : Array α) (config : BVHConfig := {}) : BVH :=
  if items.isEmpty then
    { root := BVHNode.leaf (AABB.fromMinMax Vec3.zero Vec3.one) #[]
      itemCount := 0
      config }
  else
    -- Create primitive info
    let primInfos := items.mapIdx fun i item =>
      let b := Bounded3D.bounds item
      { index := i, bounds := b, centroid := b.center : PrimInfo }
    -- Build tree
    let root := buildBVHNode primInfos config
    { root, itemCount := items.size, config }

/-- Check a BVH node for ray intersection. -/
private partial def checkBVHNode (node : BVHNode) (ray : Ray) (hitTest : Nat → Option RayHit)
    (closest : Option SpatialHit) : Option SpatialHit :=
  match node with
  | BVHNode.leaf _ indices =>
    -- Test all primitives in leaf
    indices.foldl (fun acc idx =>
      match hitTest idx with
      | none => acc
      | some hit =>
        let maxT := match acc with
          | some a => a.t
          | none => Float.infinity
        if hit.t < maxT then
          some { index := idx, t := hit.t, point := hit.point, normal := hit.normal }
        else acc
    ) closest
  | BVHNode.branch _ left right =>
    -- Get intersection intervals for both children
    let leftInt := Intersection.rayAABB ray left.bounds
    let rightInt := Intersection.rayAABB ray right.bounds
    match leftInt, rightInt with
    | none, none => closest
    | some _, none => rayCastBVHNode left ray hitTest closest
    | none, some _ => rayCastBVHNode right ray hitTest closest
    | some (lMin, _), some (rMin, _) =>
      -- Visit closer child first
      let (first, second) := if lMin < rMin then (left, right) else (right, left)
      let closest' := rayCastBVHNode first ray hitTest closest
      rayCastBVHNode second ray hitTest closest'
where
  rayCastBVHNode (node : BVHNode) (ray : Ray) (hitTest : Nat → Option RayHit)
      (closest : Option SpatialHit) : Option SpatialHit :=
    match Intersection.rayAABB ray node.bounds with
    | none => closest
    | some (tMin, _) =>
      match closest with
      | some hit => if tMin > hit.t then closest else checkBVHNode node ray hitTest closest
      | none => checkBVHNode node ray hitTest closest

/-- Cast a ray through the BVH, finding the closest hit. -/
def rayCast (bvh : BVH) (ray : Ray) (hitTest : Nat → Option RayHit) : Option SpatialHit :=
  match Intersection.rayAABB ray bvh.root.bounds with
  | none => none
  | some _ => checkBVHNode bvh.root ray hitTest none

/-- Cast all rays through BVH node. -/
private partial def rayCastAllBVHNode (node : BVHNode) (ray : Ray) (hitTest : Nat → Option RayHit)
    (acc : Array SpatialHit) : Array SpatialHit :=
  match Intersection.rayAABB ray node.bounds with
  | none => acc
  | some _ =>
    match node with
    | BVHNode.leaf _ indices =>
      indices.foldl (fun a idx =>
        match hitTest idx with
        | none => a
        | some hit => a.push { index := idx, t := hit.t, point := hit.point, normal := hit.normal }
      ) acc
    | BVHNode.branch _ left right =>
      let acc' := rayCastAllBVHNode left ray hitTest acc
      rayCastAllBVHNode right ray hitTest acc'

/-- Cast a ray through the BVH, returning all hits. -/
def rayCastAll (bvh : BVH) (ray : Ray) (hitTest : Nat → Option RayHit) : Array SpatialHit :=
  rayCastAllBVHNode bvh.root ray hitTest #[]

/-- Check if any ray hits in node. -/
private partial def rayAnyBVHNode (node : BVHNode) (ray : Ray) (maxT : Float) (hitTest : Nat → Bool) : Bool :=
  match Intersection.rayAABB ray node.bounds with
  | none => false
  | some (tMin, _) =>
    if tMin > maxT then false
    else match node with
      | BVHNode.leaf _ indices => indices.any hitTest
      | BVHNode.branch _ left right =>
        rayAnyBVHNode left ray maxT hitTest || rayAnyBVHNode right ray maxT hitTest

/-- Check if ray hits any object (shadow ray optimization). -/
def rayAny (bvh : BVH) (ray : Ray) (maxT : Float) (hitTest : Nat → Bool) : Bool :=
  rayAnyBVHNode bvh.root ray maxT hitTest

/-- Query AABB in BVH node. -/
private partial def queryAABBBVHNode (node : BVHNode) (query : AABB) (acc : Array Nat) : Array Nat :=
  if !Intersection.aabbAABB node.bounds query then acc
  else match node with
    | BVHNode.leaf _ indices => acc ++ indices
    | BVHNode.branch _ left right =>
      let acc' := queryAABBBVHNode left query acc
      queryAABBBVHNode right query acc'

/-- Query all items whose bounds intersect the given AABB. -/
def queryAABB (bvh : BVH) (query : AABB) : Array Nat :=
  queryAABBBVHNode bvh.root query #[]

/-- Query all items whose bounds intersect the given sphere. -/
def querySphere (bvh : BVH) (center : Vec3) (radius : Float) : Array Nat :=
  let sphereAABB := AABB.fromCenterExtents center (Vec3.mk radius radius radius)
  queryAABB bvh sphereAABB

/-- Refit BVH node. -/
private partial def refitBVHNode {α : Type} [Bounded3D α] (node : BVHNode) (items : Array α) : BVHNode :=
  match node with
  | BVHNode.leaf _ indices =>
    if indices.isEmpty then node
    else
      let bounds := indices.foldl (fun acc idx =>
        if h : idx < items.size then
          AABB.merge acc (Bounded3D.bounds items[idx])
        else acc
      ) (if h : indices[0]! < items.size
         then Bounded3D.bounds items[indices[0]!]
         else AABB.fromMinMax Vec3.zero Vec3.one)
      BVHNode.leaf bounds indices
  | BVHNode.branch _ left right =>
    let left' := refitBVHNode left items
    let right' := refitBVHNode right items
    let bounds := AABB.merge left'.bounds right'.bounds
    BVHNode.branch bounds left' right'

/-- Refit the BVH after items have moved (updates bounds only). -/
def refit {α : Type} [Bounded3D α] (bvh : BVH) (items : Array α) : BVH :=
  { bvh with root := refitBVHNode bvh.root items }

/-- Compute BVH node depth. -/
private partial def bvhNodeDepth : BVHNode → Nat
  | BVHNode.leaf _ _ => 0
  | BVHNode.branch _ left right => 1 + Nat.max (bvhNodeDepth left) (bvhNodeDepth right)

/-- Get tree depth. -/
def depth (bvh : BVH) : Nat :=
  bvhNodeDepth bvh.root

/-- Count BVH nodes. -/
private partial def countBVHNodes : BVHNode → Nat
  | BVHNode.leaf _ _ => 1
  | BVHNode.branch _ left right => 1 + countBVHNodes left + countBVHNodes right

/-- Count total nodes. -/
def nodeCount (bvh : BVH) : Nat :=
  countBVHNodes bvh.root

/-- Count BVH leaf nodes. -/
private partial def countBVHLeaves : BVHNode → Nat
  | BVHNode.leaf _ _ => 1
  | BVHNode.branch _ left right => countBVHLeaves left + countBVHLeaves right

/-- Count leaf nodes. -/
def leafCount (bvh : BVH) : Nat :=
  countBVHLeaves bvh.root

/-- Check if the BVH is empty. -/
def isEmpty (bvh : BVH) : Bool := bvh.itemCount == 0

end BVH

end Linalg.Spatial
