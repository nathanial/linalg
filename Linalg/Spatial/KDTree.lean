/-
  KD-Tree for efficient point queries.
  Supports k-nearest neighbor and range queries.
-/

import Linalg.Spatial.Common

namespace Linalg.Spatial

/-- KD-Tree node for 3D points. -/
inductive KDNode3D where
  /-- Branch node splitting on an axis. -/
  | branch (splitAxis : Fin 3) (splitValue : Float) (left right : KDNode3D)
  /-- Leaf node containing point indices. -/
  | leaf (indices : Array Nat)
deriving Repr, Inhabited

/-- KD-Tree node for 2D points. -/
inductive KDNode2D where
  /-- Branch node splitting on an axis. -/
  | branch (splitAxis : Fin 2) (splitValue : Float) (left right : KDNode2D)
  /-- Leaf node containing point indices. -/
  | leaf (indices : Array Nat)
deriving Repr, Inhabited

/-- 3D KD-Tree. -/
structure KDTree3D where
  /-- Root node. -/
  root : KDNode3D
  /-- Bounds of all points. -/
  bounds : AABB
  /-- Number of points. -/
  pointCount : Nat
deriving Repr, Inhabited

/-- 2D KD-Tree. -/
structure KDTree2D where
  /-- Root node. -/
  root : KDNode2D
  /-- Bounds of all points. -/
  bounds : AABB2D
  /-- Number of points. -/
  pointCount : Nat
deriving Repr, Inhabited

namespace KDTree3D

/-- Get component of Vec3 by axis index. -/
private def getAxis (v : Vec3) (axis : Fin 3) : Float :=
  match axis with
  | ⟨0, _⟩ => v.x
  | ⟨1, _⟩ => v.y
  | ⟨2, _⟩ => v.z

/-- Build KD-Tree node recursively. -/
private partial def buildNode3D (points : Array (Nat × Vec3)) (depth : Nat) (maxLeaf : Nat) : KDNode3D :=
  if points.size <= maxLeaf then
    KDNode3D.leaf (points.map (·.1))
  else
    -- Choose axis based on depth (cycling)
    let axis : Fin 3 := ⟨depth % 3, by omega⟩
    -- Sort by axis
    let sorted := points.qsort fun (_, a) (_, b) => getAxis a axis < getAxis b axis
    -- Split at median
    let mid := sorted.size / 2
    let splitValue := getAxis sorted[mid]!.2 axis
    let left := sorted.extract 0 mid
    let right := sorted.extract mid sorted.size
    if left.isEmpty || right.isEmpty then
      -- Can't split further, make leaf
      KDNode3D.leaf (points.map (·.1))
    else
      KDNode3D.branch axis splitValue
        (buildNode3D left (depth + 1) maxLeaf)
        (buildNode3D right (depth + 1) maxLeaf)

/-- Build a KD-Tree from an array of items with positions. -/
def build {α : Type} [HasPosition3D α] [Inhabited α] (items : Array α) (maxLeaf : Nat := 8) : KDTree3D :=
  if items.isEmpty then
    { root := KDNode3D.leaf #[]
      bounds := AABB.fromMinMax Vec3.zero Vec3.one
      pointCount := 0 }
  else
    -- Create array of (index, position) pairs
    let indexed := items.mapIdx fun i item => (i, HasPosition3D.position item)
    -- Compute bounds
    let bounds := indexed.foldl (fun acc (_, p) => AABB.expand acc p)
      (AABB.fromPoint indexed[0]!.2)
    -- Build tree
    let root := buildNode3D indexed 0 maxLeaf
    { root, bounds, pointCount := items.size }

/-- Find nearest neighbor in node. -/
private partial def nearestNode3D {α : Type} [HasPosition3D α] (node : KDNode3D) (items : Array α)
    (query : Vec3) (best : Option Nat × Float) : Option Nat × Float :=
  match node with
  | KDNode3D.leaf indices =>
    indices.foldl (fun (bestIdx, bestDist) idx =>
      if h : idx < items.size then
        let pos := HasPosition3D.position items[idx]
        let dist := query.distanceSquared pos
        if dist < bestDist then (some idx, dist) else (bestIdx, bestDist)
      else (bestIdx, bestDist)
    ) best
  | KDNode3D.branch axis splitValue left right =>
    let queryValue := getAxis query axis
    let (nearChild, farChild) :=
      if queryValue < splitValue then (left, right) else (right, left)
    -- Search near child first
    let (bestIdx, bestDist) := nearestNode3D nearChild items query best
    -- Check if far child could have closer points
    let distToPlane := (queryValue - splitValue) * (queryValue - splitValue)
    if distToPlane < bestDist then
      nearestNode3D farChild items query (bestIdx, bestDist)
    else
      (bestIdx, bestDist)

/-- Find the nearest neighbor to a query point. -/
def nearest {α : Type} [HasPosition3D α] (tree : KDTree3D) (items : Array α) (query : Vec3) : Option Nat :=
  if tree.pointCount == 0 then none
  else
    let (idx, _) := nearestNode3D tree.root items query (none, Float.infinity)
    idx

/-- K-nearest search in node. -/
private partial def kNearestNode3D {α : Type} [HasPosition3D α] (node : KDNode3D) (items : Array α)
    (query : Vec3) (k : Nat) (heap : MaxHeap Nat) : MaxHeap Nat :=
  match node with
  | KDNode3D.leaf indices =>
    indices.foldl (fun h idx =>
      if h_idx : idx < items.size then
        let pos := HasPosition3D.position items[idx]
        let dist := query.distanceSquared pos
        h.insert dist idx
      else h
    ) heap
  | KDNode3D.branch axis splitValue left right =>
    let queryValue := getAxis query axis
    let (nearChild, farChild) :=
      if queryValue < splitValue then (left, right) else (right, left)
    -- Search near child first
    let heap' := kNearestNode3D nearChild items query k heap
    -- Check if far child could have closer points
    let distToPlane := (queryValue - splitValue) * (queryValue - splitValue)
    if !heap'.isFull || distToPlane < heap'.maxDistance then
      kNearestNode3D farChild items query k heap'
    else
      heap'

/-- Find k nearest neighbors to a query point. -/
def kNearest {α : Type} [HasPosition3D α] (tree : KDTree3D) (items : Array α) (query : Vec3) (k : Nat) : Array Nat :=
  if k == 0 || tree.pointCount == 0 then #[]
  else
    let heap := kNearestNode3D tree.root items query k (MaxHeap.empty k)
    heap.toSortedArray

/-- Radius search in node. -/
private partial def withinRadiusNode3D {α : Type} [HasPosition3D α] (node : KDNode3D) (items : Array α)
    (query : Vec3) (radiusSq : Float) (acc : Array Nat) : Array Nat :=
  match node with
  | KDNode3D.leaf indices =>
    indices.foldl (fun a idx =>
      if h : idx < items.size then
        let pos := HasPosition3D.position items[idx]
        let dist := query.distanceSquared pos
        if dist <= radiusSq then a.push idx else a
      else a
    ) acc
  | KDNode3D.branch axis splitValue left right =>
    let queryValue := getAxis query axis
    -- Check which children could contain points within radius
    let acc' := if queryValue - Float.sqrt radiusSq <= splitValue then
      withinRadiusNode3D left items query radiusSq acc
    else acc
    if queryValue + Float.sqrt radiusSq >= splitValue then
      withinRadiusNode3D right items query radiusSq acc'
    else acc'

/-- Find all points within a radius of the query point. -/
def withinRadius {α : Type} [HasPosition3D α] (tree : KDTree3D) (items : Array α) (query : Vec3) (radius : Float) : Array Nat :=
  withinRadiusNode3D tree.root items query (radius * radius) #[]

/-- AABB search in node. -/
private partial def withinAABBNode3D {α : Type} [HasPosition3D α] (node : KDNode3D) (items : Array α)
    (query : AABB) (acc : Array Nat) : Array Nat :=
  match node with
  | KDNode3D.leaf indices =>
    indices.foldl (fun a idx =>
      if h : idx < items.size then
        let pos := HasPosition3D.position items[idx]
        if query.containsPoint pos then a.push idx else a
      else a
    ) acc
  | KDNode3D.branch axis splitValue left right =>
    let minValue := getAxis query.min axis
    let maxValue := getAxis query.max axis
    let acc' := if minValue <= splitValue then
      withinAABBNode3D left items query acc
    else acc
    if maxValue >= splitValue then
      withinAABBNode3D right items query acc'
    else acc'

/-- Find all points within an AABB. -/
def withinAABB {α : Type} [HasPosition3D α] (tree : KDTree3D) (items : Array α) (query : AABB) : Array Nat :=
  withinAABBNode3D tree.root items query #[]

/-- Compute node depth. -/
private partial def nodeDepth3D : KDNode3D → Nat
  | KDNode3D.leaf _ => 0
  | KDNode3D.branch _ _ left right => 1 + Nat.max (nodeDepth3D left) (nodeDepth3D right)

/-- Get tree depth. -/
def depth (tree : KDTree3D) : Nat :=
  nodeDepth3D tree.root

/-- Count nodes. -/
private partial def countNodes3D : KDNode3D → Nat
  | KDNode3D.leaf _ => 1
  | KDNode3D.branch _ _ left right => 1 + countNodes3D left + countNodes3D right

/-- Count total nodes. -/
def nodeCount (tree : KDTree3D) : Nat :=
  countNodes3D tree.root

end KDTree3D

namespace KDTree2D

/-- Get component of Vec2 by axis index. -/
private def getAxis (v : Vec2) (axis : Fin 2) : Float :=
  match axis with
  | ⟨0, _⟩ => v.x
  | ⟨1, _⟩ => v.y

/-- Build KD-Tree node recursively. -/
private partial def buildNode2D (points : Array (Nat × Vec2)) (depth : Nat) (maxLeaf : Nat) : KDNode2D :=
  if points.size <= maxLeaf then
    KDNode2D.leaf (points.map (·.1))
  else
    let axis : Fin 2 := ⟨depth % 2, by omega⟩
    let sorted := points.qsort fun (_, a) (_, b) => getAxis a axis < getAxis b axis
    let mid := sorted.size / 2
    let splitValue := getAxis sorted[mid]!.2 axis
    let left := sorted.extract 0 mid
    let right := sorted.extract mid sorted.size
    if left.isEmpty || right.isEmpty then
      KDNode2D.leaf (points.map (·.1))
    else
      KDNode2D.branch axis splitValue
        (buildNode2D left (depth + 1) maxLeaf)
        (buildNode2D right (depth + 1) maxLeaf)

/-- Build a KD-Tree from an array of items with positions. -/
def build {α : Type} [HasPosition2D α] [Inhabited α] (items : Array α) (maxLeaf : Nat := 8) : KDTree2D :=
  if items.isEmpty then
    { root := KDNode2D.leaf #[]
      bounds := AABB2D.fromMinMax Vec2.zero Vec2.one
      pointCount := 0 }
  else
    let indexed := items.mapIdx fun i item => (i, HasPosition2D.position item)
    let bounds := indexed.foldl (fun acc (_, p) => AABB2D.expand acc p)
      (AABB2D.fromPoint indexed[0]!.2)
    let root := buildNode2D indexed 0 maxLeaf
    { root, bounds, pointCount := items.size }

/-- Find nearest neighbor in node. -/
private partial def nearestNode2D {α : Type} [HasPosition2D α] (node : KDNode2D) (items : Array α)
    (query : Vec2) (best : Option Nat × Float) : Option Nat × Float :=
  match node with
  | KDNode2D.leaf indices =>
    indices.foldl (fun (bestIdx, bestDist) idx =>
      if h : idx < items.size then
        let pos := HasPosition2D.position items[idx]
        let dist := query.distanceSquared pos
        if dist < bestDist then (some idx, dist) else (bestIdx, bestDist)
      else (bestIdx, bestDist)
    ) best
  | KDNode2D.branch axis splitValue left right =>
    let queryValue := getAxis query axis
    let (nearChild, farChild) :=
      if queryValue < splitValue then (left, right) else (right, left)
    let (bestIdx, bestDist) := nearestNode2D nearChild items query best
    let distToPlane := (queryValue - splitValue) * (queryValue - splitValue)
    if distToPlane < bestDist then
      nearestNode2D farChild items query (bestIdx, bestDist)
    else
      (bestIdx, bestDist)

/-- Find the nearest neighbor to a query point. -/
def nearest {α : Type} [HasPosition2D α] (tree : KDTree2D) (items : Array α) (query : Vec2) : Option Nat :=
  if tree.pointCount == 0 then none
  else
    let (idx, _) := nearestNode2D tree.root items query (none, Float.infinity)
    idx

/-- K-nearest search in node. -/
private partial def kNearestNode2D {α : Type} [HasPosition2D α] (node : KDNode2D) (items : Array α)
    (query : Vec2) (k : Nat) (heap : MaxHeap Nat) : MaxHeap Nat :=
  match node with
  | KDNode2D.leaf indices =>
    indices.foldl (fun h idx =>
      if h_idx : idx < items.size then
        let pos := HasPosition2D.position items[idx]
        let dist := query.distanceSquared pos
        h.insert dist idx
      else h
    ) heap
  | KDNode2D.branch axis splitValue left right =>
    let queryValue := getAxis query axis
    let (nearChild, farChild) :=
      if queryValue < splitValue then (left, right) else (right, left)
    let heap' := kNearestNode2D nearChild items query k heap
    let distToPlane := (queryValue - splitValue) * (queryValue - splitValue)
    if !heap'.isFull || distToPlane < heap'.maxDistance then
      kNearestNode2D farChild items query k heap'
    else
      heap'

/-- Find k nearest neighbors to a query point. -/
def kNearest {α : Type} [HasPosition2D α] (tree : KDTree2D) (items : Array α) (query : Vec2) (k : Nat) : Array Nat :=
  if k == 0 || tree.pointCount == 0 then #[]
  else
    let heap := kNearestNode2D tree.root items query k (MaxHeap.empty k)
    heap.toSortedArray

/-- Radius search in node. -/
private partial def withinRadiusNode2D {α : Type} [HasPosition2D α] (node : KDNode2D) (items : Array α)
    (query : Vec2) (radiusSq : Float) (acc : Array Nat) : Array Nat :=
  match node with
  | KDNode2D.leaf indices =>
    indices.foldl (fun a idx =>
      if h : idx < items.size then
        let pos := HasPosition2D.position items[idx]
        let dist := query.distanceSquared pos
        if dist <= radiusSq then a.push idx else a
      else a
    ) acc
  | KDNode2D.branch axis splitValue left right =>
    let queryValue := getAxis query axis
    let acc' := if queryValue - Float.sqrt radiusSq <= splitValue then
      withinRadiusNode2D left items query radiusSq acc
    else acc
    if queryValue + Float.sqrt radiusSq >= splitValue then
      withinRadiusNode2D right items query radiusSq acc'
    else acc'

/-- Find all points within a radius of the query point. -/
def withinRadius {α : Type} [HasPosition2D α] (tree : KDTree2D) (items : Array α) (query : Vec2) (radius : Float) : Array Nat :=
  withinRadiusNode2D tree.root items query (radius * radius) #[]

/-- Rectangle search in node. -/
private partial def withinRectNode2D {α : Type} [HasPosition2D α] (node : KDNode2D) (items : Array α)
    (query : AABB2D) (acc : Array Nat) : Array Nat :=
  match node with
  | KDNode2D.leaf indices =>
    indices.foldl (fun a idx =>
      if h : idx < items.size then
        let pos := HasPosition2D.position items[idx]
        if query.containsPoint pos then a.push idx else a
      else a
    ) acc
  | KDNode2D.branch axis splitValue left right =>
    let minValue := getAxis query.min axis
    let maxValue := getAxis query.max axis
    let acc' := if minValue <= splitValue then
      withinRectNode2D left items query acc
    else acc
    if maxValue >= splitValue then
      withinRectNode2D right items query acc'
    else acc'

/-- Find all points within a rectangle. -/
def withinRect {α : Type} [HasPosition2D α] (tree : KDTree2D) (items : Array α) (query : AABB2D) : Array Nat :=
  withinRectNode2D tree.root items query #[]

/-- Compute node depth. -/
private partial def nodeDepth2D : KDNode2D → Nat
  | KDNode2D.leaf _ => 0
  | KDNode2D.branch _ _ left right => 1 + Nat.max (nodeDepth2D left) (nodeDepth2D right)

/-- Get tree depth. -/
def depth (tree : KDTree2D) : Nat :=
  nodeDepth2D tree.root

/-- Count nodes. -/
private partial def countNodes2D : KDNode2D → Nat
  | KDNode2D.leaf _ => 1
  | KDNode2D.branch _ _ left right => 1 + countNodes2D left + countNodes2D right

/-- Count total nodes. -/
def nodeCount (tree : KDTree2D) : Nat :=
  countNodes2D tree.root

end KDTree2D

end Linalg.Spatial
