/-
  Octree for 3D spatial partitioning.
  Recursively subdivides space into 8 octants.
-/

import Linalg.Spatial.Common
import Linalg.Geometry.Frustum
import Linalg.Geometry.Intersection

namespace Linalg.Spatial

/-- Octant index (0-7). Bits: X=1, Y=2, Z=4 -/
abbrev OctantIndex := Fin 8

/-- Get octant for a point relative to center.
    Bit 0 = X >= center.x, Bit 1 = Y >= center.y, Bit 2 = Z >= center.z -/
def octantFor (center point : Vec3) : OctantIndex :=
  let xBit : Fin 2 := if point.x >= center.x then 1 else 0
  let yBit : Fin 2 := if point.y >= center.y then 1 else 0
  let zBit : Fin 2 := if point.z >= center.z then 1 else 0
  ⟨xBit.val + yBit.val * 2 + zBit.val * 4, by
    have hx := xBit.isLt; have hy := yBit.isLt; have hz := zBit.isLt; omega⟩

/-- Helper to create a valid octant index from a Nat known to be < 8. -/
private def mkOctant (n : Nat) : OctantIndex :=
  ⟨n % 8, Nat.mod_lt n (by omega)⟩

/-- Stored item reference inside an octree leaf. -/
structure OctreeItem where
  idx : Nat
  bounds : AABB
deriving Repr, Inhabited

/-- Octree node. -/
inductive OctreeNode where
  /-- Internal node with 8 children (some may be empty). -/
  | internal (bounds : AABB) (children : Array (Option OctreeNode))
  /-- Leaf node containing item references. -/
  | leaf (bounds : AABB) (items : Array OctreeItem)
deriving Repr, Inhabited

namespace OctreeNode

/-- Get the bounds of a node. -/
def bounds : OctreeNode → AABB
  | internal b _ => b
  | leaf b _ => b

/-- Check if a node is a leaf. -/
def isLeaf : OctreeNode → Bool
  | internal _ _ => false
  | leaf _ _ => true

/-- Get child bounds for a given octant. -/
def childBounds (b : AABB) (o : OctantIndex) : AABB :=
  let c := b.center
  let minX := if o.val &&& 1 != 0 then c.x else b.min.x
  let maxX := if o.val &&& 1 != 0 then b.max.x else c.x
  let minY := if o.val &&& 2 != 0 then c.y else b.min.y
  let maxY := if o.val &&& 2 != 0 then b.max.y else c.y
  let minZ := if o.val &&& 4 != 0 then c.z else b.min.z
  let maxZ := if o.val &&& 4 != 0 then b.max.z else c.z
  AABB.fromMinMax (Vec3.mk minX minY minZ) (Vec3.mk maxX maxY maxZ)

end OctreeNode

/-- Octree for 3D spatial partitioning. -/
structure Octree where
  /-- Root node of the tree. -/
  root : OctreeNode
  /-- Overall bounds. -/
  bounds : AABB
  /-- Total number of items. -/
  itemCount : Nat
  /-- Configuration. -/
  config : TreeConfig
deriving Repr, Inhabited

namespace Octree

/-- Create an empty octree with given bounds. -/
def empty (bounds : AABB) (config : TreeConfig := {}) : Octree :=
  { root := OctreeNode.leaf bounds #[]
    bounds
    itemCount := 0
    config }

/-- Insert an item into a node, potentially subdividing. -/
private partial def insertNode {α : Type} [Bounded3D α] (node : OctreeNode) (idx : Nat) (item : α)
    (depth : Nat) (config : TreeConfig) : OctreeNode :=
  let itemBounds := Bounded3D.bounds item
  let itemRef : OctreeItem := { idx, bounds := itemBounds }
  insertNodeRef node itemRef depth config
where
  insertNodeRef (node : OctreeNode) (item : OctreeItem) (depth : Nat) (config : TreeConfig) : OctreeNode :=
    match node with
    | OctreeNode.leaf bounds items =>
      let newItems := items.push item
      let size := bounds.size
      let minDim := Float.min size.x (Float.min size.y size.z)
      -- Check if we should split
      if newItems.size > config.maxLeafItems && depth < config.maxDepth &&
         minDim > config.minNodeSize then
        -- Create internal node
        let children := #[none, none, none, none, none, none, none, none]
        let internalNode := OctreeNode.internal bounds children
        -- Re-insert all items into the new internal node
        newItems.foldl (fun n it => insertNodeRef n it depth config) internalNode
      else
        OctreeNode.leaf bounds newItems
    | OctreeNode.internal bounds children =>
      -- Find which octants the item overlaps
      let overlaps := #[0, 1, 2, 3, 4, 5, 6, 7].filter fun o =>
        let childB := OctreeNode.childBounds bounds (mkOctant o)
        Intersection.aabbAABB childB item.bounds
      -- Insert into each overlapping octant
      let newChildren := overlaps.foldl (fun cs o =>
        let childBounds := OctreeNode.childBounds bounds (mkOctant o)
        let existingChild := cs[o]!
        let child := match existingChild with
          | some c => insertNodeRef c item (depth + 1) config
          | none => insertNodeRef (OctreeNode.leaf childBounds #[]) item (depth + 1) config
        cs.set! o (some child)
      ) children
      OctreeNode.internal bounds newChildren

/-- Insert a single item into the octree. -/
def insert {α : Type} [Bounded3D α] (tree : Octree) (idx : Nat) (item : α) : Octree :=
  { tree with
    root := insertNode tree.root idx item 0 tree.config
    itemCount := tree.itemCount + 1 }

/-- Build an octree from an array of items. -/
def build {α : Type} [Bounded3D α] [Inhabited α] (items : Array α) (config : TreeConfig := {}) : Octree :=
  if items.isEmpty then
    empty (AABB.fromMinMax Vec3.zero Vec3.one) config
  else
    -- Compute overall bounds with small padding
    let bounds := items.foldl (fun acc item =>
      AABB.merge acc (Bounded3D.bounds item)
    ) (Bounded3D.bounds items[0]!)
    -- Add small padding
    let padding := bounds.size.scale 0.01
    let paddedBounds := AABB.fromMinMax (bounds.min.sub padding) (bounds.max.add padding)
    -- Insert all items
    items.foldl (fun (tree, idx) item =>
      (tree.insert idx item, idx + 1)
    ) (empty paddedBounds config, 0) |>.1

/-- Query a node for items in an AABB. -/
private partial def queryNodeAABB (node : OctreeNode) (query : AABB) : Array Nat :=
  if !Intersection.aabbAABB node.bounds query then #[]
  else match node with
    | OctreeNode.leaf _ items => items.map (·.idx)
    | OctreeNode.internal _ children =>
      children.foldl (fun acc child =>
        match child with
        | some c => acc ++ queryNodeAABB c query
        | none => acc
      ) #[]

/-- Query candidate items that may intersect the given AABB.
    Returns indices of items in leaf nodes that overlap the query region.
    Note: This is a spatial acceleration query - candidates may not actually
    intersect the AABB. Use `queryAABBExact` for precise intersection testing. -/
def queryAABB (tree : Octree) (query : AABB) : Array Nat :=
  let results := queryNodeAABB tree.root query
  -- Remove duplicates
  results.foldl (fun acc idx =>
    if acc.contains idx then acc else acc.push idx
  ) #[]

/-- Query items whose bounds actually intersect the given AABB.
    Unlike `queryAABB`, this verifies each candidate against its actual bounds. -/
def queryAABBExact {α : Type} [Bounded3D α] (tree : Octree) (items : Array α) (query : AABB) : Array Nat :=
  let candidates := tree.queryAABB query
  candidates.filter fun idx =>
    if h : idx < items.size then
      Intersection.aabbAABB (Bounded3D.bounds items[idx]) query
    else false

/-- Query a node for items in a sphere. -/
private partial def queryNodeSphere (node : OctreeNode) (center : Vec3) (radiusSq : Float) : Array Nat :=
  -- Quick reject: check if sphere's bounding box intersects node
  let radius := Float.sqrt radiusSq
  let sphereBounds := AABB.fromCenterExtents center (Vec3.mk radius radius radius)
  if !Intersection.aabbAABB node.bounds sphereBounds then #[]
  else match node with
    | OctreeNode.leaf _ items => items.map (·.idx)
    | OctreeNode.internal _ children =>
      children.foldl (fun acc child =>
        match child with
        | some c => acc ++ queryNodeSphere c center radiusSq
        | none => acc
      ) #[]

/-- Query candidate items that may intersect the given sphere.
    Returns indices of items in leaf nodes that overlap the sphere's bounding box.
    Note: This is a spatial acceleration query - use for broad-phase collision detection. -/
def querySphere (tree : Octree) (center : Vec3) (radius : Float) : Array Nat :=
  let results := queryNodeSphere tree.root center (radius * radius)
  results.foldl (fun acc idx =>
    if acc.contains idx then acc else acc.push idx
  ) #[]

/-- Query a node for items containing a point. -/
private partial def queryNodePoint (node : OctreeNode) (point : Vec3) : Array Nat :=
  if !node.bounds.containsPoint point then #[]
  else match node with
    | OctreeNode.leaf _ items => items.map (·.idx)
    | OctreeNode.internal bounds children =>
      let o := octantFor bounds.center point
      match children[o.val]! with
      | some child => queryNodePoint child point
      | none => #[]

/-- Query candidate items that may contain the given point.
    Returns indices of items stored in the leaf node containing the point.
    Note: This is a spatial acceleration query - candidates' actual bounds may not
    contain the point. Use `queryPointExact` for precise containment testing. -/
def queryPoint (tree : Octree) (point : Vec3) : Array Nat :=
  queryNodePoint tree.root point

/-- Query items whose bounds actually contain the given point.
    Unlike `queryPoint`, this verifies each candidate against its actual bounds. -/
def queryPointExact {α : Type} [Bounded3D α] (tree : Octree) (items : Array α) (point : Vec3) : Array Nat :=
  let candidates := tree.queryPoint point
  candidates.filter fun idx =>
    if h : idx < items.size then
      (Bounded3D.bounds items[idx]).containsPoint point
    else false

/-- Query a node for items in a frustum. -/
private partial def queryNodeFrustum (node : OctreeNode) (frustum : Frustum) : Array Nat :=
  -- Quick reject using AABB test
  if !frustum.isAABBVisible node.bounds then #[]
  else match node with
    | OctreeNode.leaf _ items => items.map (·.idx)
    | OctreeNode.internal _ children =>
      children.foldl (fun acc child =>
        match child with
        | some c => acc ++ queryNodeFrustum c frustum
        | none => acc
      ) #[]

/-- Query candidate items that may be visible in the frustum.
    Returns indices of items in leaf nodes that overlap the frustum.
    Note: This is a spatial acceleration query - candidates' actual bounds may not
    be fully inside the frustum. Use `queryFrustumExact` for precise visibility testing. -/
def queryFrustum (tree : Octree) (frustum : Frustum) : Array Nat :=
  let results := queryNodeFrustum tree.root frustum
  results.foldl (fun acc idx =>
    if acc.contains idx then acc else acc.push idx
  ) #[]

/-- Query items whose bounds are actually visible in the frustum.
    Unlike `queryFrustum`, this verifies each candidate against its actual bounds. -/
def queryFrustumExact {α : Type} [Bounded3D α] (tree : Octree) (items : Array α) (frustum : Frustum) : Array Nat :=
  let candidates := tree.queryFrustum frustum
  candidates.filter fun idx =>
    if h : idx < items.size then
      frustum.isAABBVisible (Bounded3D.bounds items[idx])
    else false

/-- Ray cast through a node, finding all intersecting items. -/
private partial def raycastNode (node : OctreeNode) (ray : Ray) : Array Nat :=
  -- Check if ray intersects node bounds
  match Intersection.rayAABB ray node.bounds with
  | none => #[]
  | some _ =>
    match node with
    | OctreeNode.leaf _ items => items.map (·.idx)
    | OctreeNode.internal _ children =>
      children.foldl (fun acc child =>
        match child with
        | some c => acc ++ raycastNode c ray
        | none => acc
      ) #[]

/-- Cast a ray through the octree, returning all potentially hit items. -/
def rayCast (tree : Octree) (ray : Ray) : Array Nat :=
  let results := raycastNode tree.root ray
  results.foldl (fun acc idx =>
    if acc.contains idx then acc else acc.push idx
  ) #[]

/-- Count nodes recursively. -/
private partial def countOctreeNodes : OctreeNode → Nat
  | OctreeNode.leaf _ _ => 1
  | OctreeNode.internal _ children =>
    1 + children.foldl (fun acc child =>
      match child with
      | some c => acc + countOctreeNodes c
      | none => acc
    ) 0

/-- Count total nodes in the tree. -/
def nodeCount (tree : Octree) : Nat :=
  countOctreeNodes tree.root

/-- Compute depth recursively. -/
private partial def octreeDepth : OctreeNode → Nat
  | OctreeNode.leaf _ _ => 0
  | OctreeNode.internal _ children =>
    1 + children.foldl (fun acc child =>
      match child with
      | some c => Nat.max acc (octreeDepth c)
      | none => acc
    ) 0

/-- Get the maximum depth of the tree. -/
def maxDepth (tree : Octree) : Nat :=
  octreeDepth tree.root

/-- Check if the tree is empty. -/
def isEmpty (tree : Octree) : Bool := tree.itemCount == 0

/-- Helper for k-nearest search. -/
private partial def kNearestOctreeNode {α : Type} [HasPosition3D α] (node : OctreeNode) (items : Array α)
    (point : Vec3) (k : Nat) (heap : MaxHeap Nat) : MaxHeap Nat :=
  let nodeDist := node.bounds.distanceSquared point
  if heap.isFull && nodeDist > heap.maxDistance then heap
  else match node with
    | OctreeNode.leaf _ leafItems =>
      leafItems.foldl (fun h item =>
        let idx := item.idx
        if h_idx : idx < items.size then
          let pos := HasPosition3D.position items[idx]
          let dist := point.distanceSquared pos
          h.insert dist idx
        else h
      ) heap
    | OctreeNode.internal _ children =>
      -- Process children sorted by distance for better pruning
      let sortedChildren := children.mapIdx (fun i child => (i, child))
        |>.filter (·.2.isSome)
        |>.map (fun (i, c) =>
          let childBounds := OctreeNode.childBounds node.bounds (mkOctant i)
          (i, c.get!, childBounds.distanceSquared point))
        |>.qsort (fun a b => a.2.2 < b.2.2)
      sortedChildren.foldl (fun h (_, child, _) =>
        kNearestOctreeNode child items point k h
      ) heap

/-- Find k nearest items to a point. -/
def kNearest {α : Type} [HasPosition3D α] (tree : Octree) (items : Array α) (point : Vec3) (k : Nat) : Array Nat :=
  if k == 0 || tree.isEmpty then #[]
  else
    let heap := kNearestOctreeNode tree.root items point k (MaxHeap.empty k)
    heap.toSortedArray

/-- Remove an item from a node by index. -/
private partial def removeFromNode (node : OctreeNode) (idx : Nat) (itemBounds : AABB)
    : OctreeNode × Bool :=
  -- If item bounds don't intersect this node, nothing to do
  if !Intersection.aabbAABB node.bounds itemBounds then (node, false)
  else match node with
    | OctreeNode.leaf bounds items =>
      let filtered := items.filter (·.idx != idx)
      let removed := filtered.size < items.size
      (OctreeNode.leaf bounds filtered, removed)
    | OctreeNode.internal bounds children =>
      let (newChildren, removed) := children.foldl (fun (acc, didRemove) child =>
        match child with
        | some c =>
            let (newChild, childRemoved) := removeFromNode c idx itemBounds
            (acc.push (some newChild), didRemove || childRemoved)
        | none => (acc.push none, didRemove)
      ) (#[], false)
      (OctreeNode.internal bounds newChildren, removed)

/-- Remove an item from the octree by index. Requires item bounds for efficient lookup. -/
def remove (tree : Octree) (idx : Nat) (itemBounds : AABB) : Octree :=
  let (newRoot, removed) := removeFromNode tree.root idx itemBounds
  let newCount := if removed && tree.itemCount > 0 then tree.itemCount - 1 else tree.itemCount
  { tree with
    root := newRoot
    itemCount := newCount }

end Octree

end Linalg.Spatial
