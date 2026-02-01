/-
  Quadtree for 2D spatial partitioning.
  Recursively subdivides space into 4 quadrants.
-/

import Linalg.Spatial.Common

namespace Linalg.Spatial

/-- Quadrant index: 0=SW, 1=SE, 2=NW, 3=NE -/
abbrev QuadrantIndex := Fin 4

/-- Get quadrant for a point relative to center.
    Bit 0 = X >= center.x, Bit 1 = Y >= center.y -/
def quadrantFor (center point : Vec2) : QuadrantIndex :=
  let xBit : Fin 2 := if point.x >= center.x then 1 else 0
  let yBit : Fin 2 := if point.y >= center.y then 1 else 0
  ⟨xBit.val + yBit.val * 2, by have hx := xBit.isLt; have hy := yBit.isLt; omega⟩

/-- Stored item reference inside a quadtree leaf. -/
structure QuadtreeItem where
  idx : Nat
  bounds : AABB2D
deriving Repr, Inhabited

/-- Quadtree node. -/
inductive QuadtreeNode where
  /-- Internal node with 4 children (some may be empty). -/
  | internal (bounds : AABB2D) (children : Array (Option QuadtreeNode))
  /-- Leaf node containing item references. -/
  | leaf (bounds : AABB2D) (items : Array QuadtreeItem)
deriving Repr, Inhabited

namespace QuadtreeNode

/-- Get the bounds of a node. -/
def bounds : QuadtreeNode → AABB2D
  | internal b _ => b
  | leaf b _ => b

/-- Check if a node is a leaf. -/
def isLeaf : QuadtreeNode → Bool
  | internal _ _ => false
  | leaf _ _ => true

/-- Get child bounds for a given quadrant. -/
def childBounds (b : AABB2D) (q : QuadrantIndex) : AABB2D :=
  let c := b.center
  match q with
  | ⟨0, _⟩ => AABB2D.fromMinMax b.min c  -- SW
  | ⟨1, _⟩ => AABB2D.fromMinMax (Vec2.mk c.x b.min.y) (Vec2.mk b.max.x c.y)  -- SE
  | ⟨2, _⟩ => AABB2D.fromMinMax (Vec2.mk b.min.x c.y) (Vec2.mk c.x b.max.y)  -- NW
  | ⟨3, _⟩ => AABB2D.fromMinMax c b.max  -- NE

end QuadtreeNode

/-- Quadtree for 2D spatial partitioning. -/
structure Quadtree where
  /-- Root node of the tree. -/
  root : QuadtreeNode
  /-- Overall bounds. -/
  bounds : AABB2D
  /-- Total number of items. -/
  itemCount : Nat
  /-- Configuration. -/
  config : TreeConfig
deriving Repr, Inhabited

namespace Quadtree

/-- Create an empty quadtree with given bounds. -/
def empty (bounds : AABB2D) (config : TreeConfig := {}) : Quadtree :=
  { root := QuadtreeNode.leaf bounds #[]
    bounds
    itemCount := 0
    config }

/-- Helper to create a valid quadrant index from a Nat known to be < 4. -/
private def mkQuadrant (n : Nat) : QuadrantIndex :=
  ⟨n % 4, Nat.mod_lt n (by omega)⟩

/-- Insert an item into a node, potentially subdividing. -/
private partial def insertNode {α : Type} [Bounded2D α] (node : QuadtreeNode) (idx : Nat) (item : α)
    (depth : Nat) (config : TreeConfig) : QuadtreeNode :=
  let itemBounds := Bounded2D.bounds item
  let itemRef : QuadtreeItem := { idx, bounds := itemBounds }
  insertNodeRef node itemRef depth config
where
  insertNodeRef (node : QuadtreeNode) (item : QuadtreeItem) (depth : Nat) (config : TreeConfig) : QuadtreeNode :=
    match node with
    | QuadtreeNode.leaf bounds items =>
      let newItems := items.push item
      -- Check if we should split
      if newItems.size > config.maxLeafItems && depth < config.maxDepth &&
         bounds.width > config.minNodeSize && bounds.height > config.minNodeSize then
        -- Create internal node
        let children := #[none, none, none, none]
        let internalNode := QuadtreeNode.internal bounds children
        -- Re-insert all items into the new internal node
        newItems.foldl (fun n it => insertNodeRef n it depth config) internalNode
      else
        QuadtreeNode.leaf bounds newItems
    | QuadtreeNode.internal bounds children =>
      -- Find which quadrants the item overlaps
      let overlaps := #[0, 1, 2, 3].filter fun q =>
        let childB := QuadtreeNode.childBounds bounds (mkQuadrant q)
        childB.intersects item.bounds
      -- Insert into each overlapping quadrant
      let newChildren := overlaps.foldl (fun cs q =>
        let childBounds := QuadtreeNode.childBounds bounds (mkQuadrant q)
        let existingChild := cs[q]!
        let child := match existingChild with
          | some c => insertNodeRef c item (depth + 1) config
          | none => insertNodeRef (QuadtreeNode.leaf childBounds #[]) item (depth + 1) config
        cs.set! q (some child)
      ) children
      QuadtreeNode.internal bounds newChildren

/-- Insert a single item into the quadtree. -/
def insert {α : Type} [Bounded2D α] (tree : Quadtree) (idx : Nat) (item : α) : Quadtree :=
  { tree with
    root := insertNode tree.root idx item 0 tree.config
    itemCount := tree.itemCount + 1 }

/-- Build a quadtree from an array of items. -/
def build {α : Type} [Bounded2D α] [Inhabited α] (items : Array α) (config : TreeConfig := {}) : Quadtree :=
  if items.isEmpty then
    empty (AABB2D.fromMinMax Vec2.zero Vec2.one) config
  else
    -- Compute overall bounds with small padding
    let bounds := items.foldl (fun acc item =>
      AABB2D.merge acc (Bounded2D.bounds item)
    ) (Bounded2D.bounds items[0]!)
    -- Add small padding to avoid edge cases
    let padding := bounds.size.scale 0.01
    let paddedBounds := AABB2D.fromMinMax (bounds.min.sub padding) (bounds.max.add padding)
    -- Insert all items
    items.foldl (fun (tree, idx) item =>
      (tree.insert idx item, idx + 1)
    ) (empty paddedBounds config, 0) |>.1

/-- Query a node for items in a rectangle. -/
private partial def queryNodeRect (node : QuadtreeNode) (query : AABB2D) : Array Nat :=
  if !node.bounds.intersects query then #[]
  else match node with
    | QuadtreeNode.leaf _ items => items.map (·.idx)
    | QuadtreeNode.internal _ children =>
      children.foldl (fun acc child =>
        match child with
        | some c => acc ++ queryNodeRect c query
        | none => acc
      ) #[]

/-- Query candidate items that may intersect the given rectangle.
    Returns indices of items in leaf nodes that overlap the query region.
    Note: This is a spatial acceleration query - candidates may not actually
    intersect the rectangle. Use `queryRectExact` for precise intersection testing. -/
def queryRect (tree : Quadtree) (query : AABB2D) : Array Nat :=
  let results := queryNodeRect tree.root query
  -- Remove duplicates (items may be in multiple nodes)
  results.foldl (fun acc idx =>
    if acc.contains idx then acc else acc.push idx
  ) #[]

/-- Query items whose bounds actually intersect the given rectangle.
    Unlike `queryRect`, this verifies each candidate against its actual bounds. -/
def queryRectExact {α : Type} [Bounded2D α] (tree : Quadtree) (items : Array α) (query : AABB2D) : Array Nat :=
  let candidates := tree.queryRect query
  candidates.filter fun idx =>
    if h : idx < items.size then
      (Bounded2D.bounds items[idx]).intersects query
    else false

/-- Query a node for items in a circle. -/
private partial def queryNodeCircle (node : QuadtreeNode) (center : Vec2) (radiusSq : Float) : Array Nat :=
  -- Quick reject: check if circle's bounding box intersects node
  let radius := Float.sqrt radiusSq
  let circleBounds := AABB2D.fromCenterExtents center (Vec2.mk radius radius)
  if !node.bounds.intersects circleBounds then #[]
  else match node with
    | QuadtreeNode.leaf _ items => items.map (·.idx)
    | QuadtreeNode.internal _ children =>
      children.foldl (fun acc child =>
        match child with
        | some c => acc ++ queryNodeCircle c center radiusSq
        | none => acc
      ) #[]

/-- Query candidate items that may intersect the given circle.
    Returns indices of items in leaf nodes that overlap the circle's bounding box.
    Note: This is a spatial acceleration query - use for broad-phase collision detection. -/
def queryCircle (tree : Quadtree) (center : Vec2) (radius : Float) : Array Nat :=
  let results := queryNodeCircle tree.root center (radius * radius)
  results.foldl (fun acc idx =>
    if acc.contains idx then acc else acc.push idx
  ) #[]

/-- Query a node for items containing a point. -/
private partial def queryNodePoint (node : QuadtreeNode) (point : Vec2) : Array Nat :=
  if !node.bounds.containsPoint point then #[]
  else match node with
    | QuadtreeNode.leaf _ items => items.map (·.idx)
    | QuadtreeNode.internal bounds children =>
      let q := quadrantFor bounds.center point
      match children[q.val]! with
      | some child => queryNodePoint child point
      | none => #[]

/-- Query candidate items that may contain the given point.
    Returns indices of items stored in the leaf node containing the point.
    Note: This is a spatial acceleration query - candidates' actual bounds may not
    contain the point. Use `queryPointExact` for precise containment testing. -/
def queryPoint (tree : Quadtree) (point : Vec2) : Array Nat :=
  queryNodePoint tree.root point

/-- Query items whose bounds actually contain the given point.
    Unlike `queryPoint`, this verifies each candidate against its actual bounds. -/
def queryPointExact {α : Type} [Bounded2D α] (tree : Quadtree) (items : Array α) (point : Vec2) : Array Nat :=
  let candidates := tree.queryPoint point
  candidates.filter fun idx =>
    if h : idx < items.size then
      (Bounded2D.bounds items[idx]).containsPoint point
    else false

/-- Count nodes recursively. -/
private partial def countNodes : QuadtreeNode → Nat
  | QuadtreeNode.leaf _ _ => 1
  | QuadtreeNode.internal _ children =>
    1 + children.foldl (fun acc child =>
      match child with
      | some c => acc + countNodes c
      | none => acc
    ) 0

/-- Count total nodes in the tree. -/
def nodeCount (tree : Quadtree) : Nat :=
  countNodes tree.root

/-- Compute depth recursively. -/
private partial def treeDepth : QuadtreeNode → Nat
  | QuadtreeNode.leaf _ _ => 0
  | QuadtreeNode.internal _ children =>
    1 + children.foldl (fun acc child =>
      match child with
      | some c => Nat.max acc (treeDepth c)
      | none => acc
    ) 0

/-- Get the maximum depth of the tree. -/
def maxDepth (tree : Quadtree) : Nat :=
  treeDepth tree.root

/-- Check if the tree is empty. -/
def isEmpty (tree : Quadtree) : Bool := tree.itemCount == 0

/-- Helper for k-nearest search. -/
private partial def kNearestNode {α : Type} [HasPosition2D α] (node : QuadtreeNode) (items : Array α)
    (point : Vec2) (k : Nat) (heap : MaxHeap Nat) : MaxHeap Nat :=
  -- Check if node could contain closer points
  let nodeDist := node.bounds.distanceSquared point
  if heap.isFull && nodeDist > heap.maxDistance then heap
  else match node with
    | QuadtreeNode.leaf _ leafItems =>
      leafItems.foldl (fun h item =>
        let idx := item.idx
        if h_idx : idx < items.size then
          let pos := HasPosition2D.position items[idx]
          let dist := point.distanceSquared pos
          h.insert dist idx
        else h
      ) heap
    | QuadtreeNode.internal _ children =>
      -- Sort children by distance to point for better pruning
      let sortedChildren := children.mapIdx (fun i child => (i, child))
        |>.filter (·.2.isSome)
        |>.map (fun (i, c) => (i, c.get!, c.get!.bounds.distanceSquared point))
        |>.qsort (fun a b => a.2.2 < b.2.2)
      sortedChildren.foldl (fun h (_, child, _) =>
        kNearestNode child items point k h
      ) heap

/-- Find k nearest items to a point. -/
def kNearest {α : Type} [HasPosition2D α] (tree : Quadtree) (items : Array α) (point : Vec2) (k : Nat) : Array Nat :=
  if k == 0 || tree.isEmpty then #[]
  else
    -- Start with infinite radius, progressively shrink
    let heap := kNearestNode tree.root items point k (MaxHeap.empty k)
    heap.toSortedArray

/-- Remove an item from a node by index. -/
private partial def removeFromNode (node : QuadtreeNode) (idx : Nat) (itemBounds : AABB2D)
    : QuadtreeNode × Bool :=
  -- If item bounds don't intersect this node, nothing to do
  if !node.bounds.intersects itemBounds then (node, false)
  else match node with
    | QuadtreeNode.leaf bounds items =>
      let filtered := items.filter (·.idx != idx)
      let removed := filtered.size < items.size
      (QuadtreeNode.leaf bounds filtered, removed)
    | QuadtreeNode.internal bounds children =>
      let (newChildren, removed) := children.foldl (fun (acc, didRemove) child =>
        match child with
        | some c =>
            let (newChild, childRemoved) := removeFromNode c idx itemBounds
            (acc.push (some newChild), didRemove || childRemoved)
        | none => (acc.push none, didRemove)
      ) (#[], false)
      (QuadtreeNode.internal bounds newChildren, removed)

/-- Remove an item from the quadtree by index. Requires item bounds for efficient lookup. -/
def remove (tree : Quadtree) (idx : Nat) (itemBounds : AABB2D) : Quadtree :=
  let (newRoot, removed) := removeFromNode tree.root idx itemBounds
  let newCount := if removed && tree.itemCount > 0 then tree.itemCount - 1 else tree.itemCount
  { tree with
    root := newRoot
    itemCount := newCount }

end Quadtree

end Linalg.Spatial
