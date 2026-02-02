/-
  Dynamic AABB Tree (DBVH) for broad-phase collision detection.

  Provides incremental insert/remove/update for moving objects and AABB queries.
-/

import Linalg.Spatial.Common
import Linalg.Geometry.Intersection

namespace Linalg.Spatial

/-- Configuration for dynamic AABB tree. -/
structure DynamicAABBTreeConfig where
  /-- Extra margin added to leaf bounds to reduce update churn. -/
  fatMargin : Float := 0.1
deriving Repr, BEq, Inhabited

/-- Dynamic AABB tree node. -/
structure DynamicAABBNode where
  parent : Option Nat := none
  left : Option Nat := none
  right : Option Nat := none
  bounds : AABB := AABB.fromMinMax Vec3.zero Vec3.zero
  height : Nat := 0
  item : Option Nat := none
deriving Repr, Inhabited

namespace DynamicAABBNode

/-- Check if node is a leaf. -/
def isLeaf (n : DynamicAABBNode) : Bool :=
  n.left.isNone && n.right.isNone

end DynamicAABBNode

/-- Dynamic AABB tree. -/
structure DynamicAABBTree where
  nodes : Array DynamicAABBNode
  root : Option Nat
  freeList : Array Nat
  size : Nat
  config : DynamicAABBTreeConfig
deriving Repr, Inhabited

namespace DynamicAABBTree

/-- Create an empty dynamic tree. -/
def empty (config : DynamicAABBTreeConfig := {}) : DynamicAABBTree :=
  { nodes := #[], root := none, freeList := #[], size := 0, config }

/-- Number of allocated nodes (includes internal nodes). -/
def nodeCount (t : DynamicAABBTree) : Nat := t.size

/-- Get item stored at a leaf proxy. -/
def getItem (t : DynamicAABBTree) (proxy : Nat) : Option Nat :=
  if proxy < t.nodes.size then t.nodes[proxy]!.item else none

/-- Get bounds stored at a proxy. -/
def getBounds (t : DynamicAABBTree) (proxy : Nat) : Option AABB :=
  if proxy < t.nodes.size then some t.nodes[proxy]!.bounds else none

private def fatten (b : AABB) (margin : Float) : AABB :=
  let delta := Vec3.one.scale margin
  AABB.fromMinMax (b.min - delta) (b.max + delta)

private def allocNode (t : DynamicAABBTree) : DynamicAABBTree × Nat := Id.run do
  let mut nodes := t.nodes
  let mut free := t.freeList
  let mut idx := 0
  if free.isEmpty then
    idx := nodes.size
    nodes := nodes.push default
  else
    idx := free[free.size - 1]!
    free := free.pop
  let t := { t with nodes := nodes, freeList := free, size := t.size + 1 }
  return (t, idx)

private def freeNode (t : DynamicAABBTree) (idx : Nat) : DynamicAABBTree :=
  { t with freeList := t.freeList.push idx, size := t.size - 1 }

private def nodeHeight (t : DynamicAABBTree) (idx : Nat) : Nat :=
  t.nodes[idx]!.height

private def balance (t : DynamicAABBTree) (iA : Nat) : DynamicAABBTree × Nat := Id.run do
  let mut tree := t
  let mut nodes := tree.nodes
  let nodeA := nodes[iA]!
  if nodeA.isLeaf || nodeA.height < 2 then
    return (tree, iA)

  let iB := nodeA.left.get!
  let iC := nodeA.right.get!
  let nodeB := nodes[iB]!
  let nodeC := nodes[iC]!
  let balance := (Int.ofNat nodeC.height) - (Int.ofNat nodeB.height)

  if balance > 1 then
    let iF := nodeC.left.get!
    let iG := nodeC.right.get!
    let nodeF := nodes[iF]!
    let nodeG := nodes[iG]!

    -- Rotate left: C becomes parent of A
    let parent := nodeA.parent
    let nodeC' := { nodeC with parent := parent, left := some iA }
    nodes := nodes.set! iC nodeC'
    let nodeA' := { nodeA with parent := some iC }
    nodes := nodes.set! iA nodeA'

    match parent with
    | none =>
        tree := { tree with root := some iC }
    | some p =>
        let pnode := nodes[p]!
        let pnode' := if pnode.left == some iA then { pnode with left := some iC }
                      else { pnode with right := some iC }
        nodes := nodes.set! p pnode'

    if nodeF.height > nodeG.height then
      let nodeC'' := { nodes[iC]! with right := some iF }
      nodes := nodes.set! iC nodeC''
      let nodeA'' := { nodes[iA]! with right := some iG }
      nodes := nodes.set! iA nodeA''
      nodes := nodes.set! iF { nodeF with parent := some iC }
      nodes := nodes.set! iG { nodeG with parent := some iA }
    else
      let nodeC'' := { nodes[iC]! with right := some iG }
      nodes := nodes.set! iC nodeC''
      let nodeA'' := { nodes[iA]! with right := some iF }
      nodes := nodes.set! iA nodeA''
      nodes := nodes.set! iF { nodeF with parent := some iA }
      nodes := nodes.set! iG { nodeG with parent := some iC }

    let nodeA3 := nodes[iA]!
    let leftA := nodes[nodeA3.left.get!]!
    let rightA := nodes[nodeA3.right.get!]!
    let boundsA := AABB.merge leftA.bounds rightA.bounds
    let heightA := 1 + max leftA.height rightA.height
    nodes := nodes.set! iA { nodeA3 with bounds := boundsA, height := heightA }

    let nodeC3 := nodes[iC]!
    let leftC := nodes[nodeC3.left.get!]!
    let rightC := nodes[nodeC3.right.get!]!
    let boundsC := AABB.merge leftC.bounds rightC.bounds
    let heightC := 1 + max leftC.height rightC.height
    nodes := nodes.set! iC { nodeC3 with bounds := boundsC, height := heightC }

    tree := { tree with nodes := nodes }
    return (tree, iC)

  if balance < -1 then
    let iD := nodeB.left.get!
    let iE := nodeB.right.get!
    let nodeD := nodes[iD]!
    let nodeE := nodes[iE]!

    -- Rotate right: B becomes parent of A
    let parent := nodeA.parent
    let nodeB' := { nodeB with parent := parent, right := some iA }
    nodes := nodes.set! iB nodeB'
    let nodeA' := { nodeA with parent := some iB }
    nodes := nodes.set! iA nodeA'

    match parent with
    | none =>
        tree := { tree with root := some iB }
    | some p =>
        let pnode := nodes[p]!
        let pnode' := if pnode.left == some iA then { pnode with left := some iB }
                      else { pnode with right := some iB }
        nodes := nodes.set! p pnode'

    if nodeD.height > nodeE.height then
      let nodeB'' := { nodes[iB]! with left := some iD }
      nodes := nodes.set! iB nodeB''
      let nodeA'' := { nodes[iA]! with left := some iE }
      nodes := nodes.set! iA nodeA''
      nodes := nodes.set! iD { nodeD with parent := some iB }
      nodes := nodes.set! iE { nodeE with parent := some iA }
    else
      let nodeB'' := { nodes[iB]! with left := some iE }
      nodes := nodes.set! iB nodeB''
      let nodeA'' := { nodes[iA]! with left := some iD }
      nodes := nodes.set! iA nodeA''
      nodes := nodes.set! iD { nodeD with parent := some iA }
      nodes := nodes.set! iE { nodeE with parent := some iB }

    let nodeA3 := nodes[iA]!
    let leftA := nodes[nodeA3.left.get!]!
    let rightA := nodes[nodeA3.right.get!]!
    let boundsA := AABB.merge leftA.bounds rightA.bounds
    let heightA := 1 + max leftA.height rightA.height
    nodes := nodes.set! iA { nodeA3 with bounds := boundsA, height := heightA }

    let nodeB3 := nodes[iB]!
    let leftB := nodes[nodeB3.left.get!]!
    let rightB := nodes[nodeB3.right.get!]!
    let boundsB := AABB.merge leftB.bounds rightB.bounds
    let heightB := 1 + max leftB.height rightB.height
    nodes := nodes.set! iB { nodeB3 with bounds := boundsB, height := heightB }

    tree := { tree with nodes := nodes }
    return (tree, iB)

  return (tree, iA)

private def fixUpwards (t : DynamicAABBTree) (start : Option Nat) : DynamicAABBTree := Id.run do
  let mut tree := t
  let mut current := start
  while current.isSome do
    let idx := current.get!
    let (tree', newIdx) := balance tree idx
    tree := tree'
    let node := tree.nodes[newIdx]!
    match node.left, node.right with
    | some l, some r =>
        let left := tree.nodes[l]!
        let right := tree.nodes[r]!
        let bounds := AABB.merge left.bounds right.bounds
        let height := 1 + max left.height right.height
        let node' := { node with bounds := bounds, height := height }
        tree := { tree with nodes := tree.nodes.set! newIdx node' }
        current := node'.parent
    | _, _ =>
        current := node.parent
  return tree

private def insertLeaf (t : DynamicAABBTree) (leaf : Nat) : DynamicAABBTree := Id.run do
  let mut tree := t
  let mut nodes := tree.nodes
  if tree.root.isNone then
    let node := nodes[leaf]!
    nodes := nodes.set! leaf { node with parent := none }
    tree := { tree with nodes := nodes, root := some leaf }
    return tree

  let leafBounds := nodes[leaf]!.bounds
  let mut index := tree.root.get!
  while !(nodes[index]!).isLeaf do
    let node := nodes[index]!
    let leftIdx := node.left.get!
    let rightIdx := node.right.get!
    let left := nodes[leftIdx]!
    let right := nodes[rightIdx]!
    let combined := AABB.merge node.bounds leafBounds
    let cost := 2.0 * combined.surfaceArea
    let inheritanceCost := 2.0 * (combined.surfaceArea - node.bounds.surfaceArea)
    let leftMerge := AABB.merge left.bounds leafBounds
    let rightMerge := AABB.merge right.bounds leafBounds
    let costLeft := (if left.isLeaf then leftMerge.surfaceArea
                     else leftMerge.surfaceArea - left.bounds.surfaceArea) + inheritanceCost
    let costRight := (if right.isLeaf then rightMerge.surfaceArea
                      else rightMerge.surfaceArea - right.bounds.surfaceArea) + inheritanceCost
    if cost < costLeft && cost < costRight then
      break
    index := if costLeft < costRight then leftIdx else rightIdx

  let sibling := index
  let parent := nodes[sibling]!.parent
  let (tree', newParentIdx) := allocNode tree
  tree := tree'
  nodes := tree.nodes
  let newParent := {
    parent := parent
    left := some sibling
    right := some leaf
    bounds := AABB.merge nodes[sibling]!.bounds leafBounds
    height := nodes[sibling]!.height + 1
    item := none
  }
  nodes := nodes.set! newParentIdx newParent
  nodes := nodes.set! sibling { nodes[sibling]! with parent := some newParentIdx }
  nodes := nodes.set! leaf { nodes[leaf]! with parent := some newParentIdx }
  match parent with
  | none =>
      tree := { tree with root := some newParentIdx }
  | some p =>
      let pnode := nodes[p]!
      let pnode' := if pnode.left == some sibling then { pnode with left := some newParentIdx }
                    else { pnode with right := some newParentIdx }
      nodes := nodes.set! p pnode'
  tree := { tree with nodes := nodes }
  tree := fixUpwards tree (some newParentIdx)
  return tree

private def removeLeaf (t : DynamicAABBTree) (leaf : Nat) : DynamicAABBTree := Id.run do
  let mut tree := t
  let mut nodes := tree.nodes
  if tree.root == some leaf then
    let node := nodes[leaf]!
    nodes := nodes.set! leaf { node with parent := none }
    tree := { tree with nodes := nodes, root := none }
    return tree

  let parent := nodes[leaf]!.parent.get!
  let grand := nodes[parent]!.parent
  let sibling := if nodes[parent]!.left == some leaf then nodes[parent]!.right.get!
                 else nodes[parent]!.left.get!

  match grand with
  | none =>
      nodes := nodes.set! sibling { nodes[sibling]! with parent := none }
      tree := { tree with nodes := nodes, root := some sibling }
  | some g =>
      let gnode := nodes[g]!
      let gnode' := if gnode.left == some parent then { gnode with left := some sibling }
                    else { gnode with right := some sibling }
      nodes := nodes.set! g gnode'
      nodes := nodes.set! sibling { nodes[sibling]! with parent := some g }
      tree := { tree with nodes := nodes }
      tree := fixUpwards tree (some g)

  tree := freeNode tree parent
  return tree

/-- Insert a new item with bounds. Returns updated tree and proxy id. -/
def insert (t : DynamicAABBTree) (bounds : AABB) (item : Nat) : DynamicAABBTree × Nat := Id.run do
  let (tree, proxy) := allocNode t
  let fat := fatten bounds tree.config.fatMargin
  let node := { bounds := fat, height := 0, item := some item }
  let tree := { tree with nodes := tree.nodes.set! proxy node }
  let tree := insertLeaf tree proxy
  return (tree, proxy)

/-- Remove an item by proxy id. -/
def remove (t : DynamicAABBTree) (proxy : Nat) : DynamicAABBTree := Id.run do
  if proxy >= t.nodes.size then
    return t
  let mut tree := removeLeaf t proxy
  let mut nodes := tree.nodes
  let node := nodes[proxy]!
  nodes := nodes.set! proxy { node with parent := none, left := none, right := none, item := none, height := 0 }
  tree := { tree with nodes := nodes }
  tree := freeNode tree proxy
  return tree

/-- Update a proxy with new bounds, reinserting if necessary. -/
def update (t : DynamicAABBTree) (proxy : Nat) (bounds : AABB) : DynamicAABBTree := Id.run do
  if proxy >= t.nodes.size then
    return t
  let node := t.nodes[proxy]!
  let fat := fatten bounds t.config.fatMargin
  if AABB.containsAABB node.bounds fat then
    return t
  let mut tree := removeLeaf t proxy
  let mut nodes := tree.nodes
  nodes := nodes.set! proxy { node with bounds := fat }
  tree := { tree with nodes := nodes }
  tree := insertLeaf tree proxy
  return tree

/-- Query the tree for items whose bounds intersect the query AABB. -/
def queryAABB (t : DynamicAABBTree) (query : AABB) : Array Nat := Id.run do
  let mut results : Array Nat := #[]
  match t.root with
  | none => return results
  | some root =>
      let mut stack : Array Nat := #[root]
      while !stack.isEmpty do
        let idx := stack[stack.size - 1]!
        stack := stack.pop
        let node := t.nodes[idx]!
        if !Intersection.aabbAABB node.bounds query then
          continue
        if node.isLeaf then
          match node.item with
          | some item => results := results.push item
          | none => ()
        else
          match node.left with
          | some l => stack := stack.push l
          | none => ()
          match node.right with
          | some r => stack := stack.push r
          | none => ()
      return results

/-- Query the tree for items hit by a ray (broad-phase). -/
def queryRay (t : DynamicAABBTree) (ray : Ray) (maxT : Float := Float.infinity) : Array Nat := Id.run do
  let mut results : Array Nat := #[]
  match t.root with
  | none => return results
  | some root =>
      let mut stack : Array Nat := #[root]
      while !stack.isEmpty do
        let idx := stack[stack.size - 1]!
        stack := stack.pop
        let node := t.nodes[idx]!
        match Intersection.rayAABB ray node.bounds with
        | none => continue
        | some (tMin, tMax) =>
            if tMin > maxT || tMax < 0.0 then
              continue
            if node.isLeaf then
              match node.item with
              | some item => results := results.push item
              | none => ()
            else
              match node.left with
              | some l => stack := stack.push l
              | none => ()
              match node.right with
              | some r => stack := stack.push r
              | none => ()
      return results

/-- Collect broad-phase overlap pairs (item id pairs). -/
def broadPhasePairs (t : DynamicAABBTree) : Array (Nat × Nat) := Id.run do
  let mut pairs : Array (Nat × Nat) := #[]
  for i in [:t.nodes.size] do
    let node := t.nodes[i]!
    if node.isLeaf then
      match node.item with
      | none => ()
      | some item =>
          let hits := t.queryAABB node.bounds
          for other in hits do
            if item < other then
              pairs := pairs.push (item, other)
  return pairs

end DynamicAABBTree

end Linalg.Spatial
