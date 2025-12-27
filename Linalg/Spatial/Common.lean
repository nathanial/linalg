/-
  Common types for spatial data structures.
  Provides typeclasses for extracting bounds/positions and configuration structs.
-/

import Linalg.Vec2
import Linalg.Vec3
import Linalg.Geometry.AABB
import Linalg.Geometry.AABB2D
import Linalg.Geometry.Ray

namespace Linalg.Spatial

/-- Typeclass for items that have a 3D bounding box. -/
class Bounded3D (α : Type) where
  bounds : α → AABB

/-- Typeclass for items that have a 2D bounding box. -/
class Bounded2D (α : Type) where
  bounds : α → AABB2D

/-- Typeclass for items that have a point position in 3D. -/
class HasPosition3D (α : Type) where
  position : α → Vec3

/-- Typeclass for items that have a point position in 2D. -/
class HasPosition2D (α : Type) where
  position : α → Vec2

-- Instance: Vec3 has position (itself)
instance : HasPosition3D Vec3 where
  position := id

-- Instance: Vec2 has position (itself)
instance : HasPosition2D Vec2 where
  position := id

-- Instance: AABB has bounds (itself)
instance : Bounded3D AABB where
  bounds := id

-- Instance: AABB2D has bounds (itself)
instance : Bounded2D AABB2D where
  bounds := id

-- Instance: Vec3 can be treated as a point-sized AABB
instance : Bounded3D Vec3 where
  bounds := AABB.fromPoint

-- Instance: Vec2 can be treated as a point-sized AABB2D
instance : Bounded2D Vec2 where
  bounds := AABB2D.fromPoint

/-- Configuration for tree-based spatial structures (Quadtree, Octree). -/
structure TreeConfig where
  /-- Maximum tree depth to prevent infinite subdivision. -/
  maxDepth : Nat := 16
  /-- Maximum items per leaf before splitting. -/
  maxLeafItems : Nat := 8
  /-- Minimum node dimension (stop subdividing below this size). -/
  minNodeSize : Float := 0.001
deriving Repr, BEq, Inhabited

namespace TreeConfig

/-- Default configuration. -/
def default : TreeConfig := {}

/-- Configuration for dense scenes (more items per leaf, less depth). -/
def dense : TreeConfig := { maxDepth := 12, maxLeafItems := 16 }

/-- Configuration for sparse scenes (fewer items per leaf, more depth). -/
def sparse : TreeConfig := { maxDepth := 20, maxLeafItems := 4 }

end TreeConfig

/-- Configuration for BVH construction. -/
structure BVHConfig where
  /-- Maximum items per leaf node. -/
  maxLeafItems : Nat := 4
  /-- Number of buckets for SAH evaluation. -/
  sahBuckets : Nat := 12
  /-- Cost of traversing a node (relative to intersection). -/
  traversalCost : Float := 1.0
  /-- Cost of testing intersection with an item. -/
  intersectionCost : Float := 1.0
deriving Repr, BEq, Inhabited

namespace BVHConfig

/-- Default configuration. -/
def default : BVHConfig := {}

/-- Configuration optimized for triangles (common in ray tracing). -/
def triangles : BVHConfig := { maxLeafItems := 4, sahBuckets := 16 }

/-- Configuration optimized for fewer, larger objects. -/
def objects : BVHConfig := { maxLeafItems := 1, sahBuckets := 8 }

end BVHConfig

/-- Configuration for uniform grids. -/
structure GridConfig where
  /-- Size of each cell. If none, computed from bounds. -/
  cellSize : Option Float := none
  /-- Target number of cells along the longest axis. Used if cellSize is none. -/
  targetCells : Nat := 32
deriving Repr, BEq, Inhabited

/-- Result of a spatial ray query. -/
structure SpatialHit where
  /-- Index of the hit item. -/
  index : Nat
  /-- Distance along ray. -/
  t : Float
  /-- Hit point. -/
  point : Vec3
  /-- Surface normal at hit. -/
  normal : Vec3
deriving Repr, BEq, Inhabited

/-- Result of a 2D spatial query (for Quadtree). -/
structure SpatialHit2D where
  /-- Index of the hit item. -/
  index : Nat
  /-- Distance to item. -/
  distance : Float
  /-- Closest point on item. -/
  point : Vec2
deriving Repr, BEq, Inhabited

/-- A simple priority queue for k-nearest neighbor queries. -/
structure MaxHeap (α : Type) where
  items : Array (Float × α)
  maxSize : Nat
deriving Repr, Inhabited

namespace MaxHeap

variable {α : Type}

/-- Create an empty max heap with given capacity. -/
def empty (maxSize : Nat) : MaxHeap α := { items := #[], maxSize }

/-- Check if the heap is full. -/
def isFull (h : MaxHeap α) : Bool := h.items.size >= h.maxSize

/-- Get the maximum distance in the heap. -/
def maxDistance (h : MaxHeap α) : Float :=
  if h.items.isEmpty then Float.infinity
  else h.items.foldl (fun acc (d, _) => Float.max acc d) 0.0

/-- Insert an item with given distance. If full and distance >= max, ignore. -/
def insert (h : MaxHeap α) (distance : Float) (item : α) : MaxHeap α :=
  if h.isFull then
    -- Find and replace the maximum if new distance is smaller
    let maxIdx := h.items.foldl (fun (bestIdx, bestDist, idx) (d, _) =>
      if d > bestDist then (idx, d, idx + 1) else (bestIdx, bestDist, idx + 1)
    ) (0, 0.0, 0)
    if distance < maxIdx.2.1 then
      { h with items := h.items.set! maxIdx.1 (distance, item) }
    else
      h
  else
    { h with items := h.items.push (distance, item) }

/-- Get all items sorted by distance (ascending). -/
def toSortedArray (h : MaxHeap α) : Array α :=
  let sorted := h.items.qsort (fun (d1, _) (d2, _) => d1 < d2)
  sorted.map (fun (_, item) => item)

end MaxHeap

end Linalg.Spatial
