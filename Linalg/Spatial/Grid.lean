/-
  Uniform spatial grid for O(1) cell lookup.
  Supports both 2D and 3D variants.
-/

import Linalg.Spatial.Common
import Std.Data.HashMap

namespace Linalg.Spatial

/-- 3D cell coordinate. -/
structure CellCoord3D where
  x : Int
  y : Int
  z : Int
deriving Repr, BEq, Inhabited, Hashable

/-- 2D cell coordinate. -/
structure CellCoord2D where
  x : Int
  y : Int
deriving Repr, BEq, Inhabited, Hashable

/-- Convert Float to Int by flooring. -/
private def floatToInt (f : Float) : Int :=
  Int.ofNat f.floor.toUInt64.toNat

/-- 3D uniform spatial grid. -/
structure Grid3D where
  /-- Size of each cell. -/
  cellSize : Float
  /-- Origin for cell coordinate calculation. -/
  origin : Vec3
  /-- Map from cell coordinate to item indices. -/
  cells : Std.HashMap CellCoord3D (Array Nat)
  /-- Total number of items inserted. -/
  itemCount : Nat
deriving Repr, Inhabited

/-- 2D uniform spatial grid. -/
structure Grid2D where
  /-- Size of each cell. -/
  cellSize : Float
  /-- Origin for cell coordinate calculation. -/
  origin : Vec2
  /-- Map from cell coordinate to item indices. -/
  cells : Std.HashMap CellCoord2D (Array Nat)
  /-- Total number of items inserted. -/
  itemCount : Nat
deriving Repr, Inhabited

namespace Grid3D

/-- Create an empty 3D grid. -/
def empty (cellSize : Float) (origin : Vec3 := Vec3.zero) : Grid3D :=
  { cellSize, origin, cells := {}, itemCount := 0 }

/-- Get the cell coordinate for a point. -/
def cellFor (g : Grid3D) (p : Vec3) : CellCoord3D :=
  let rel := p.sub g.origin
  { x := floatToInt (rel.x / g.cellSize),
    y := floatToInt (rel.y / g.cellSize),
    z := floatToInt (rel.z / g.cellSize) }

/-- Get the cell coordinate range for an AABB. -/
def cellRangeFor (g : Grid3D) (b : AABB) : CellCoord3D × CellCoord3D :=
  (g.cellFor b.min, g.cellFor b.max)

/-- Get all items in a specific cell. -/
def itemsInCell (g : Grid3D) (cell : CellCoord3D) : Array Nat :=
  g.cells.getD cell #[]

/-- Helper to iterate over a range of integers. -/
private def intRange (start finish : Int) : List Int :=
  if start > finish then []
  else (List.range ((finish - start).toNat + 1)).map (fun i => start + Int.ofNat i)

/-- Insert an item index into the grid based on its bounds. -/
def insert {α : Type} [Bounded3D α] (g : Grid3D) (idx : Nat) (item : α) : Grid3D :=
  let bounds := Bounded3D.bounds item
  let (minCell, maxCell) := g.cellRangeFor bounds
  -- Insert into all cells the item overlaps
  let cells := (intRange minCell.x maxCell.x).foldl (fun cells x =>
    (intRange minCell.y maxCell.y).foldl (fun cells y =>
      (intRange minCell.z maxCell.z).foldl (fun cells z =>
        let coord : CellCoord3D := ⟨x, y, z⟩
        let existing := cells.getD coord #[]
        cells.insert coord (existing.push idx)
      ) cells
    ) cells
  ) g.cells
  { g with cells, itemCount := g.itemCount + 1 }

/-- Remove an item from the grid. -/
def remove (g : Grid3D) (idx : Nat) (itemBounds : AABB) : Grid3D :=
  let (minCell, maxCell) := g.cellRangeFor itemBounds
  let cells := (intRange minCell.x maxCell.x).foldl (fun cells x =>
    (intRange minCell.y maxCell.y).foldl (fun cells y =>
      (intRange minCell.z maxCell.z).foldl (fun cells z =>
        let coord : CellCoord3D := ⟨x, y, z⟩
        match cells.get? coord with
        | some arr =>
          let filtered := arr.filter (· != idx)
          if filtered.isEmpty then cells.erase coord
          else cells.insert coord filtered
        | none => cells
      ) cells
    ) cells
  ) g.cells
  { g with cells, itemCount := g.itemCount - 1 }

/-- Build a grid from an array of items. -/
def build {α : Type} [Bounded3D α] (items : Array α) (cellSize : Float) (origin : Vec3 := Vec3.zero) : Grid3D :=
  items.foldl (fun (g, idx) item => (g.insert idx item, idx + 1)) (empty cellSize origin, 0) |>.1

/-- Build a grid with automatic cell size based on bounds. -/
def buildAuto {α : Type} [Bounded3D α] [Inhabited α] (items : Array α) (config : GridConfig := {}) : Grid3D :=
  if items.isEmpty then empty 1.0
  else
    -- Compute overall bounds
    let bounds := items.foldl (fun acc item =>
      AABB.merge acc (Bounded3D.bounds item)
    ) (Bounded3D.bounds items[0]!)
    let size := bounds.size
    let maxDim := Float.max size.x (Float.max size.y size.z)
    let cellSize := match config.cellSize with
      | some cs => cs
      | none => maxDim / config.targetCells.toFloat
    build items cellSize bounds.min

/-- Query all items in cells intersecting the given AABB. Returns unique indices. -/
def queryAABB (g : Grid3D) (query : AABB) : Array Nat :=
  let (minCell, maxCell) := g.cellRangeFor query
  let indices := (intRange minCell.x maxCell.x).foldl (fun acc x =>
    (intRange minCell.y maxCell.y).foldl (fun acc y =>
      (intRange minCell.z maxCell.z).foldl (fun acc z =>
        let coord : CellCoord3D := ⟨x, y, z⟩
        let cellItems := g.itemsInCell coord
        cellItems.foldl (fun acc idx => acc.push idx) acc
      ) acc
    ) acc
  ) #[]
  -- Remove duplicates
  indices.foldl (fun acc idx =>
    if acc.contains idx then acc else acc.push idx
  ) #[]

/-- Query all items in cells intersecting the given sphere. Returns unique indices. -/
def querySphere (g : Grid3D) (center : Vec3) (radius : Float) : Array Nat :=
  -- Use AABB query as broad phase
  let aabb := AABB.fromCenterExtents center (Vec3.mk radius radius radius)
  queryAABB g aabb

/-- Get the number of non-empty cells. -/
def cellCount (g : Grid3D) : Nat := g.cells.size

/-- Check if the grid is empty. -/
def isEmpty (g : Grid3D) : Bool := g.itemCount == 0

end Grid3D

namespace Grid2D

/-- Create an empty 2D grid. -/
def empty (cellSize : Float) (origin : Vec2 := Vec2.zero) : Grid2D :=
  { cellSize, origin, cells := {}, itemCount := 0 }

/-- Get the cell coordinate for a point. -/
def cellFor (g : Grid2D) (p : Vec2) : CellCoord2D :=
  let rel := p.sub g.origin
  { x := floatToInt (rel.x / g.cellSize),
    y := floatToInt (rel.y / g.cellSize) }

/-- Get the cell coordinate range for an AABB2D. -/
def cellRangeFor (g : Grid2D) (b : AABB2D) : CellCoord2D × CellCoord2D :=
  (g.cellFor b.min, g.cellFor b.max)

/-- Get all items in a specific cell. -/
def itemsInCell (g : Grid2D) (cell : CellCoord2D) : Array Nat :=
  g.cells.getD cell #[]

/-- Helper to iterate over a range of integers. -/
private def intRange (start finish : Int) : List Int :=
  if start > finish then []
  else (List.range ((finish - start).toNat + 1)).map (fun i => start + Int.ofNat i)

/-- Insert an item index into the grid based on its bounds. -/
def insert {α : Type} [Bounded2D α] (g : Grid2D) (idx : Nat) (item : α) : Grid2D :=
  let bounds := Bounded2D.bounds item
  let (minCell, maxCell) := g.cellRangeFor bounds
  let cells := (intRange minCell.x maxCell.x).foldl (fun cells x =>
    (intRange minCell.y maxCell.y).foldl (fun cells y =>
      let coord : CellCoord2D := ⟨x, y⟩
      let existing := cells.getD coord #[]
      cells.insert coord (existing.push idx)
    ) cells
  ) g.cells
  { g with cells, itemCount := g.itemCount + 1 }

/-- Remove an item from the grid. -/
def remove (g : Grid2D) (idx : Nat) (itemBounds : AABB2D) : Grid2D :=
  let (minCell, maxCell) := g.cellRangeFor itemBounds
  let cells := (intRange minCell.x maxCell.x).foldl (fun cells x =>
    (intRange minCell.y maxCell.y).foldl (fun cells y =>
      let coord : CellCoord2D := ⟨x, y⟩
      match cells.get? coord with
      | some arr =>
        let filtered := arr.filter (· != idx)
        if filtered.isEmpty then cells.erase coord
        else cells.insert coord filtered
      | none => cells
    ) cells
  ) g.cells
  { g with cells, itemCount := g.itemCount - 1 }

/-- Build a grid from an array of items. -/
def build {α : Type} [Bounded2D α] (items : Array α) (cellSize : Float) (origin : Vec2 := Vec2.zero) : Grid2D :=
  items.foldl (fun (g, idx) item => (g.insert idx item, idx + 1)) (empty cellSize origin, 0) |>.1

/-- Build a grid with automatic cell size based on bounds. -/
def buildAuto {α : Type} [Bounded2D α] [Inhabited α] (items : Array α) (config : GridConfig := {}) : Grid2D :=
  if items.isEmpty then empty 1.0
  else
    -- Compute overall bounds
    let bounds := items.foldl (fun acc item =>
      AABB2D.merge acc (Bounded2D.bounds item)
    ) (Bounded2D.bounds items[0]!)
    let size := bounds.size
    let maxDim := Float.max size.x size.y
    let cellSize := match config.cellSize with
      | some cs => cs
      | none => maxDim / config.targetCells.toFloat
    build items cellSize bounds.min

/-- Query all items in cells intersecting the given AABB2D. Returns unique indices. -/
def queryRect (g : Grid2D) (query : AABB2D) : Array Nat :=
  let (minCell, maxCell) := g.cellRangeFor query
  let indices := (intRange minCell.x maxCell.x).foldl (fun acc x =>
    (intRange minCell.y maxCell.y).foldl (fun acc y =>
      let coord : CellCoord2D := ⟨x, y⟩
      let cellItems := g.itemsInCell coord
      cellItems.foldl (fun acc idx => acc.push idx) acc
    ) acc
  ) #[]
  -- Remove duplicates
  indices.foldl (fun acc idx =>
    if acc.contains idx then acc else acc.push idx
  ) #[]

/-- Query all items in cells intersecting the given circle. Returns unique indices. -/
def queryCircle (g : Grid2D) (center : Vec2) (radius : Float) : Array Nat :=
  -- Use rect query as broad phase
  let aabb := AABB2D.fromCenterExtents center (Vec2.mk radius radius)
  queryRect g aabb

/-- Get the number of non-empty cells. -/
def cellCount (g : Grid2D) : Nat := g.cells.size

/-- Check if the grid is empty. -/
def isEmpty (g : Grid2D) : Bool := g.itemCount == 0

end Grid2D

end Linalg.Spatial
