/-
  Sweep-and-Prune broad-phase for dynamic scenes.

  Maintains a sorted list of AABB intervals along a chosen axis and provides
  overlap pair generation and queries.
-/

import Std.Data.HashMap
import Linalg.Spatial.Common
import Linalg.Geometry.Intersection

namespace Linalg.Spatial

/-- Axis used for sweep-and-prune ordering. -/
inductive SweepAxis
  | x
  | y
  | z
deriving Repr, BEq, Inhabited

/-- Sweep-and-prune entry with cached interval. -/
structure SAPEntry where
  id : Nat
  min : Float
  max : Float
deriving Repr, Inhabited

/-- Sweep-and-prune broad-phase structure. -/
structure SweepAndPrune where
  entries : Array SAPEntry
  aabbs : Std.HashMap Nat AABB
  axis : SweepAxis
deriving Repr, Inhabited

namespace SweepAndPrune

/-- Create an empty sweep-and-prune structure. -/
def empty (axis : SweepAxis := .x) : SweepAndPrune :=
  { entries := #[], aabbs := {}, axis }

private def axisValue (v : Vec3) (axis : SweepAxis) : Float :=
  match axis with
  | .x => v.x
  | .y => v.y
  | .z => v.z

private def interval (b : AABB) (axis : SweepAxis) : Float × Float :=
  (axisValue b.min axis, axisValue b.max axis)

private def sortEntries (entries : Array SAPEntry) : Array SAPEntry :=
  entries.qsort (fun a b => a.min < b.min)

private def rebuildEntries (t : SweepAndPrune) (axis : SweepAxis) : Array SAPEntry := Id.run do
  let mut entries : Array SAPEntry := #[]
  for e in t.entries do
    match Std.HashMap.get? t.aabbs e.id with
    | some b =>
        let (min, max) := interval b axis
        entries := entries.push { id := e.id, min := min, max := max }
    | none => ()
  return sortEntries entries

private def overlapCount (entries : Array SAPEntry) : Nat := Id.run do
  let mut count := 0
  let mut active : Array SAPEntry := #[]
  for e in entries do
    active := active.filter (fun a => a.max >= e.min)
    count := count + active.size
    active := active.push e
  return count

private def chooseAxis (t : SweepAndPrune) : SweepAxis :=
  if t.entries.isEmpty then
    t.axis
  else
    let entriesX := rebuildEntries t .x
    let entriesY := rebuildEntries t .y
    let entriesZ := rebuildEntries t .z
    let countX := overlapCount entriesX
    let countY := overlapCount entriesY
    let countZ := overlapCount entriesZ
    if countX <= countY && countX <= countZ then .x
    else if countY <= countZ then .y
    else .z

/-- Rebuild with the given axis (re-sorts entries). -/
def rebuildAxis (t : SweepAndPrune) (axis : SweepAxis) : SweepAndPrune :=
  { t with entries := rebuildEntries t axis, axis := axis }

/-- Rebuild with an automatically chosen axis (minimizes overlap candidates). -/
def rebuildAutoAxis (t : SweepAndPrune) : SweepAndPrune :=
  let axis := chooseAxis t
  t.rebuildAxis axis

/-- Insert an AABB with item id. -/
def insert (t : SweepAndPrune) (id : Nat) (bounds : AABB) : SweepAndPrune :=
  let (min, max) := interval bounds t.axis
  let entry : SAPEntry := { id := id, min := min, max := max }
  let entries := sortEntries (t.entries.push entry)
  let aabbs := t.aabbs.insert id bounds
  { t with entries := entries, aabbs := aabbs }

/-- Remove an item by id. -/
def remove (t : SweepAndPrune) (id : Nat) : SweepAndPrune :=
  let entries := t.entries.filter (fun e => e.id != id)
  let aabbs := t.aabbs.erase id
  { t with entries := entries, aabbs := aabbs }

/-- Update an item's bounds (re-sorts entries). -/
def update (t : SweepAndPrune) (id : Nat) (bounds : AABB) : SweepAndPrune :=
  let (min, max) := interval bounds t.axis
  let entries := t.entries.map (fun e => if e.id == id then { e with min := min, max := max } else e)
  let entries := sortEntries entries
  let aabbs := t.aabbs.insert id bounds
  { t with entries := entries, aabbs := aabbs }

/-- Query items whose bounds intersect the given AABB. -/
def queryAABB (t : SweepAndPrune) (query : AABB) : Array Nat := Id.run do
  let (qmin, qmax) := interval query t.axis
  let mut results : Array Nat := #[]
  for e in t.entries do
    if e.max >= qmin && e.min <= qmax then
      match Std.HashMap.get? t.aabbs e.id with
      | some b =>
          if Intersection.aabbAABB b query then
            results := results.push e.id
      | none => ()
  return results

/-- Generate overlap pairs (item id pairs). -/
def broadPhasePairs (t : SweepAndPrune) : Array (Nat × Nat) := Id.run do
  let mut pairs : Array (Nat × Nat) := #[]
  let mut active : Array SAPEntry := #[]
  for e in t.entries do
    active := active.filter (fun a => a.max >= e.min)
    for a in active do
      match Std.HashMap.get? t.aabbs a.id, Std.HashMap.get? t.aabbs e.id with
      | some ba, some bb =>
          if Intersection.aabbAABB ba bb then
            if a.id < e.id then
              pairs := pairs.push (a.id, e.id)
            else if e.id < a.id then
              pairs := pairs.push (e.id, a.id)
      | _, _ => ()
    active := active.push e
  return pairs

/-- Generate overlap pairs using an auto-selected axis. -/
def broadPhasePairsAuto (t : SweepAndPrune) : Array (Nat × Nat) :=
  (t.rebuildAutoAxis).broadPhasePairs

/-- Query with an auto-selected axis. -/
def queryAABBAuto (t : SweepAndPrune) (query : AABB) : Array Nat :=
  (t.rebuildAutoAxis).queryAABB query

end SweepAndPrune

end Linalg.Spatial
