/-
  Intersection tests between geometric primitives.
-/

import Linalg.Geometry.Ray
import Linalg.Geometry.AABB
import Linalg.Geometry.Sphere
import Linalg.Geometry.Plane
import Linalg.Geometry.Triangle
import Linalg.Geometry.OBB
import Linalg.Geometry.Capsule

namespace Linalg

/-- Result of a ray intersection test. -/
structure RayHit where
  t : Float        -- Parameter along ray (distance if direction is normalized)
  point : Vec3     -- Intersection point
  normal : Vec3    -- Surface normal at intersection
deriving Repr, Inhabited

namespace Intersection

/-- Ray-Sphere intersection. Returns the nearest hit if any. -/
def raySphere (ray : Ray) (sphere : Sphere) : Option RayHit :=
  let oc := ray.origin.sub sphere.center
  let a := ray.direction.dot ray.direction
  let halfB := oc.dot ray.direction
  let c := oc.dot oc - sphere.radius * sphere.radius
  let discriminant := halfB * halfB - a * c

  if discriminant < 0.0 then none
  else
    let sqrtD := Float.sqrt discriminant
    -- Try the nearer root first
    let t := (-halfB - sqrtD) / a
    if t < 0.0 then
      -- Try the farther root (ray origin inside sphere)
      let t2 := (-halfB + sqrtD) / a
      if t2 < 0.0 then none
      else
        let point := ray.pointAt t2
        let normal := (point.sub sphere.center).normalize
        some ⟨t2, point, normal⟩
    else
      let point := ray.pointAt t
      let normal := (point.sub sphere.center).normalize
      some ⟨t, point, normal⟩

/-- Ray-AABB intersection using the slab method.
    Returns (tMin, tMax) interval if intersection exists. -/
def rayAABB (ray : Ray) (aabb : AABB) : Option (Float × Float) :=
  let invDx := 1.0 / ray.direction.x
  let invDy := 1.0 / ray.direction.y
  let invDz := 1.0 / ray.direction.z

  let tx1 := (aabb.min.x - ray.origin.x) * invDx
  let tx2 := (aabb.max.x - ray.origin.x) * invDx
  let ty1 := (aabb.min.y - ray.origin.y) * invDy
  let ty2 := (aabb.max.y - ray.origin.y) * invDy
  let tz1 := (aabb.min.z - ray.origin.z) * invDz
  let tz2 := (aabb.max.z - ray.origin.z) * invDz

  let tMinX := Float.min tx1 tx2
  let tMaxX := Float.max tx1 tx2
  let tMinY := Float.min ty1 ty2
  let tMaxY := Float.max ty1 ty2
  let tMinZ := Float.min tz1 tz2
  let tMaxZ := Float.max tz1 tz2

  let tMin := Float.max tMinX (Float.max tMinY tMinZ)
  let tMax := Float.min tMaxX (Float.min tMaxY tMaxZ)

  if tMax < 0.0 || tMin > tMax then none
  else some (tMin, tMax)

/-- Ray-AABB intersection returning a RayHit. -/
def rayAABBHit (ray : Ray) (aabb : AABB) : Option RayHit :=
  match rayAABB ray aabb with
  | none => none
  | some (tMin, _tMax) =>
    let t := if tMin >= 0.0 then tMin else 0.0  -- Clamp to origin if inside
    let point := ray.pointAt t
    -- Calculate normal based on which face was hit
    let center := aabb.center
    let extents := aabb.extents
    let localPt := point.sub center
    let nx := localPt.x / extents.x
    let ny := localPt.y / extents.y
    let nz := localPt.z / extents.z
    let ax := Float.abs' nx
    let ay := Float.abs' ny
    let az := Float.abs' nz
    let normal :=
      if ax > ay && ax > az then
        if nx > 0.0 then Vec3.unitX else Vec3.left
      else if ay > az then
        if ny > 0.0 then Vec3.unitY else Vec3.down
      else
        if nz > 0.0 then Vec3.unitZ else Vec3.forward
    some ⟨t, point, normal⟩

/-- Ray-Plane intersection. -/
def rayPlane (ray : Ray) (plane : Plane) : Option RayHit :=
  let denom := plane.normal.dot ray.direction
  if Float.abs' denom < Float.epsilon then none  -- Ray parallel to plane
  else
    let t := (plane.distance - plane.normal.dot ray.origin) / denom
    if t < 0.0 then none  -- Intersection behind ray origin
    else
      let point := ray.pointAt t
      -- Normal faces toward the ray
      let normal := if denom < 0.0 then plane.normal else plane.normal.neg
      some ⟨t, point, normal⟩

/-- Sphere-Sphere intersection test. -/
def sphereSphere (a b : Sphere) : Bool :=
  let distSq := a.center.distanceSquared b.center
  let radiusSum := a.radius + b.radius
  distSq <= radiusSum * radiusSum

/-- Sphere-Sphere penetration depth (negative if not intersecting). -/
def sphereSpherePenetration (a b : Sphere) : Float :=
  let dist := a.center.distance b.center
  (a.radius + b.radius) - dist

/-- AABB-AABB intersection test. -/
def aabbAABB (a b : AABB) : Bool :=
  a.min.x <= b.max.x && a.max.x >= b.min.x &&
  a.min.y <= b.max.y && a.max.y >= b.min.y &&
  a.min.z <= b.max.z && a.max.z >= b.min.z

/-- Sphere-AABB intersection test. -/
def sphereAABB (sphere : Sphere) (aabb : AABB) : Bool :=
  -- Find closest point on AABB to sphere center
  let closest := aabb.closestPoint sphere.center
  sphere.center.distanceSquared closest <= sphere.radius * sphere.radius

/-- Sphere-Plane intersection test (checks if sphere intersects plane). -/
def spherePlane (sphere : Sphere) (plane : Plane) : Bool :=
  plane.distanceToPoint sphere.center <= sphere.radius

/-- Sphere-Plane classification.
    Returns: -1 = entirely behind, 0 = intersecting, 1 = entirely in front -/
def spherePlaneClassify (sphere : Sphere) (plane : Plane) : Int :=
  let dist := plane.signedDistance sphere.center
  if dist > sphere.radius then 1
  else if dist < -sphere.radius then -1
  else 0

/-- AABB-Plane intersection test. -/
def aabbPlane (aabb : AABB) (plane : Plane) : Bool :=
  -- Project half-extents onto plane normal
  let center := aabb.center
  let extents := aabb.extents
  let r := extents.x * Float.abs' plane.normal.x +
           extents.y * Float.abs' plane.normal.y +
           extents.z * Float.abs' plane.normal.z
  Float.abs' (plane.signedDistance center) <= r

/-- Point containment tests (convenience re-exports). -/
def pointInSphere (p : Vec3) (s : Sphere) : Bool := s.containsPoint p
def pointInAABB (p : Vec3) (b : AABB) : Bool := b.containsPoint p
def pointOnPlane (p : Vec3) (plane : Plane) (eps : Float := Float.epsilon) : Bool :=
  plane.containsPoint p eps

/-- Closest points between two line segments.
    Returns (point on segment 1, point on segment 2, squared distance). -/
def closestPointsSegmentSegment (p1 q1 p2 q2 : Vec3) : Vec3 × Vec3 × Float :=
  let d1 := q1.sub p1
  let d2 := q2.sub p2
  let r := p1.sub p2
  let a := d1.dot d1
  let e := d2.dot d2
  let f := d2.dot r

  -- Check if either or both segments degenerate to points
  if a < Float.epsilon && e < Float.epsilon then
    -- Both segments degenerate to points
    (p1, p2, p1.distanceSquared p2)
  else if a < Float.epsilon then
    -- First segment degenerates to point
    let t := Float.clamp (f / e) 0.0 1.0
    let closestOnSeg2 := p2.add (d2.scale t)
    (p1, closestOnSeg2, p1.distanceSquared closestOnSeg2)
  else
    let c := d1.dot r
    if e < Float.epsilon then
      -- Second segment degenerates to point
      let s := Float.clamp (-c / a) 0.0 1.0
      let closestOnSeg1 := p1.add (d1.scale s)
      (closestOnSeg1, p2, closestOnSeg1.distanceSquared p2)
    else
      -- General case
      let b := d1.dot d2
      let denom := a * e - b * b

      let s := if Float.abs' denom >= Float.epsilon then
        Float.clamp ((b * f - c * e) / denom) 0.0 1.0
      else 0.0

      let t := (b * s + f) / e
      let (s', t') :=
        if t < 0.0 then
          (Float.clamp (-c / a) 0.0 1.0, 0.0)
        else if t > 1.0 then
          (Float.clamp ((b - c) / a) 0.0 1.0, 1.0)
        else (s, t)

      let closestOnSeg1 := p1.add (d1.scale s')
      let closestOnSeg2 := p2.add (d2.scale t')
      (closestOnSeg1, closestOnSeg2, closestOnSeg1.distanceSquared closestOnSeg2)

/-- Ray-Triangle intersection using the Moller-Trumbore algorithm.
    Returns hit info including barycentric coordinates. -/
def rayTriangle (ray : Ray) (tri : Triangle) (cullBackface : Bool := false) : Option RayHit := Id.run do
  let edge1 := tri.edge01
  let edge2 := tri.edge02
  let h := ray.direction.cross edge2
  let a := edge1.dot h

  -- Check if ray is parallel to triangle
  if Float.abs' a < Float.epsilon then
    return none

  -- Backface culling (if enabled)
  if cullBackface && a < 0.0 then
    return none

  let f := 1.0 / a
  let s := ray.origin.sub tri.v0
  let u := f * s.dot h

  -- Check u bounds
  if u < 0.0 || u > 1.0 then
    return none

  let q := s.cross edge1
  let v := f * ray.direction.dot q

  -- Check v bounds and u+v <= 1
  if v < 0.0 || u + v > 1.0 then
    return none

  -- Compute t to find intersection point
  let t := f * edge2.dot q

  if t < Float.epsilon then
    return none  -- Line intersection but not ray intersection

  let point := ray.pointAt t
  let normal := tri.unitNormal
  -- Flip normal to face ray if needed
  let normal := if ray.direction.dot normal > 0.0 then normal.neg else normal
  return some ⟨t, point, normal⟩

/-- Ray-Triangle hit with barycentric coordinates. -/
structure TriangleHit where
  t : Float
  point : Vec3
  normal : Vec3
  u : Float  -- barycentric u
  v : Float  -- barycentric v
  w : Float  -- barycentric w (= 1 - u - v)
  deriving Repr

/-- Ray-Triangle intersection returning barycentric coordinates.
    Useful for texture mapping and vertex attribute interpolation. -/
def rayTriangleBarycentric (ray : Ray) (tri : Triangle) (cullBackface : Bool := false) : Option TriangleHit := Id.run do
  let edge1 := tri.edge01
  let edge2 := tri.edge02
  let h := ray.direction.cross edge2
  let a := edge1.dot h

  if Float.abs' a < Float.epsilon then
    return none

  if cullBackface && a < 0.0 then
    return none

  let f := 1.0 / a
  let s := ray.origin.sub tri.v0
  let u := f * s.dot h

  if u < 0.0 || u > 1.0 then
    return none

  let q := s.cross edge1
  let v := f * ray.direction.dot q

  if v < 0.0 || u + v > 1.0 then
    return none

  let t := f * edge2.dot q

  if t < Float.epsilon then
    return none

  let point := ray.pointAt t
  let normal := tri.unitNormal
  let normal := if ray.direction.dot normal > 0.0 then normal.neg else normal
  let w := 1.0 - u - v
  -- Note: u,v,w here are barycentric coords where point = u*v0 + v*v1 + w*v2
  return some { t := t, point := point, normal := normal, u := w, v := u, w := v }

/-- Ray-OBB intersection.
    Transforms the ray to OBB local space and uses slab method. -/
def rayOBB (ray : Ray) (obb : OBB) : Option RayHit := Id.run do
  -- Transform ray to OBB local space
  let localOrigin := obb.worldToLocal ray.origin
  let localDir := obb.orientation.inverse.rotateVec3 ray.direction

  -- Use slab method in local space (axis-aligned)
  let mut tMin := Float.negInfinity
  let mut tMax := Float.infinity
  let mut normalIndex := 0
  let mut normalSign := 1.0

  -- X slab
  if Float.abs' localDir.x < Float.epsilon then
    if localOrigin.x < -obb.halfExtents.x || localOrigin.x > obb.halfExtents.x then
      return none
  else
    let invD := 1.0 / localDir.x
    let mut t1 := (-obb.halfExtents.x - localOrigin.x) * invD
    let mut t2 := (obb.halfExtents.x - localOrigin.x) * invD
    let mut sign := -1.0
    if t1 > t2 then
      let tmp := t1; t1 := t2; t2 := tmp; sign := 1.0
    if t1 > tMin then
      tMin := t1; normalIndex := 0; normalSign := sign
    if t2 < tMax then
      tMax := t2
    if tMin > tMax then return none

  -- Y slab
  if Float.abs' localDir.y < Float.epsilon then
    if localOrigin.y < -obb.halfExtents.y || localOrigin.y > obb.halfExtents.y then
      return none
  else
    let invD := 1.0 / localDir.y
    let mut t1 := (-obb.halfExtents.y - localOrigin.y) * invD
    let mut t2 := (obb.halfExtents.y - localOrigin.y) * invD
    let mut sign := -1.0
    if t1 > t2 then
      let tmp := t1; t1 := t2; t2 := tmp; sign := 1.0
    if t1 > tMin then
      tMin := t1; normalIndex := 1; normalSign := sign
    if t2 < tMax then
      tMax := t2
    if tMin > tMax then return none

  -- Z slab
  if Float.abs' localDir.z < Float.epsilon then
    if localOrigin.z < -obb.halfExtents.z || localOrigin.z > obb.halfExtents.z then
      return none
  else
    let invD := 1.0 / localDir.z
    let mut t1 := (-obb.halfExtents.z - localOrigin.z) * invD
    let mut t2 := (obb.halfExtents.z - localOrigin.z) * invD
    let mut sign := -1.0
    if t1 > t2 then
      let tmp := t1; t1 := t2; t2 := tmp; sign := 1.0
    if t1 > tMin then
      tMin := t1; normalIndex := 2; normalSign := sign
    if t2 < tMax then
      tMax := t2
    if tMin > tMax then return none

  -- Check if intersection is behind ray
  if tMax < 0.0 then return none

  let t := if tMin >= 0.0 then tMin else tMax
  let point := ray.pointAt t

  -- Compute normal in world space
  let localNormal := match normalIndex with
    | 0 => Vec3.mk normalSign 0.0 0.0
    | 1 => Vec3.mk 0.0 normalSign 0.0
    | _ => Vec3.mk 0.0 0.0 normalSign
  let normal := obb.orientation.rotateVec3 localNormal

  return some ⟨t, point, normal⟩

/-- Ray-Capsule intersection.
    Tests ray against the cylindrical body and hemispherical caps. -/
def rayCapsule (ray : Ray) (capsule : Capsule) : Option RayHit := Id.run do
  let ab := capsule.segment
  let ao := ray.origin.sub capsule.a

  let abab := ab.dot ab
  let abao := ab.dot ao
  let abrd := ab.dot ray.direction

  -- Coefficients for quadratic in t: at² + bt + c = 0
  let a := abab - abrd * abrd
  let b := abab * (ao.dot ray.direction) - abao * abrd
  let c := abab * (ao.dot ao) - abao * abao - capsule.radius * capsule.radius * abab

  let mut bestT := Float.infinity
  let mut bestNormal := Vec3.zero

  -- Check infinite cylinder
  if Float.abs' a > Float.epsilon then
    let discriminant := b * b - a * c
    if discriminant >= 0.0 then
      let sqrtD := Float.sqrt discriminant
      for sign in #[-1.0, 1.0] do
        let t := (-b + sign * sqrtD) / a
        if t >= 0.0 && t < bestT then
          -- Check if hit is within the cylinder body (not caps)
          let hitPoint := ray.pointAt t
          let projection := (hitPoint.sub capsule.a).dot ab / abab
          if projection >= 0.0 && projection <= 1.0 then
            -- Hit is on the cylinder body
            let axisPoint := capsule.a.add (ab.scale projection)
            let normal := (hitPoint.sub axisPoint).normalize
            bestT := t
            bestNormal := normal

  -- Check hemisphere at endpoint A
  let ocA := ray.origin.sub capsule.a
  let aA := ray.direction.dot ray.direction
  let halfBA := ocA.dot ray.direction
  let cA := ocA.dot ocA - capsule.radius * capsule.radius
  let discA := halfBA * halfBA - aA * cA
  if discA >= 0.0 then
    let sqrtD := Float.sqrt discA
    for sign in #[-1.0, 1.0] do
      let t := (-halfBA + sign * sqrtD) / aA
      if t >= 0.0 && t < bestT then
        let hitPoint := ray.pointAt t
        -- Check that the hit is on the hemisphere (not beyond the cylinder)
        let toHit := hitPoint.sub capsule.a
        if toHit.dot ab <= 0.0 then
          let normal := toHit.normalize
          bestT := t
          bestNormal := normal

  -- Check hemisphere at endpoint B
  let ocB := ray.origin.sub capsule.b
  let halfBB := ocB.dot ray.direction
  let cB := ocB.dot ocB - capsule.radius * capsule.radius
  let discB := halfBB * halfBB - aA * cB
  if discB >= 0.0 then
    let sqrtD := Float.sqrt discB
    for sign in #[-1.0, 1.0] do
      let t := (-halfBB + sign * sqrtD) / aA
      if t >= 0.0 && t < bestT then
        let hitPoint := ray.pointAt t
        -- Check that the hit is on the hemisphere (not beyond the cylinder)
        let toHit := hitPoint.sub capsule.b
        if toHit.dot ab >= 0.0 then
          let normal := toHit.normalize
          bestT := t
          bestNormal := normal

  if bestT < Float.infinity then
    return some ⟨bestT, ray.pointAt bestT, bestNormal⟩
  else
    return none

/-- OBB-Sphere intersection test. -/
def obbSphere (obb : OBB) (sphere : Sphere) : Bool :=
  obb.intersectsSphere sphere.center sphere.radius

/-- Capsule-Sphere intersection test. -/
def capsuleSphere (capsule : Capsule) (sphere : Sphere) : Bool :=
  capsule.intersectsSphere sphere

/-- Capsule-Capsule intersection test. -/
def capsuleCapsule (c1 c2 : Capsule) : Bool :=
  c1.intersectsCapsule c2

/-- OBB-OBB intersection test. -/
def obbOBB (a b : OBB) : Bool :=
  a.intersectsOBB b

/-- OBB-AABB intersection test. -/
def obbAABB (obb : OBB) (aabb : AABB) : Bool :=
  obb.intersectsAABB aabb

/-- Capsule-AABB intersection test. -/
def capsuleAABB (capsule : Capsule) (aabb : AABB) : Bool :=
  capsule.intersectsAABB aabb

/-- Helper: Project triangle vertices onto an axis and return (min, max). -/
private def projectTriangle (tri : Triangle) (axis : Vec3) : Float × Float :=
  let p0 := tri.v0.dot axis
  let p1 := tri.v1.dot axis
  let p2 := tri.v2.dot axis
  (Float.min p0 (Float.min p1 p2), Float.max p0 (Float.max p1 p2))

/-- Helper: Check if projections on an axis are separated. -/
private def isSeparatedOnAxis (tri1 tri2 : Triangle) (axis : Vec3) : Bool :=
  -- Skip degenerate axes (near-zero length)
  if axis.lengthSquared < Float.epsilon then false
  else
    let (min1, max1) := projectTriangle tri1 axis
    let (min2, max2) := projectTriangle tri2 axis
    max1 < min2 || max2 < min1

/-- Helper: Check if a point is on the positive side of an edge (for coplanar test).
    Edge from v0 to v1, point p, normal n. -/
private def pointOnPositiveSide (v0 v1 p n : Vec3) : Bool :=
  let edge := v1.sub v0
  let edgeNormal := edge.cross n
  let toPoint := p.sub v0
  toPoint.dot edgeNormal >= 0.0

/-- Helper: Check if point is inside triangle (for coplanar case). -/
private def pointInTriangleCoplanar (p : Vec3) (tri : Triangle) (n : Vec3) : Bool :=
  pointOnPositiveSide tri.v0 tri.v1 p n &&
  pointOnPositiveSide tri.v1 tri.v2 p n &&
  pointOnPositiveSide tri.v2 tri.v0 p n

/-- Helper: Check if two 2D segments intersect (for coplanar triangle test).
    Projects 3D points to 2D by dropping the axis with largest normal component. -/
private def segmentsIntersect2D (a0 a1 b0 b1 : Vec3) (dropAxis : Nat) : Bool :=
  -- Project to 2D by dropping one coordinate
  let proj := fun (v : Vec3) =>
    match dropAxis with
    | 0 => (v.y, v.z)  -- Drop X
    | 1 => (v.x, v.z)  -- Drop Y
    | _ => (v.x, v.y)  -- Drop Z

  let (ax0, ay0) := proj a0
  let (ax1, ay1) := proj a1
  let (bx0, by0) := proj b0
  let (bx1, by1) := proj b1

  -- 2D cross product
  let cross2D := fun (ux uy vx vy : Float) => ux * vy - uy * vx

  let d1x := ax1 - ax0; let d1y := ay1 - ay0
  let d2x := bx1 - bx0; let d2y := by1 - by0

  let cross := cross2D d1x d1y d2x d2y

  if Float.abs' cross < Float.epsilon then
    -- Parallel segments - could check for overlap but rare case, skip for now
    false
  else
    let dx := bx0 - ax0; let dy := by0 - ay0
    let t := cross2D dx dy d2x d2y / cross
    let u := cross2D dx dy d1x d1y / cross
    t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0

/-- Helper: Get the axis to drop for 2D projection (largest normal component). -/
private def getDropAxis (n : Vec3) : Nat :=
  let ax := Float.abs' n.x
  let ay := Float.abs' n.y
  let az := Float.abs' n.z
  if ax >= ay && ax >= az then 0
  else if ay >= az then 1
  else 2

/-- Helper: Check if coplanar triangles intersect. -/
private def coplanarTrianglesIntersect (tri1 tri2 : Triangle) (n : Vec3) : Bool :=
  let dropAxis := getDropAxis n

  -- Check if any vertex of tri2 is inside tri1
  if pointInTriangleCoplanar tri2.v0 tri1 n then true
  else if pointInTriangleCoplanar tri2.v1 tri1 n then true
  else if pointInTriangleCoplanar tri2.v2 tri1 n then true
  -- Check if any vertex of tri1 is inside tri2
  else if pointInTriangleCoplanar tri1.v0 tri2 n then true
  else if pointInTriangleCoplanar tri1.v1 tri2 n then true
  else if pointInTriangleCoplanar tri1.v2 tri2 n then true
  else
    -- Check edge-edge intersections
    let edges1 := #[(tri1.v0, tri1.v1), (tri1.v1, tri1.v2), (tri1.v2, tri1.v0)]
    let edges2 := #[(tri2.v0, tri2.v1), (tri2.v1, tri2.v2), (tri2.v2, tri2.v0)]
    Id.run do
      for (a0, a1) in edges1 do
        for (b0, b1) in edges2 do
          if segmentsIntersect2D a0 a1 b0 b1 dropAxis then
            return true
      return false

/-- Triangle-Triangle intersection test using the Separating Axis Theorem.
    Tests 11 potential separating axes:
    - 2 triangle normals
    - 9 edge-edge cross products (3 edges × 3 edges)
    Also handles coplanar triangles. -/
def triangleTriangle (tri1 tri2 : Triangle) : Bool := Id.run do
  let n1 := tri1.normal
  let n2 := tri2.normal

  -- Get edges (edge20 = v0 - v2 = -edge02)
  let e1 := #[tri1.edge01, tri1.edge12, tri1.edge02.neg]
  let e2 := #[tri2.edge01, tri2.edge12, tri2.edge02.neg]

  -- Test triangle normals as separating axes
  if isSeparatedOnAxis tri1 tri2 n1 then return false
  if isSeparatedOnAxis tri1 tri2 n2 then return false

  -- Test edge-edge cross products
  for edge1 in e1 do
    for edge2 in e2 do
      let axis := edge1.cross edge2
      if isSeparatedOnAxis tri1 tri2 axis then
        return false

  -- Check if triangles are coplanar
  let d1 := n1.dot (tri2.v0.sub tri1.v0)
  let d2 := n1.dot (tri2.v1.sub tri1.v0)
  let d3 := n1.dot (tri2.v2.sub tri1.v0)

  -- If all vertices of tri2 are on the same side of tri1's plane and close to it
  if Float.abs' d1 < Float.epsilon && Float.abs' d2 < Float.epsilon && Float.abs' d3 < Float.epsilon then
    -- Coplanar case
    return coplanarTrianglesIntersect tri1 tri2 n1

  -- No separating axis found - triangles intersect
  return true

/-- Triangle-Triangle intersection with contact information.
    Returns the intersection as a line segment if they intersect. -/
def triangleTriangleContact (tri1 tri2 : Triangle) : Option (Vec3 × Vec3) := Id.run do
  if !triangleTriangle tri1 tri2 then
    return none

  -- Find the intersection line by computing where tri2's plane intersects tri1's edges
  let n2 := tri2.unitNormal
  let d2 := n2.dot tri2.v0

  let signedDist := fun (v : Vec3) => n2.dot v - d2

  let d0 := signedDist tri1.v0
  let d1 := signedDist tri1.v1
  let d2v := signedDist tri1.v2

  let mut points : Array Vec3 := #[]

  -- Check edge v0-v1
  if d0 * d1 < 0.0 then
    let t := d0 / (d0 - d1)
    points := points.push (Vec3.lerp tri1.v0 tri1.v1 t)

  -- Check edge v1-v2
  if d1 * d2v < 0.0 then
    let t := d1 / (d1 - d2v)
    points := points.push (Vec3.lerp tri1.v1 tri1.v2 t)

  -- Check edge v2-v0
  if d2v * d0 < 0.0 then
    let t := d2v / (d2v - d0)
    points := points.push (Vec3.lerp tri1.v2 tri1.v0 t)

  if points.size >= 2 then
    return some (points[0]!, points[1]!)
  else
    -- Coplanar or touching case - return centroid as single point
    let center := (tri1.centroid.add tri2.centroid).scale 0.5
    return some (center, center)

end Intersection

end Linalg
