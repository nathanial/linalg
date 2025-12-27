# Linalg Roadmap

Potential features and enhancements for the linalg library.

## Geometric Primitives - COMPLETED

- [x] **Triangle** - Three vertices with barycentric coordinate support
- [x] **OBB (Oriented Bounding Box)** - Center, half-extents, and rotation
- [x] **Capsule** - Line segment with radius (useful for character collision)
- [x] **Frustum** - Six planes for view frustum culling
- [x] **Polygon2D** - 2D polygon with winding order and area calculation
- [x] **Circle** - 2D circle primitive
- [x] **Line2D/Segment2D** - Infinite line and finite segment types

## Intersection Tests

- [x] **Ray-Triangle** - Moller-Trumbore algorithm
- [x] **Ray-OBB** - Oriented bounding box intersection
- [x] **Ray-Capsule** - Capsule intersection for character selection
- [x] **Sphere-Plane** - Sphere vs plane intersection/classification
- [x] **AABB-Plane** - Box vs plane classification
- [x] **Frustum-AABB** - Frustum culling test
- [x] **Frustum-Sphere** - Sphere visibility test
- [ ] **Triangle-Triangle** - Mesh collision detection
- [x] **Point-Triangle** - Point in triangle test (2D and 3D)
- [x] **Closest point queries** - Closest point on line/plane/triangle to point
- [x] **OBB-OBB** - Oriented bounding box collision (SAT)
- [x] **OBB-Sphere** - OBB vs sphere intersection
- [x] **Capsule-Capsule** - Capsule vs capsule collision
- [x] **Capsule-Sphere** - Capsule vs sphere intersection

## Curves and Splines

- [x] **Bezier curves** - Quadratic and cubic Bezier
- [x] **Catmull-Rom splines** - Smooth interpolation through control points
- [ ] **B-splines** - Basis splines for smooth curves
- [x] **Arc-length parameterization** - Constant-speed traversal
- [ ] **Bezier patches** - 2D surface patches

## Animation and Interpolation

- [x] **Easing functions** - EaseIn, EaseOut, EaseInOut variants (quad, cubic, elastic, bounce)
- [x] **Spring interpolation** - Damped spring for smooth UI motion
- [ ] **SmoothDamp** - Unity-style smooth following
- [ ] **Hermite interpolation** - Position and tangent interpolation
- [ ] **Animation curves** - Keyframe-based value curves

## Transform Utilities

- [x] **Transform type** - Combined position, rotation, scale
- [ ] **Transform hierarchy** - Parent-child relationships with local/world space
- [ ] **Dual quaternions** - Rigid body transforms (rotation + translation)
- [x] **Matrix decomposition** - Extract TRS from arbitrary matrix
- [ ] **Euler angle utilities** - Different rotation orders (XYZ, ZYX, etc.)
- [x] **LookAt quaternion** - Quaternion that looks at target

## 2D Game Math

- [x] **Affine2D (Mat2x3)** - 2D affine transform matrix
- [x] **Rotation2D** - 2D rotation type (angle-based)
- [ ] **Polygon operations** - Convex hull, triangulation, point-in-polygon
- [ ] **SAT (Separating Axis Theorem)** - 2D convex collision detection
- [ ] **GJK algorithm** - General convex shape intersection

## Spatial Structures

- [ ] **BVH (Bounding Volume Hierarchy)** - Acceleration structure for ray tracing
- [ ] **Octree** - 3D spatial partitioning
- [ ] **Quadtree** - 2D spatial partitioning
- [ ] **KD-Tree** - K-dimensional spatial queries
- [ ] **Grid** - Uniform spatial grid

## Noise and Random

- [x] **Perlin noise** - Classic gradient noise (1D, 2D, 3D)
- [x] **Simplex noise** - Improved gradient noise
- [x] **Fractal noise** - FBM, turbulence, ridged multifractal
- [x] **Domain warping** - Distort coordinates with noise for organic patterns
- [ ] **Worley noise** - Cellular/Voronoi noise
- [ ] **Value noise** - Simple interpolated noise
- [ ] **Random in shapes** - Random point in circle, sphere, cone, etc.
- [ ] **PCG random** - Fast, high-quality PRNG

## Physics Helpers

- [ ] **Velocity integration** - Euler, Verlet, RK4
- [ ] **Collision response** - Impulse-based resolution
- [ ] **Inertia tensors** - Common shape inertia calculations
- [ ] **Swept collision** - Moving shape vs static shape
- [ ] **Continuous collision** - Time of impact calculation

## Color (Optional)

- [ ] **Color3/Color4** - RGB/RGBA color types
- [ ] **HSV/HSL conversion** - Color space transforms
- [ ] **Color blending** - Standard blend modes
- [ ] **Gradient** - Linear and radial gradients

*Note: May integrate with `tincture` library instead*

## Performance

- [ ] **SIMD operations** - FFI to native SIMD (NEON on ARM, SSE on x86)
- [ ] **Batch operations** - Process arrays of vectors/matrices
- [ ] **Memory layout** - SoA (Structure of Arrays) variants
- [ ] **Compile-time evaluation** - More `@[reducible]` annotations

## API Improvements

- [ ] **Swizzling** - Vector component swizzling (e.g., `v.xyz`, `v.xxy`)
- [ ] **Named constructors** - `Vec3.right`, `Vec3.forward`, `Vec3.up`
- [ ] **Conversion instances** - `Coe Vec2 Vec3`, etc.
- [ ] **Format instances** - Pretty printing with configurable precision
- [ ] **Decidable equality** - For vectors with exact float comparison
- [ ] **JSON serialization** - ToJson/FromJson instances

## Documentation

- [ ] **API docs** - Docstrings for all public functions
- [ ] **Tutorial** - Getting started guide
- [ ] **Examples** - Common game math recipes
- [ ] **Comparison** - Performance comparison with other implementations

## Integration

- [ ] **Afferent integration** - Replace afferent's internal math types
- [ ] **Protobuf codecs** - Serialize vectors/matrices to protobuf
- [ ] **Ledger attributes** - Store geometry in fact database

---

## Priority Suggestions

### High Priority (Most Useful) - COMPLETED
1. ~~Triangle primitive and ray-triangle intersection~~
2. ~~Frustum for view culling~~
3. ~~Easing functions for animation~~
4. ~~Transform type combining position/rotation/scale~~

### Medium Priority - COMPLETED
5. ~~Bezier curves~~ ✓
6. ~~OBB and capsule primitives~~ ✓
7. ~~Perlin/Simplex noise~~ ✓
8. ~~2D affine transforms~~ ✓

### Nice to Have
9. Spatial structures (BVH, octree)
10. SIMD optimization
11. Physics helpers
