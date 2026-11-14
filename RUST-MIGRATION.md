# DGraphX Rust Migration Plan

## Overview

Rewrite DGraphX in Rust, preserving and enhancing core math, 3D/2D transformation, and rendering logic. Modernization goals: better architecture, testability, cross-platform support, GPU acceleration option.

## Core Philosophy

*   **Software Rendering:** Avoid heavy GPU-bound engines. Keep the logic for rasterization, 3D projection, and pixel manipulation in the application code.
*   **Zero-Dependency Logic:** The core math and rendering algorithms should remain independent of external math libraries where possible.
*   **Modern Ergonomics:** Replace manual X11 event loops with a clean, message-passing architecture.

## Proposed Tech Stack

| Component | C Original | Rust Modernization | Reason |
| :--- | :--- | :--- | :--- |
| **Language** | C (C99) | **Rust (Edition 2024)** | Memory safety, better abstractions. |
| **Build System** | Makefile | **Cargo** | Dependency management, cross-compilation. |
| **Windowing** | Xlib (Linux only) | **winit** | Cross-platform (Linux, Windows, macOS, Web). |
| **Framebuffer** | XCreatePixmap | **softbuffer** | Direct pixel access without GPU overhead. |
| **Math** | `matrix.c` | **glam** or **nalgebra** | SIMD-accelerated math (or pure Rust port). |
| **Concurrency** | Single-threaded | **rayon** | Parallelize surface calculation (3D plots). |
| **Scripting** | Hardcoded C functions | **evalexpr** | Runtime function input (already in C version via tinyexpr). |

---

## Phase 1: Core Math Library

**Goal**: Extract and modernize math foundations. Fully testable without rendering.

```
dgraphx-math/
├── Cargo.toml
├── src/
│   ├── lib.rs
│   ├── matrix.rs          # Matrix4x4, operations, inverse
│   ├── vector.rs          # Vec2, Vec3, Vec4 with SIMD potential
│   ├── transform.rs       # Rotation, translation, scale, projection
│   ├── interpolation.rs   # Hermite, spline (2D)
│   ├── solve.rs           # Seidel-Gauss, Richardson iteration
│   └── dtoa.rs           # Number to string (modernize)
└── tests/
```

**Deliverables**:
- 100% test coverage on all algorithms
- Benchmarks vs C implementation
- No external math dependencies (pure Rust)

---

## Phase 2: Software Renderer

**Goal**: Port X11 rasterization to backend-agnostic software renderer.

```
dgraphx-renderer/
├── src/
│   ├── lib.rs
│   ├── color.rs           # Color types, palette
│   ├── canvas.rs          # Pixel buffer, 2D primitives
│   ├── rasterizer.rs      # Scanline rasterization with z-buffer
│   ├── rasterizer_2d.rs   # 2D graph rasterization
│   └── font.rs            # Bitmap font rendering
```

**Key decisions**:
- Use `image` crate for pixel buffer
- Export to PNG for testing (no display required)
- Eventually swap to GPU rendering

---

## Phase 3: Windowing & Event Loop

**Goal**: Replace X11 with cross-platform windowing.

```
dgraphx-core/
├── src/
│   ├── lib.rs
│   ├── window.rs          # Window creation, resize handling
│   ├── input.rs           # Keyboard/mouse event handling
│   └── app.rs             # Main loop, expose/redraw events
```

**Windowing crate comparison**:

| Crate | Pros | Cons |
|-------|------|------|
| `winit` | Cross-platform, simple, raw-window-handle compatible | Raw, need backend |
| `macroquad` | Simple game loop | Opinionated |
| `ggez` | Game-focused | May be too heavy |

**Recommendation**: Start with `winit` + `raw-window-handle` for flexibility.

---

## Phase 4: GPU Acceleration (Optional Enhancement)

**Goal**: Modernize rendering with GPU.

```
dgraphx-gpu/
├── src/
│   ├── lib.rs
│   ├── shader.rs          # GLSL/WGSL vertex + fragment shaders
│   └── pipeline.rs        # Render pipeline setup
```

**GPU crate comparison**:

| Crate | Pros | Cons |
|-------|------|------|
| `wgpu` | Cross-platform, safe, WebGPU standard | Verbose |
| `vulkano` | Idiomatic Rust | Smaller ecosystem |
| `miniquad` | Simple, mobile support | Less flexible |

**Recommendation**: `wgpu` for cross-platform GPU rendering.

---

## Phase 5: Integration & Features

**Goal**: Build the complete application from core components.

```
dgraphx/
├── Cargo.toml (workspace)
├── crates/
│   ├── dgraphx-math/       # Phase 1
│   ├── dgraphx-renderer/   # Phase 2
│   ├── dgraphx-core/       # Phase 3
│   └── dgraphx-gpu/        # Phase 4 (optional)
├── examples/               # Demos for each plot type
├── benches/                # Benchmarks
└── src/
    ├── lib.rs
    ├── graph2d.rs          # 2D plotting
    ├── graph3d.rs           # 3D surface plots
    ├── parametric.rs        # Parametric curves
    └── ui.rs               # GUI controls
```

---

## Migration Timeline

| Phase | Description | Estimated Time |
|-------|-------------|---------------|
| Phase 1 | Math Library | 4-6 weeks |
| Phase 2 | Software Renderer | 3-4 weeks |
| Phase 3 | Windowing & Event Loop | 2-3 weeks |
| Phase 4 | GPU Acceleration | 4-6 weeks (optional) |
| Phase 5 | Integration & Features | 4-6 weeks |

**Total**: 17-25 weeks (implementation only)

---

## Feature Enhancements

### Planned Additions

1. **Function registration** - Parse function strings at runtime instead of recompiling
2. **Multiple plot views** - Split screen with multiple graphs
3. **Animation** - Time-based function evaluation for animated plots
4. **Export** - PNG, SVG, PDF output
5. **Mouse support** - Pan, zoom, rotate with pointer
6. **Scripting** - Lua/Rhai integration for custom functions

### Improvements Over Original

| Aspect | C Version | Rust Version |
|--------|-----------|--------------|
| Testing | None | Unit + integration tests |
| Memory safety | Manual | Compile-time guarantees |
| Error handling | printf/return codes | `Result` types |
| Cross-platform | X11 only | Windows, macOS, Linux, Web |
| GPU rendering | None | Optional wgpu backend |
| Function definition | Recompile | Runtime parsing |
| Documentation | Basic README | Full API docs |

---

## Development Strategy

### Branching Model

```
c/branch          # Original C version (archived, reference)
main or rust/     # New Rust implementation
```

### Validation Approach

1. Keep C version in `c/` branch for reference
2. Develop Rust in `main` (or `rust` branch)
3. **Pixel-diff validation**: Compare Rust output pixel-by-pixel with C version during development
4. Run benchmarks to validate performance parity
5. Property-based testing with `proptest` for math functions

### Dependencies Philosophy

- **Math core**: Zero external dependencies (pure Rust)
- **Rendering**: `image` crate for pixel buffers
- **Windowing**: `winit` for cross-platform window management
- **GPU**: `wgpu` for GPU rendering
- **Serialization**: `serde` for config/file formats
- **Testing**: `proptest`, `criterion` for benchmarks

---

## C Algorithm Reference

### Matrix Operations (matrix.c)

```c
void MxV(double* V, double* v);      // Matrix-vector multiply
void MxM(double* V, double* N);      // Matrix-matrix multiply
int InvMatrix(int n, double* A, double* B);  // Matrix inverse
void cross(double* a, double* b, double *c);  // Cross product
void normalize(double *a, int n);    // Vector normalization
void MatrixId(double *N);            // Identity matrix
double inner(double *v, double* u, int n);  // Dot product
```

### Interpolation (interpolation.c)

```c
void interpH(double* x, int n, double *m);  // Hermite interpolation
void interpS(double* x, int n, double *m); // Spline interpolation
```

### 3D Transformations (DX.c)

```c
void RotateX(int mode);  // Rotation matrices
void RotateY(int mode);
void RotateZ(int mode);
void IdMatrix();
void ChangeCameraPosition(double valueOfShift);
void ChangeScreenPosition(double valueOfShift);
void ObliqueProjection(double x, double y, double z);
```

### Solvers (Interp2D.h / DGraphX.c)

```c
void SeidelGauss(double *a, int n, int *niter);
void RichardsonIter(double *a, int n, double min, double max, int *niter);
```

---

## File Structure After Migration

```
DGraphX/
├── RUST-MIGRATION.md      # This file
├── AGENTS.md              # Agent documentation
├── c/                      # Original C version (archived)
│   ├── DGraphX.c
│   ├── DGraphX.h
│   ├── Xfun.c
│   ├── matrix.c / matrix.h
│   ├── interpolation.c
│   ├── dtoa.c
│   ├── DX.c
│   ├── Interp2D.h
│   ├── Makefile
│   └── data.txt
├── crates/                 # Rust workspace
│   ├── dgraphx-math/
│   ├── dgraphx-renderer/
│   ├── dgraphx-core/
│   └── dgraphx-gpu/       # Optional
├── examples/
├── benches/
└── src/                   # Main application
    ├── main.rs
    ├── graph2d.rs
    ├── graph3d.rs
    └── ui.rs
```

---

## Initial Implementation Blueprint

A sample structure for the core math port:

```rust
// src/math/matrix.rs
pub struct Matrix4 {
    pub data: [f64; 16],
}

impl Matrix4 {
    pub fn multiply_vec(&self, v: [f64; 3]) -> [f64; 3] {
        // Ported logic from matrix.c
    }
}

// src/renderer/canvas.rs
pub struct Canvas {
    pub pixels: Vec<u32>,
    pub width: usize,
    pub height: usize,
}

impl Canvas {
    pub fn draw_line(&mut self, x0: i32, y0: i32, x1: i32, y1: i32, color: u32) {
        // Implementation of line drawing
    }
}
```

## Next Steps

1.  **Pilot:** Port the `matrix.c` logic to a standalone Rust crate and verify results against the C version.
2.  **Prototype:** Build a minimal window that renders a rotating 3D wireframe cube using the ported math.
