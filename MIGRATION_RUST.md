# DGraphX Rust Migration Plan

This document outlines the strategy for migrating the **DGraphX** project from C/X11 to a modern Rust-based architecture. The goal is to preserve the "from scratch" rendering philosophy while gaining memory safety, high-performance parallelism, and cross-platform support.

## Core Philosophy
*   **Software Rendering:** Avoid heavy GPU-bound engines. Keep the logic for rasterization, 3D projection, and pixel manipulation in the application code.
*   **Zero-Dependency Logic:** The core math and rendering algorithms should remain independent of external math libraries where possible
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
| **Scripting** | Hardcoded C functions | **evalexpr** | Allow users to input functions at runtime. |

---

## Migration Phases

### Phase 1: Foundation & Windowing
*   Initialize a new Cargo project.
*   Set up `winit` for the event loop and window management.
*   Integrate `softbuffer` to create a raw pixel buffer (Canvas).
*   **Goal:** A window that displays a solid color or a simple test pattern (e.g., a static rectangle).

### Phase 2: Math & Transformation Logic
*   **Port Matrices:** Port the 4x4 matrix multiplication and transformation logic from `matrix.c` and `DX.c`. Use nalgebra for matrix ops
*   **Projection Pipeline:** Create a `Camera` struct that handles:
    *   Rotation (Euler angles or Quaternions).
    *   Projection (Perspective or Orthographic).
    *   Viewport mapping (3D coordinates to pixel coordinates).
*   **Interpolation:** Port the Hermite and Spline algorithms to a `math::interpolation` module.

### Phase 3: The "From Scratch" Renderer
*   Implement a `Canvas` trait or struct with methods for:
    *   `set_pixel(x, y, color)`
    *   `draw_line(p1, p2, color)` (Bresenham's or similar).
    *   `fill_polygon(vertices, color)` (for 3D faces).
*   Port the Z-buffering or painter's algorithm logic to handle depth in 3D views.

### Phase 4: Application State & Interaction
*   Model the application state using Enums:
    ```rust
    enum AppMode {
        Graph2D(Variant2D),
        Graph3D(Variant3D),
    }
    ```
*   Implement keyboard handling to replicate current DGraphX controls (Rotation, Zoom, Mode switching).
*   Replace `data.txt` parsing with a structured format like JSON or TOML using `serde`.
*   Add mouse support

### Phase 5: Modern Enhancements
*   **Parallel Plotting:** Use `rayon` to calculate the heightmap of 3D surfaces across multiple CPU cores.
*   **Live Evaluation:** Integrate a math expression parser so users can type functions in the terminal or a GUI text box.
*   **WebAssembly (WASM):** Compile the project for the web so it can run in any modern browser via a `<canvas>`.

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
