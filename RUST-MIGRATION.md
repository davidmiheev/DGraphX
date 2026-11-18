# DGraphX Rust Migration Plan

Status: **proposal**. Checked against the C code after the z-buffer renderer, runtime
expressions and trackball rotation landed (PR #2). Update the "Current C architecture" and
"Feature parity checklist" sections whenever the C version changes before the port is done.

## Goals and non-goals

**Goals**
- Reach **feature parity** with the C version first, then add features.
- Memory and thread safety; typed state instead of global flags.
- Testability: the math, plotting and rendering crates run headless.
- Cross-platform: Linux (X11 and Wayland), macOS, Windows. Web is a stretch goal.
- Keep the project's character: **software rendering**, small dependency footprint, runs on
  slow machines without a GPU.

**Non-goals**
- GPU rendering (`wgpu` etc.). It contradicts the software-rendering principle and doubles
  the renderer surface. Revisit only after parity, with a separate decision record, if the
  software renderer cannot meet the performance targets below.
- Game engines or GUI frameworks that own the main loop (`macroquad`, `ggez`, ...).

## Why Rust: concrete wins for this codebase

- **Parallel tessellation.** The C expression wrappers evaluate through shared globals
  (`cur_x/cur_y/cur_t/cur_s`) and are not reentrant, so tessellation is single-threaded. In
  Rust a compiled expression is an immutable `Fn(f64, f64) -> f64 + Send + Sync`, so surface
  and tube tessellation can run on `rayon`.
- **Race-free rasterization.** The C rasterizer gives each thread a band of rows by index
  arithmetic. In Rust each worker receives a disjoint `&mut` slice of the colour and depth
  buffers (`par_chunks_mut`), so the compiler rules out data races, including the shared
  scratch-state race fixed during review of PR #2.
- **Typed state.** `key_1` ... `key_12`, `variant`, `dim`, `mode`, `count` become an
  `AppState` struct plus a `Command` enum and a key-binding table.
- **Errors as values.** X11 failures, `data.txt` parsing and expression errors become
  `Result`s with messages, instead of mixed `printf` and silent fall-backs.

## Current C architecture (what must be ported)

Sizes: `DX.c` ~1,350 lines, `DGraphX.c` ~910, `Xfun.c` ~700, `interpolation.c` ~280,
`matrix.c` ~190. **Most of the work is in `DX.c` and the UI behaviour in `DGraphX.c`**, not in
the matrix code.

| Area | C location | Notes for the port |
|---|---|---|
| Camera | `DX.c`: `M`, `Projection`, `ProjectW`, `InitCameraPosition` | 3x3 rotation `M` + perspective divide `(cameraPos - screenPos) / w`. Camera space: x = right, **y = depth** (`w = cameraPos - y`), z = vertical. No 4x4 matrices are needed. |
| Rotation | `elemRotation` (`M = M*N`), `camRotation`/`RotateDrag` (`M = N*M`) | Keyboard rotates about **world** axes, mouse drag about **screen** axes (trackball). `M` is never re-orthonormalized in C: fix in the port. |
| Near-plane clipping | `DX.c`: `rasterTriA` | Clip in camera space against `w = cameraPos - screenPos + 0.1`, lerp attributes, fan-triangulate. |
| Rasterizer | `Xfun.c`: `FB_Begin`, `FB_Tri`, `FB_TriA`, `FB_Seg`, `FB_Flush` | Commands recorded per frame, replayed per row band on up to 8 threads. Bounding-box + edge-function (barycentric) rasterization, not scanline. Depth = `1/w`, larger is nearer. Attributes are **perspective-correct** (`a/w` over `1/w`). |
| Shading | `pallette`, `shadePhong`, `flatShade`, `shadeHeight` | 64x32 colour table indexed by diffuse (`n.light`) and specular (`r.view`). Per-pixel normals and view vectors (Phong), per-face (flat), or colour by height. |
| Tessellation and caches | `DrawGraph3DX` (`gv`/`gm`), `ParametricGraph3D` (`pgv`), `ParametricCurve3D` (`rings`), `DrawPolytope`, `DrawLinear3D` | Cached per `renderEpoch`; camera moves only re-project. Stride-2 tessellation while interacting. |
| Tube along a curve | `tubeTangent`, `transportFrame`, `tubeRingAt` | Parallel-transport frames (no twist), one-sided differences at domain edges, NaN spans leave gaps instead of erasing the curve, backface culling. |
| Overlays | `DrawAxes`, `DrawFloorGrid`, `DrawFloorContours`, `DrawSurfaceCurves` | Floor grid; marching-squares contours on the floor and iso-curves on the surface (variant 1). |
| 2D plotting | `DX.c`: `map`, `invmap`, `scale`, `xshift`, `yshift`, `DrawGrid`, `DrawAxes2D`, `DrawGraph2D`, `DrawParametric2D` | Drawn with Xlib lines directly (not through the framebuffer). |
| Interpolation | `interpolation.c`: `interpH`, `interpS`; node placement in `DrawWindowContent` | Hermite and cubic spline; nodes uniform, random insert/remove, or Chebyshev. `df()` is a finite difference of the runtime `f()`. |
| Numerics | `solveEq` (BVP by shooting + bisection, `shoot`), `SeidelGauss`, `RichardsonIter` (5 Chebyshev parameters) | Variants 5-7 use a manufactured solution (`func` built from `solution`), so the displayed error is the true max error, O(h^2). `max_iter` caps sweeps for both iterative solvers. |
| Expressions | `tinyexpr.c`, `compile_functions`, `XK_i` handler | 10 named slots (`f, u, v, w, cx, cy, cz, px, py, pz`) + `iters`. Invalid input keeps the previous expression. Undefined points evaluate to NaN and are skipped. |
| Event loop | `Xfun.c`: `DrawWindow` | Batches up to 64 events, coalesces motion, caps redraws at ~30 fps (skipped redraws carried over), coarse frame during input then a full frame after 40 ms of quiet. |
| Fonts | `SetFont` | X core fonts, cached per name, fallback to `fixed`. |

**Do not port** (dead or replaced by the standard library): `ObliqueProjection`, `PrintMatrix`,
`comp_1`, `LineSightControlAxes`, `DrawVector`, `ufun`, `dtoa.c`/`itoa`/`scat` (use `format!`),
the general `InvMatrix` (it only inverts rotations, whose inverse is the transpose), the
commented-out older `solveEq`, and the never-reached stdin prompts in `main()` (`dim` starts
at 3, so they never run; replace with CLI flags).

## Feature parity checklist

Every row needs a Rust acceptance check (unit, golden or manual) before the C version is
retired.

**Modes and plots**
- [ ] 3D title screen (variant 0) with help text; `Ctrl` toggles 2D/3D and resets state.
- [ ] 3D variant 1: surface `z = u(x,y)` on [-5,5]^2. (`v` and `w` are editable but **not drawn**:
      `DrawGraph3DX` is called with `n = 1`. Decide whether to keep that.)
- [ ] 3D variant 2: curve `(cx, cy, cz)(t)`, t in [0, 10π], drawn as a lit tube.
- [ ] 3D variant 3: parametric surface `(px, py, pz)(t, s)` on [-π,π]^2.
- [ ] 3D variant 4: polyhedron from `data.txt` (faces <= 100, vertices per face <= 10, bad
      file falls back to the title screen with a message).
- [ ] 3D variants 5/6/7: exact solution / Seidel-Gauss / Richardson on an `ngrid` x `ngrid`
      grid with grid size, error, iteration count and time labels.
- [ ] 2D variants 1-5: Hermite interpolation, spline interpolation, graph of `f`,
      parametric curve `(cx, cy)(t)`, BVP solution with the equation label. (Known C bug:
      equation 1 is labelled `exp(x)` but `shoot` solves with `100 sin(20x)`; fix, don't copy.
      `solveEq` also writes `res.txt` and prints to stdout every frame; drop that.)

**Keyboard**
- [ ] `q` quit, `0` title, `1`-`7` / `F1`-`F7` select variant.
- [ ] 3D: `←`/`→` rotate about the selected world axis (`X`/`Y`/`Z` select), or move the
      light in light mode (`l`, only in shaded colour mode); `↑`/`↓` change camera distance
      (`s`) or screen distance (`f`); `C`/`V` zoom; `Enter` axes; `L` colour mode
      (shaded / by height); `P` Phong; `G` floor grid; `R` progressive refinement.
- [ ] 2D: arrows pan, `C`/`V` zoom, `n`/`b` add/remove nodes, `8`/`9` uniform/Chebyshev
      nodes, `Enter` node markers.
- [ ] Both: `Backspace` hide/show plot, `Tab` cycle coefficient (3D: 1 / 0.1 / 0.01;
      2D variant 5: equation), `+`/`-` solver grid size (10..100), `i` edit an expression.

**Mouse**
- [ ] 3D: left-drag trackball rotation, right-drag camera distance, wheel zoom.
- [ ] 2D: left-drag pan, right-drag and wheel zoom about the view centre.

**Behaviour**
- [ ] Coarse frame while interacting, full frame after input stops; never left coarse.
- [ ] Camera moves do not re-evaluate expressions (tessellation cache).
- [ ] NaN/undefined regions leave gaps; nothing is drawn at a fabricated 0.
- [ ] Non-zero exit status when no window can be created.

## Target architecture

```
dgraphx/
├── Cargo.toml              # workspace (edition 2024, Rust >= 1.85)
├── c/                      # C reference, built and run in CI until parity
├── crates/
│   ├── dgraphx-math/       # camera, rotations, interpolation, solvers (glam f64 types)
│   ├── dgraphx-expr/       # expression slots -> Fn + Send + Sync, error positions
│   ├── dgraphx-render/     # Canvas (u32 colour + f32 depth), clip, raster, shading LUT,
│   │                       # lines, text; PNG export behind a feature
│   ├── dgraphx-plot/       # tessellation per plot type + caches; pure, no window
│   └── dgraphx-app/        # winit + softbuffer shell, AppState, Command, keymap, loop
├── tests/golden/           # numeric and image references produced from the C build
└── benches/
```

Core loop:

```rust
enum Command {
    SelectVariant(u8), ToggleDimension, RotateWorld { axis: Axis, dir: i8 },
    Drag { button: Button, dx: f64, dy: f64 }, Zoom(f64), Toggle(Flag),
    EditExpr { slot: Slot, source: String }, Quit, /* ... */
}

struct AppState { mode: Mode, variant: Variant, camera: Camera, flags: Flags,
                  exprs: Exprs, epoch: u64, /* ... */ }

fn update(state: &mut AppState, cmd: Command) -> Redraw;          // pure, unit-testable
fn render(state: &AppState, quality: Quality, canvas: &mut Canvas); // headless-testable
```

- Input events are mapped to `Command`s by a key-binding table; the shell coalesces motion
  events, applies the redraw cap and drives `Quality::Coarse` / `Quality::Full`.
- **Camera**: store a rotation (`DMat3` re-orthonormalized after each update, or `DQuat`).
  Offer **trackball** (current behaviour) and **turntable** (horizontal drag about world Z,
  keeps the floor level) modes.
- **Expression input** moves into the window (a one-line text field drawn on the canvas).
  Reading stdin freezes the C window and does not exist on Windows GUI apps or the web.
- **Text**: an embedded bitmap font (for example `font8x8`) or `fontdue` with a bundled TTF.
  Check the font licence.
- **Precision**: f64 for geometry and math, f32 for depth and interpolated attributes (as in C).

## Dependencies

| Need | Choice | Notes |
|---|---|---|
| Window + input | `winit` | Evaluate `minifb` too: one crate for window + pixel buffer, simpler, less flexible. |
| Presenting pixels | `softbuffer` | CPU buffer to window, no GPU. Handle HiDPI scale factors. |
| Vector/matrix math | `glam` (`DVec3`, `DMat3`, `DQuat`) | Small and well tested; replaces "zero math dependencies", which would mean re-implementing and testing the same types. |
| Parallelism | `rayon` | Feature-gated: wasm needs `wasm-bindgen-rayon` plus cross-origin isolation, so keep a single-threaded fallback. |
| Expressions | **decide by benchmark** | Candidates: `exmex` (can compute partial derivatives, giving exact normals and `df` instead of finite differences), `evalexpr` (dynamic typing, likely too slow for hot loops), or a small hand-written parser compiling to closures (what tinyexpr does). Criteria: speed on ~330k evaluations per surface rebuild, `Send + Sync`, error positions, maintenance status. |
| PNG export / tests | `image` | Behind an `export` feature and in dev-dependencies only. |
| Testing | `proptest`, `criterion`, `cargo-fuzz` | Property tests, benchmarks, parser fuzzing. |
| Config / CLI | `serde`, `clap` | Phase 5. |

## Validation strategy

Pixel-for-pixel equality with C is not realistic (floating-point differences, fonts, line
rasterization). Compare data exactly and images with tolerance:

1. **Numeric goldens from C**: tessellation vertices for each plot type, interpolation
   coefficients, solver grids and iteration counts. Compare with a relative tolerance (~1e-9).
2. **Image goldens**: the C build renders each variant under `Xvfb` in CI; the Rust build
   renders the same scene headless. Pass if, for example, <= 1% of pixels differ by more than
   8 levels, with text regions masked.
3. **Property tests**: rotations stay orthonormal after many random drags; perspective-correct
   interpolation reproduces affine 3D attributes exactly; clipped triangles never reach
   `w < near`.
4. **Convergence tests**: solver error falls ~4x when `ngrid - 1` doubles (O(h^2)), as
   verified on the C version.
5. **Fuzzing** of the expression parser.
6. **Benchmarks**: frame time per variant at 800x700, C vs Rust, coarse and full quality.
   Initial target: >= 30 fps full quality for variants 1 and 3 on 4 cores.

## Phases (vertical slices)

Each phase ends with something runnable and an exit criterion. Effort is relative
(S < M < L); calibrate in weeks after Phase 0.

**Phase 0 - Prepare the C reference (S, in C).**
Turn the checklist above into tracked issues. Add a headless dump to C (for example
`DGRAPHX_DUMP=<dir>` writes the framebuffer as PPM and the tessellation/solver data as text),
an `Xvfb` smoke test in CI, and baseline benchmarks.
*Exit*: CI publishes numeric and image goldens.

**Phase 1 - Foundations (S-M).**
Workspace, `dgraphx-math` (camera, interpolation, solvers), `dgraphx-expr` after the
benchmark decision.
*Exit*: numeric goldens for interpolation and solvers pass.

**Phase 2 - Renderer (M).**
Canvas, near-plane clipping, triangle and segment raster, perspective-correct attributes,
shading table, banded `rayon` flush, PNG export.
*Exit*: image goldens for fixed test meshes pass.

**Phase 3 - First vertical slice (M).**
`dgraphx-app` shell, `AppState`/`Command`/keymap, trackball and turntable camera, event
coalescing, progressive refinement, 3D variant 1 interactive with in-window expression editing.
*Exit*: the variant-1, keyboard, mouse and behaviour rows of the checklist pass.

**Phase 4 - Parity (L).**
3D variants 2-7, 2D mode (graphs, interpolation, parametric, BVP), overlays, text, `data.txt`.
*Exit*: the whole checklist passes; tag the last C commit (`c-final`) and remove `c/`.

**Phase 5 - New features (open-ended, only after parity).**
Export (PNG, then SVG/PDF), multiple views, animation (time variable), CLI and config file,
web build.

## Risks

- Expression evaluation speed (mitigated by the Phase 1 benchmark).
- Text rendering and HiDPI differences across platforms.
- Wayland: no global window positioning; scale factors change buffer sizes.
- Web: threads need cross-origin isolation; no filesystem for `data.txt` (bundle or upload).
- Scope creep: new features before parity. Hold them until Phase 5.

## Open decisions

- Expression crate (Phase 1 benchmark).
- `winit` + `softbuffer` vs `minifb`.
- Default camera mode: trackball (current) or turntable.
- Whether variant 1 should draw `v` and `w` as well as `u`.
- Replacement for the unreachable 2D startup prompts (CLI flags vs in-app settings).
