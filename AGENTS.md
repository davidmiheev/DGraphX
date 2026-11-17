# AGENTS.md - DGraphX Codebase Guide

## Project Overview

DGraphX is a lightweight 2D/3D graphing application written in C, using X11 for graphics rendering. It requires no external dependencies beyond X11 (Xlib), pthreads and a C compiler. 3D plots are rendered by a software z-buffer rasterizer (no GPU).

## Build Commands

```sh
make          # Build the DGraphX executable
make clean    # Remove object files and executable
```

- **Compiler**: GCC (configured in Makefile)
- **Linker**: GCC with `-lX11 -lm -lpthread`
- **X11 includes**: `/opt/X11/include` (X11 lib at `/opt/X11/lib`)
- **Optimization**: `-O2` by default (adjustable via `OPTIMAZE_COMMON`)

## Project Structure

```
DGraphX/
├── DGraphX.c          # Main entry point, drawing functions, keyboard handling
├── DGraphX.h          # Main header - global state, types, drawing API
├── Xfun.c             # X11 window, event loop, X11 primitives, software framebuffer (FB_*)
├── matrix.c/h         # Matrix operations (MxV, MxM, cross, normalize, etc.)
├── interpolation.c    # 2D interpolation (Hermite, spline)
├── dtoa.c             # Number-to-string conversion utility
├── tinyexpr.c/h       # Vendored expression parser (Zlib) - runtime function input
├── DX.c               # 2D mapping, 3D camera/rotation, tessellation, shading
├── Interp2D.h         # Interpolation function declarations
├── data.txt           # Polyhedra vertex data (read at runtime)
├── Makefile           # Build configuration
├── AGENTS.md          # This guide
├── RUST-MIGRATION.md  # Notes/plan for a possible Rust rewrite
└── .github/workflows/c-cpp.yml  # CI - installs libx11-dev, runs `make` on push/PR to master
```

## Code Patterns

### Global State
All global variables declared in `DGraphX.h` and defined in `Xfun.c`:
```c
extern COLORS rColors;
extern int width, height;
extern double xc, yc, xmin, xmax, ymin, ymax;
extern int x_0, y_0;
```

### Drawing Architecture
- `DrawWindow()` takes three callbacks:
  - `DrawWindowContent` - renders the current frame
  - `KeyPressFunction` - handles keyboard input, returns `KEY_PRESS_QUIT`, `KEY_PRESS_EXPOSE`, or `KEY_PRESS_NOTHING`
  - `MouseEventFunc` - handles ButtonPress/ButtonRelease/MotionNotify, same return contract
- Drawing happens to an internal `draw` pixmap (3000x3000), copied to the window with `XCopyArea`
- 2D mode draws directly with Xlib primitives (`WDrawLine`, `DrawLine`, ...); 3D mode goes
  through the software framebuffer described below

### Runtime function expressions (tinyexpr)
- Plottable functions are tinyexpr strings compiled by `compile_functions()` in `DGraphX.c`
- Press `i` in the window to enter a new expression via stdin (blocks the event loop; needs a tty)
- Evaluation goes through globals `cur_x/cur_y/cur_t/cur_s` - the wrappers are NOT reentrant
- `df()` is a central finite difference of `f()`, used by Hermite/spline interpolation

### 3D rendering pipeline (DX.c -> Xfun.c)
- **Camera**: one 3x3 rotation `M` plus a perspective divide. Camera space is x = screen right,
  **y = depth** (`w = cameraPos - y`), z = screen vertical. `cameraPos`/`screenPos` set distance
  and zoom (`ChangeCameraPosition`, `ChangeScreenPosition`).
- **Rotation**: keyboard rotation (`RotateX/Y/Z`, `elemRotation`) right-multiplies `M = M*N`,
  i.e. rotates about **world** axes. Mouse drag (`RotateDrag` -> `camRotation`) left-multiplies
  `M = N*M`, i.e. rotates about **screen** axes (trackball): vertical drag about camera X,
  horizontal drag about camera Z. `M` is never re-orthonormalized.
- **Z-buffer, not painter's algorithm**: each draw function (`DrawGraph3DX`, `ParametricCurve3D`,
  `ParametricGraph3D`, `DrawPolytope`, `DrawLinear3D`) calls `FB_Begin(bg)`, emits primitives,
  then `FB_Flush()`. There is no depth sorting.
  - `rasterTriA()` clips triangles against the near plane in camera space, projects them and
    records `FB_Tri` (flat colour) or `FB_TriA` (per-vertex float attributes + a shade callback).
    Axis/grid/contour lines are recorded with `FB_Seg` (axis segments come from `zm[]`/`axq[]`).
  - Depth is `q = 1/w` per vertex (interpolated affinely in screen space, larger = nearer).
    Attributes are interpolated **perspective-correctly** (`a/w` divided by interpolated `1/w`).
  - `FB_Flush` replays the recorded commands on up to 8 threads, each owning a band of rows,
    into `fbPix` (colour) / `fbZ` (depth), then `XPutImage`s into the pixmap. Shade callbacks run
    on worker threads: they must not touch shared mutable state (use locals; `ShadeCtx` is
    copied into each command).
- **Shading**: `pallette()` precomputes a 64x32 colour table indexed by diffuse (`n.light`) and
  specular (`r.view`) terms. `shadePhong` = per-pixel normals (default, `P` toggles),
  `flatShade` = per-face, `shadeHeight` = colour by height (`L` toggles colour mode).
- **Overlays**: floor grid (`DrawFloorGrid`, `G` toggles `gridEnabled`); variant 1 also draws
  marching-squares contours on the floor and iso-curves on the surface.
- **Caching**: tessellation is cached per `renderEpoch` (bumped only on geometry-affecting
  changes: `i` edits, `+`/`-`, `Tab`, `data.txt` reload). Cached: `pgv` (parametric surface
  verts), `gv`/`gm` (z=f(x,y) grids), `rings` (tube cross-sections). Camera moves only
  re-project and re-rasterize.

### Event loop (Xfun.c)
- Queued X events are processed in batches (cap 64) and the frame is redrawn once per batch -
  don't add per-event drawing calls. Motion events are coalesced.
- Redraws are capped at ~30 fps; a skipped redraw is carried over (`pendingRedraw`).
- Progressive refinement (`R` toggles `refineEnabled`): frames triggered by input render with
  `interactive = 1` (draw functions tessellate at stride 2). After any batch, if a coarse frame
  is on screen and no input arrives for 40 ms, a full-resolution frame replaces it.
- `Backspace` (`key_6`) hides the current plot in both 2D and 3D.

### Color Macros
Colors accessed via global `rColors` struct:
```c
WSetColor(BLUE);   // not a direct value - macro expands to rColors.blue
```

### Key Code Patterns
- Keyboard uses X11 `KeySym` (e.g., `XK_Q`, `XK_Left`, `XK_Return`)
- State tracked via static variables in `DGraphX.c` (e.g., `static int variant`, `static int key_1`)
- 2D/3D mode toggle via `Ctrl` key, resets state

## Header Dependencies

```
DGraphX.h (main)
├── X11/keysym.h, X11/Xlib.h
├── Interp2D.h
└── matrix.h

Xfun.c
├── DGraphX.h
└── Interp2D.h

DGraphX.c
├── DGraphX.h
├── Interp2D.h
├── matrix.h
└── tinyexpr.h
```

## Object File Dependencies

```makefile
Xfun.o       : Xfun.c DGraphX.h Interp2D.h
DGraphX.o     : DGraphX.c DGraphX.h Interp2D.h matrix.h
interpolation.o: interpolation.c Interp2D.h DGraphX.h
dtoa.o        : dtoa.c Interp2D.h DGraphX.h
DX.o          : DX.c DGraphX.h matrix.h
matrix.o      : matrix.c DGraphX.h matrix.h
tinyexpr.o    : tinyexpr.c tinyexpr.h
```

## Known Issues / Gotchas

1. **No test suite**: No automated tests exist. Verify manually with `./DGraphX`.
2. **Runtime data dependency**: `data.txt` must exist for polyhedra (variant 4) to load.
3. **Makefile typo**: Variable is named `OPTIMAZE_*` (not `OPTIMIZE_*`).
4. **Unicode in printf**: Some printf strings contain Cyrillic characters in comments.
5. **X11 display required**: Program exits with status `X11_ERR_1` (1) if no X server is available (use `Xvfb` for headless testing).
6. **Fonts**: Adobe Helvetica (`HELVETICA16/12`) is often not installed; `SetFont` falls back to
   another Adobe font or the core `fixed` font. Fonts are cached per name and freed on quit.
7. **Solver demo (variants 5-7)**: `func()` in `interpolation.c` is built from the exact solution
   `sin(2πx)sin(2πy)/2`, so the displayed error is the true max error (O(h²)). If you change
   `func`, change `solution()` in `DGraphX.c` to match or the error becomes meaningless.
   `max_iter` (`i` -> `iters`) caps individual sweeps for both Seidel-Gauss and Richardson.

## CI Configuration

GitHub Actions runs on every push to `master` and on pull requests. It installs `libx11-dev` and executes `make` on `ubuntu-latest`. It only checks that the code compiles.

## Style Notes

- Hungarian notation: `rColors`, `prDisplay`, `nWnd`, `sString`
- Function naming: `CamelCase` for X11 wrappers (`WDrawString`, `WSetColor`), `snake_case` for math utilities
- Comments: C-style `//` and `/* */`
- No clang-format or linter configured

## Improvement Opportunities

### High Priority

1. **Add automated tests** - No test suite exists. Even basic unit tests for matrix operations, interpolation, and dtoa would prevent regressions.

2. **Extract magic numbers into constants** - `NMAX`, `STEP`, `LAMBDA`, `EPS`, etc. are hardcoded in `DGraphX.c`. Move to a config/params header.

3. **Refactor XK_i handler** - The 10-way strcmp chain for named expressions is duplicated for update/revert; a `{name, expr, target}` table would collapse it.

### Medium Priority

4. **Refactor DGraphX.c** - This file mixes math function definitions, drawing logic, keyboard handling, and main(). Could be split into separate modules.

5. **Fix Makefile typos** - `OPTIMAZE_*` variables should be `OPTIMIZE_*`.

6. **Add stricter warning flags** - Currently `-W -Wall`. Consider adding `-Wextra`, `-Wpedantic`.

7. **Add sanitizers for development** - ASan, UBSan, MSan configs for debugging builds.

### Lower Priority

8. **Improve error handling** - Some X11 errors aren't caught (no `XSetErrorHandler`); most diagnostics still go to stdout.

9. **Re-orthonormalize the view matrix** - `M` accumulates rounding error over many drags.

10. **Add clang-format config** - No code formatting standard enforced.

11. **Consider XRM (X11 Resource Manager)** - For theme/font preferences instead of hardcoded Adobe Helvetica font names.

### From the old GEMINI.md roadmap

12. **CLI arguments** - Flags for initial mode (2D/3D), data file path, camera settings.

13. **Makefile portability** - Replace hardcoded `/opt/X11` paths with `pkg-config` for cross-distro builds.

14. **Context struct for global state** - Encapsulate the `extern` globals into a struct to improve testability.

15. **Headless smoke tests** - Use `Xvfb` in CI to run basic startup/render checks.
