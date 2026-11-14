# AGENTS.md - DGraphX Codebase Guide

## Project Overview

DGraphX is a lightweight 2D/3D graphing application written in C, using X11 for graphics rendering. It requires no external dependencies beyond X11 and a C compiler.

## Build Commands

```sh
make          # Build the DGraphX executable
make clean    # Remove object files and executable
```

- **Compiler**: GCC (configured in Makefile)
- **Linker**: GCC with `-lm` for math library
- **X11 includes**: `/opt/X11/include` (X11 lib at `/opt/X11/lib`)
- **Optimization**: `-O2` by default (adjustable via `OPTIMAZE_COMMON`)

## Project Structure

```
DGraphX/
├── DGraphX.c          # Main entry point, drawing functions, keyboard handling
├── DGraphX.h          # Main header - global state, types, drawing API
├── Xfun.c             # X11 window management, event loop, X11 primitives
├── matrix.c/h         # Matrix operations (MxV, MxM, cross, normalize, etc.)
├── interpolation.c    # 2D interpolation (Hermite, spline)
├── dtoa.c             # Number-to-string conversion utility
├── tinyexpr.c/h       # Vendored expression parser (Zlib) - runtime function input
├── DX.c               # 3D transformations and projections
├── Interp2D.h         # Interpolation function declarations
├── data.txt           # Polyhedra vertex data (read at runtime)
├── Makefile           # Build configuration
└── .github/workflows/c-cpp.yml  # CI - runs `make` on push/PR to master
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
- Drawing happens to an internal `draw` pixmap (3000x3000), copied to window on expose

### Runtime function expressions (tinyexpr)
- Plottable functions are tinyexpr strings compiled by `compile_functions()` in `DGraphX.c`
- Press `i` in the window to enter a new expression via stdin (blocks the event loop; needs a tty)
- Evaluation goes through globals `cur_x/cur_y/cur_t/cur_s` - the wrappers are NOT reentrant
- `df()` is a central finite difference of `f()`, used by Hermite/spline interpolation

### Rendering pipeline (DX.c)
- Painter's algorithm: primitives (axis segments in `z[]`/`zm[]`, triangles in `zb[]`) are
  depth-sorted by `inner(view, centroid)` and drawn far-to-near. `radixSortPairs()` does the
  sort (quantized uint32 keys, in-place on `Pair.z`); `polygon` entries are sorted indirectly
  via `z[]` pairs with `num = ~index` (negative marks polys, non-negative is a `zm[]` line).
- Tessellation is cached per `renderEpoch` (global, bumped only on geometry-affecting keys:
  `i` edits, `+`/`-`, `Tab`). Cached: `pgv` (parametric surface verts), `gv`/`gm` (z=f(x,y)
  grids), `rings` (tube cross-sections). Camera moves only re-project + re-sort.
- The event loop in `Xfun.c` processes queued X events in batches (cap 64) and redraws once
  per batch — don't add per-event drawing calls.

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
└── matrix.h
```

## Object File Dependencies

```makefile
Xfun.o       : Xfun.c DGraphX.h Interp2D.h
DGraphX.o     : DGraphX.c DGraphX.h Interp2D.h matrix.h
interpolation.o: interpolation.c Interp2D.h DGraphX.h
dtoa.o        : dtoa.c Interp2D.h DGraphX.h
DX.o          : DX.c DGraphX.h matrix.h
matrix.o      : matrix.c DGraphX.h matrix.h
```

## Known Issues / Gotchas

1. **No test suite**: No automated tests exist. Verify manually with `./DGraphX`.
2. **Runtime data dependency**: `data.txt` must exist for polyhedra (variant 4) to load.
3. **Makefile typo**: Variable is named `OPTIMAZE_*` (not `OPTIMIZE_*`).
4. **Unicode in printf**: Some printf strings contain Cyrillic characters in comments.
5. **X11 display required**: Program will exit with `X11_ERR_1` if no X server available (use `Xvfb` for headless testing).

## CI Configuration

GitHub Actions runs on every push to `master` and on pull requests. It simply executes `make` on `ubuntu-latest`.

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

4. **Refactor DGraphX.c** (666 lines) - This file mixes math function definitions, drawing logic, keyboard handling, and main(). Could be split into separate modules.

5. **Fix Makefile typos** - `OPTIMAZE_*` variables should be `OPTIMIZE_*`.

6. **Add stricter warning flags** - Currently `-W -Wall`. Consider adding `-Wextra`, `-Wpedantic`.

7. **Add sanitizers for development** - ASan, UBSan, MSan configs for debugging builds.

### Lower Priority

8. **Extract plottable functions to config** - Instead of editing `DGraphX.c` and recompiling, support a function definition file.

9. **Improve error handling** - Font loading failures print to stdout rather than proper error reporting. Some X11 errors aren't caught.

10. **Mouse support exists** (drag rotate/pan, right-drag zoom, scroll wheel) with event batching - could add a redraw-rate cap for very fast scroll floods.

11. **Add clang-format config** - No code formatting standard enforced.

12. **Consider XRM (X11 Resource Manager)** - For theme/font preferences instead of hardcoded Adobe Helvetica font names.
