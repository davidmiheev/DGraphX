# DGraphX Project Guide

## Overview
DGraphX is a lightweight, dependency-free 2D and 3D graphing application written in C. It uses the X11 windowing system for rendering and is designed to run efficiently even on older hardware without graphics acceleration.

### Core Technologies
- **Language:** C (C99/C11 compatible)
- **Graphics:** X11 (Xlib)
- **Build System:** Makefile
- **Math:** Custom matrix operations and interpolation logic

## Getting Started

### Prerequisites
- X11 development headers and libraries (typically `libx11-dev` on Linux)
- A C compiler (GCC or Clang)

### Building and Running
```bash
# Build the project
make

# Run the application
./DGraphX

# Clean build artifacts
make clean
```

## Usage and Interaction
The application is entirely keyboard-driven. There is currently no mouse support.

### 3D Mode Controls
- **Rotation:** `<-` and `->` keys.
- **Switch Rotation Axis:** `X`, `Y`, or `Z` keys.
- **Camera Distance:** `Up` and `Down` keys.
- **Zoom:** `C` (In), `V` (Out).
- **Toggle Axes:** `Enter`.
- **Light Source:** `l` to enter/exit light adjustment mode, then `<-` and `->` to move.
- **Color Method:** `L` (Shift + L).
- **Quit:** `q`.

### 2D Mode Controls
- **Toggle Mode:** `Ctrl` switches between 2D and 3D modes.
- **Views:** Press `1`, `2`, `3`, or `4` to switch between different 2D visualization variants (interpolation, functions, parametric curves).

## Development Guide

### Changing Plotted Functions
To change the function being graphed, you must modify the source code and recompile. Look for the following functions in `DGraphX.c`:
- `f(x)`: Standard 2D function.
- `u(x, y)`, `v(x, y)`, `w(x, y)`: 3D surfaces.
- `parx(t, s)`, `pary(t, s)`, `parz(t, s)`: Parametric 3D surfaces.

### Project Structure
- `DGraphX.c`: Main logic, drawing functions, and keyboard event handling.
- `Xfun.c`: X11 window management and drawing primitives.
- `DGraphX.h`: Global state declarations, color macros, and type definitions.
- `matrix.c/h`: Linear algebra utilities for 3D transformations.
- `interpolation.c`: 2D interpolation algorithms (Splines, Hermite).
- `DX.c`: 3D projection and transformation logic.
- `data.txt`: Data for polyhedra vertices (required for certain 3D variants).

### Coding Conventions
- **Naming:** 
  - `CamelCase` for X11 wrapper functions (e.g., `WDrawString`, `WSetColor`).
  - `snake_case` for internal logic and math utilities.
  - Hungarian-style prefixes are sometimes used for global pointers (e.g., `prDisplay`, `nWnd`).
- **Global State:** Global variables are declared in `DGraphX.h` as `extern` and defined in `Xfun.c`.
- **Colors:** Use the macros defined in `DGraphX.h` (e.g., `RED`, `BLUE`) which map to the `rColors` global struct.

### Testing
There is no automated test suite. Validation is performed manually by running the executable and verifying the visual output.

## Roadmap & Future Improvements

### User Experience
- **Runtime Function Input:** Integrate a math expression parser (e.g., tinyexpr) to allow plotting functions without recompilation.
- **Mouse Support:** Implement X11 mouse event handling for intuitive 3D rotation and 2D panning.
- **CLI Arguments:** Support flags for initial mode (2D/3D), data file paths, and camera settings.

### Technical Debt & Code Quality
- **Makefile Portability:** Replace hardcoded X11 paths with `pkg-config` for cross-distro compatibility.
- **State Management:** Encapsulate global variables into a central context struct to improve testability and maintainability.
- **Robust Error Handling:** Add input validation for `data.txt` and other external resources.

### Infrastructure
- **Automated Testing:** Develop a suite for math utilities (`matrix.c`, `interpolation.c`).
- **Headless CI:** Use `Xvfb` in GitHub Actions to perform basic smoke tests on the application.
- **Cross-Platform Support:** Consider an abstraction layer (SDL2/GLFW) to move beyond X11-only systems.
