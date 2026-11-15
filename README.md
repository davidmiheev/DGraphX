# DGraphX
This program is a simple and cute 2D and 3D grapher (to draw various 2d and 3d plots) from scratch. It doesn't have any dependencies at all. 
For creating graphics it uses only homebrewed ```C``` сode and ```X11``` (main windowing system on UNIX-like systems). 
Grapher ```DGraphX``` is a very lightweight program and is capable of working on very slow devices (without any graphics accelerators).

## Requirements
* X11
* С compiler (e.g gcc or clang)

## Install
Clone this repo and make the executable file:
```sh
git clone https://github.com/davidmiheev/DGraphX
cd DGraphX
make
```

## Usage
Run executable and enjoy:

```sh
./DGraphX
```

After that, you will be dealing with the GUI. The GUI is very modest, but (I hope) simple to understand.
Press **i** at any time to type a new function expression in the terminal (tinyexpr syntax, e.g. `sin(x)*x`).
Important keyboard commands for the 3D mode:
* Use **<-** and **->** keys for rotation about one of the coordinate axes 
* By default, the plot rotates about the Z axis, but you can change the axis: **X** key change to rotation about the X axis, **Y** key rotation about the Y axis, and **Z** key rotation about the Z axis. 
* If you want to see axes, just press **Enter**
* Use **Up** and **Down** for changing the camera distance from the coordinate origin.
* **C** key for zoom in, **V** key for zoom out.
* You can also change the light source position: first press **l** key and change the position of the light source with **<-** and **->** keys, then again press **l** key
* You can also change the color method with **L** key (not **l** key)
* **R** key toggles progressive refinement (half-resolution rendering while interacting, full resolution when idle)
* **P** key toggles Phong shading (smooth per-pixel normals) vs flat facet shading
* **G** key toggles the floor grid (coordinate level curves under the 3D plot)
* To quit press **q**

Mouse commands for the 3D mode:
* **Left-drag** rotates the plot
* **Right-drag** changes the camera distance
* **Scroll wheel** zooms in/out

Important keyboard commands for the 2D mode:
Coming soon...

Mouse commands for the 2D mode:
* **Left-drag** pans the view
* **Right-drag** or **scroll wheel** zooms in/out around the view center

## Further Comments

Functions to plot can be changed at runtime: press **i** in the window and enter a new expression in the terminal (no recompilation needed).
There is also some wish to rewrite the program with the more modern technologies...
