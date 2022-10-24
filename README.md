# DGraphX
This program is simple and cute 2D and 3D grapher (to draw various 2d and 3d plots) from scratch. It doesn't have any dependencies at all. 
For creating graphics it uses only homebrewed ```C``` сode and ```X11``` (main windowing system on UNIX-like systems). 
```DGraphX``` is a very lightweight program and is capable work on a very slow devices (without any graphics accelerators). 
## Requirements
* ```X11```
* ```С``` compiler (e.g gcc or clang)
## Install
Clone this repo and make the executable file:
```
git clone https://github.com/DavidOSX/DGraphX
cd DGraphX
make
```
## Usage
Run executable and enjoy:
```
./DGraphX
```
After that you will be dealing with the GUI. The GUI is very modest, but (hope) simple for understanding.
Important keybord commands (program yet don't support mouse) for the 3D mode:
* Use **<-** and **->** keys for rotation about one of the coordinate axes 
* By default the plot rotates about Z axis, but you can change axis: **X** key change to rotation about X axis, **Y** key rotation about Y axis, **Z** key rotation about Z axis. 
* If you want to see axes, just press **Enter**
* Use **Up** and **Down** for changing the camera distance from the coordinate origin.
* **C** key for zoom in, **V** key for zoom out.
* You can also change the light source position: first press **l** key and change position of the light source with **<-** and **->** keys, then again press **l** key
* You can also change the color method with **L** key (not **l** key)
* To quit press **q**

Important keybord commands (program yet don't support mouse) for the 2D mode:
Coming soon...

## Further Comments
Sad but for choose what specific plot you want to draw, you need to correct code of functions in ```DGraphX.c``` and recompiling, it need to fix. 
Also there is some wish to add the mouse support and rewritten program with the more modern technologies...

