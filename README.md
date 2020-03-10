# Solar System Simulator
A realistic simulator of our solar system based on Kepler's laws of planetary motion.

Features:
 * Realistic modeling of celestial objects:
   * Physical characteristics such as mean radius and axial tilt
   * Orbital characteristics such as eccentricity, inclination, rotation period, and orbital period
 * Point representation of Kuiper Belt Objects (KBOs)
 * Planet textures
 * Complete time flow control
 * Third-person camera that can be used to observe planets:
   * from a fixed point
   * following the planet's rotation

This project constituted a semester assignment for the "Computer Graphics and Virtual Reality" course at the Department of Electrical and Computer Engineering of the University of Patras.

## Usage
### Prerequisites
 * [MathGeoLib](https://github.com/juj/MathGeoLib) (tested with commit 2103ced)
 * [SOIL](https://github.com/paralin/soil) (tested with commit 8bb18a9)

### Compilation
```shell
$ MATHGEOLIB_INCLUDE_DIR=/path/to/MathGeoLib/include MATHGEOLIB_LIBRARY=/path/to/MathGeoLib/lib SOIL_INCLUDE_DIR=/path/to/SOIL/include SOIL_LIBRARY=/path/to/SOIL/lib make -j$(nproc)
```

### Execution
```shell
$ ./ss
```

### Cleanup
```shell
$ make clean
```

## Demo
![Demo](../assets/demos/demo.gif)

## Screenshots
|![](../assets/screenshots/screenshot_1.png?raw=true)|![](../assets/screenshots/screenshot_2.png?raw=true)|![](../assets/screenshots/screenshot_3.png?raw=true)|
|---|---|---|
|![](../assets/screenshots/screenshot_4.png?raw=true)|![](../assets/screenshots/screenshot_5.png?raw=true)|![](../assets/screenshots/screenshot_6.png?raw=true)|
|![](../assets/screenshots/screenshot_7.png?raw=true)|![](../assets/screenshots/screenshot_8.png?raw=true)|![](../assets/screenshots/screenshot_9.png?raw=true)|
