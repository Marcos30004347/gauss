# Gauss Symbolic Math Library

A Computer Algebra Library written on C++ and ported to Javascript using WebAssembly.

![](assets/gauss.jpg)

## Installation

### Native
The recomended way of installing this library natively currently by [building it from source](#Building-from-source)

### NPM
```
npm i gauss-js
```

## Building from source

### Native build

Dependencies:
- [CMake](https://cmake.org/)
- [Make](https://www.gnu.org/software/make/)

To build native binaries from source call:

```
git clone git@github.com:Marcos30004347/gauss.git
cd gauss
make
make run-tests
```

Those commands will create a build/ folder on the root of the repository with all the compiled binaries.

### WASM/JS

Dependencies:
- [CMake](https://cmake.org/)
- [Make](https://www.gnu.org/software/make/)
- [Emscripten](https://emscripten.org/)

To build javascript and webassembly code from source call:

```
git clone git@github.com:Marcos30004347/gauss.git
cd gauss
make emsdk_path=PATH_TO_YOUR_EMSCRIPTEN_SDK wasm-binaries
```

Those commands will create a build-wasm/ folder on the root of the repository with all the generated code.

The library defines a abstraction upon the emscripten generated javascript code. You can get the full Javascript library
by calling:
```
make emsdk_path=PATH_TO_YOUR_EMSCRIPTEN_SDK release-wasm-js
```
This command will create a releases/wasm-js folder, there you will have a javascript usable library.


## Documentation

### Javascript

The JavaScript documentation is generated using jsdoc.

For javascript you will find a website documentation inside [js docs](gaussjs/docs/index.html)
you can just open that file on the browser of your choice.

### C++

The C++ documentation is generated using doxygen.

For C++ you will find a documentation website on [c++ docs](docs/gauss/html/index.html)

You can also find a documentation pdf file on [c++ pfd docs](docs/gauss/latex/refman.pdf)

For C++ documentation is still very poor and more elaborated descriptions and examples should
be added on the future.
