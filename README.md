# Gauss Symbolic Math Library

A Computer Algebra Library written on C++ and ported to the Web using Javescript and WebAssembly. I developed this project during the COVID-19 pandemic as a way to learn how algebraic systems work. The project implement state of the art algorithms and its very similar to the Sympy library available for python. But implemented in C++.

This project is currently being used in another project of mine that is my implementation of online algebraic systems such as Wolfram alpha called mathemagic.

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

For documentation check https://marcos30004347.github.io/gauss/

There you will find links to the gauss-js docs page, and the links to the gauss(c++) docs page and a a docs pdf.

### Javascript

The JavaScript documentation is generated using jsdoc.

Link to the gauss-js docs webpage.

https://marcos30004347.github.io/gauss/gaussjs/index.html

### C++

The C++ documentation is generated using doxygen.

Link to the gauss(c++) docs pdf.
https://marcos30004347.github.io/gauss/gauss/docs.pdf

Link to the gauss(c++) docs webpage.
https://marcos30004347.github.io/gauss/gauss/html/index.html
