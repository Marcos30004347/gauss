all: environment binaries

build_type ?= Debug

environment:
	if [ ! -d "./build" ]; then mkdir build; fi
	cd build && \
	cmake .. -DCMAKE_BUILD_TYPE=$(build_type) -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DBUILD_TESTS=ON
	if [ -f "./compile_commands.json" ]; then rm -rf ./compile_commands.json; fi
	if [ -f "./build/compile_commands.json" ]; then ln ./build/compile_commands.json .; fi

binaries: SHELL:=/bin/bash
binaries:
	cmake --build ./build

run-tests:
	cd build && ctest -C Debug

wasm-binaries:
	if [ ! -d build-wasm ]; then \
		mkdir build-wasm; \
		cd build-wasm && \
		cmake .. -DCMAKE_BUILD_TYPE=$(build_type) -DBUILD_WASM=ON \
		-DEMSDK_PATH=$(emsdk_path) \
-DCMAKE_TOOLCHAIN_FILE=$(emsdk_path)/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake \
		-DCMAKE_CROSSCOMPILING_EMULATOR=$(emsdk_path)/node/14.18.2_64bit/bin/node; \
	fi
	cmake --build ./build-wasm --config $(build_type)

clean:
	if [ -d "./build" ]; then      rm -rf ./build;      fi
	if [ -d "./build-wasm" ]; then rm -rf ./build-wasm; fi
	if [ -d "./releases" ]; then   rm -rf ./releases;   fi

release-wasm-js: wasm-binaries
	if [ ! -d "./releases/wasm-js" ]; then mkdir -p releases/wasm-js; fi

	cp ./build-wasm/gaussjs.js     ./releases/wasm-js/gaussjs.js
	cp ./build-wasm/gaussjs.wasm   ./releases/wasm-js/gaussjs.wasm
	cp ./gaussjs/package.json      ./releases/wasm-js/package.json
	cp ./gaussjs/package-lock.json ./releases/wasm-js/package-lock.json
	cp ./gaussjs/README.md         ./releases/wasm-js/README.md
	cp ./gaussjs/index.js          ./releases/wasm-js/index.js

	rm -rf build-wasm

release-linux: environment binaries
	if [ ! -d "./releases/linux" ]; then mkdir -p releases/linux; fi

	cp -r ./build/libgauss.a  ./releases/linux/libgauss.a
	cp -r ./build/include     ./releases/linux/include

	rm -rf build

release-windows: environment binaries
	if [ ! -d "./releases/windows" ]; then mkdir -p releases/windows; fi

	cp -r ./build/x64  		./releases/windows/x64
	cp -r ./build/include   ./releases/windows/include

	rm -rf build

release-macos: environment binaries
	if [ ! -d "./releases/macos" ]; then mkdir -p releases/macos; fi

	cp -r ./build/libgauss.a  ./releases/macos/libgauss.a
	cp -r ./build/include     ./releases/macos/include

	rm -rf build


# releases: release-linux-x86 release-wasm
