all: environment binaries

build_type ?= Debug

environment:
	if [ ! -d "./build" ]; then mkdir build; fi
	cd build && \
	cmake .. -DCMAKE_BUILD_TYPE=$(build_type) -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DBUILD_TESTS=ON
	if [ -f "./compile_commands.json" ]; then rm -rf ./compile_commands.json; fi
	if [ -f "./build/compile_commands.json" ]; then ln ./build/compile_commands.json .; fi

binaries:
	cmake --build ./build --config $(build_type)

run-tests:
	cd build && ctest

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

release-native: environment binaries
	if [ ! -d "./releases/native" ]; then mkdir -p releases/native; fi

	cp -r ./build/libgauss.a  ./releases/native/libgauss.a
	cp -r ./build/include     ./releases/native/include

	rm -rf build

releases: release-native release-wasm
