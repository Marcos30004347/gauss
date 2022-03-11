all: environment binaries

environment:
	if [ ! -d "./build" ]; then mkdir build; fi
	cd build && \
	cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DBUILD_TESTS=ON
	if [ -f "./compile_commands.json" ]; then rm -rf ./compile_commands.json; fi
	ln ./build/compile_commands.json .

binaries:
	cd build && make

run-tests:
	cd build && ctest

wasm-binaries:
	if [ ! -d build-wasm ]; then \
		mkdir build-wasm; \
		cd build-wasm && \
		cmake .. -DBUILD_WASM=ON \
		-DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
		-DEMSDK_PATH=$(emsdk_path) \
		-DCMAKE_TOOLCHAIN_FILE=$(emsdk_path)/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake \
		-DCMAKE_CROSSCOMPILING_EMULATOR=$(emsdk_path)/node/14.18.2_64bit/bin/node; \
	fi

	cd build-wasm && make

	cp ./build-wasm/gaussjs.js ./gaussjs/gaussjs.js
	cp ./build-wasm/gaussjs.wasm ./gaussjs/gaussjs.wasm

clean:
	if [ -d "./build" ]; then rm -rf ./build; fi
	if [ -d "./build-wasm" ]; then rm -rf ./build-wasm; fi
