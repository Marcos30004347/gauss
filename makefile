all: build-conf

build-conf:
	if [ -d "./build" ]; then rm -rf ./build; fi
	mkdir build
	cd build && \
	cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1
	if [ -f "./compile_commands.json" ]; then rm -rf ./compile_commands.json; fi
	ln ./build/compile_commands.json .

wasm:
	$(info VAR="$(WASI_SDK_PATH)")

	if [ -d "./build" ]; then rm -rf ./build; fi
	mkdir build
	cd build && \
	cmake .. -DWASI_SDK_PATH=$(wasi_sdk_path) -DBUILD_WASM=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=1

gaussjs:
	if [ -d "./build" ]; then rm -rf ./build; fi
	mkdir build
	cd build && \
  cmake .. -DBUILD_WASM=ON \
	-DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
	-DEMSDK_PATH=$(emsdk_path) \
	-DCMAKE_TOOLCHAIN_FILE=$(emsdk_path)/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake \
	-DCMAKE_CROSSCOMPILING_EMULATOR=$(emsdk_path)/node/14.18.2_64bit/bin/node
	if [ -f "./compile_commands.json" ]; then rm -rf ./compile_commands.json; fi
	ln ./build/compile_commands.json .
	cd build && make
	cd build && mkdir gaussjs && cp gaussjs.wasm.js gaussjs/gaussjs_wasm.js && cp gaussjs.wasm.wasm gaussjs/gaussjs_wasm.wasm
	cp ./src/WebAssembly/js/gauss.js ./build/gaussjs/gauss.js
	cp ./src/WebAssembly/js/server.js ./build/gaussjs/server.js
	cp ./src/WebAssembly/js/index.html ./build/gaussjs/index.html
