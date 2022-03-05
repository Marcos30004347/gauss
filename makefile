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
	if [ ! -d gauss_js_build ]; then \
		mkdir gauss_js_build; \
		cd gauss_js_build && \
		cmake .. -DBUILD_WASM=ON \
		-DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
		-DEMSDK_PATH=$(emsdk_path) \
		-DCMAKE_TOOLCHAIN_FILE=$(emsdk_path)/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake \
		-DCMAKE_CROSSCOMPILING_EMULATOR=$(emsdk_path)/node/14.18.2_64bit/bin/node; \
	fi

	cd gauss_js_build && make

	cp ./src/WebAssembly/js/gauss.js ./gauss_js_build/gauss.js
	cp ./src/WebAssembly/js/server.js ./gauss_js_build/server.js
	cp ./src/WebAssembly/js/index.html ./gauss_js_build/index.html
#cp ./src/WebAssembly/js/package.json ./gauss_js_build/package.json
clean:
	if [ -d "./build" ]; then rm -rf ./build; fi
	if [ -d "./gauss_js_build" ]; then rm -rf ./gauss_js_build; fi
