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

clean:
	if [ -d "./build" ]; then rm -rf ./build; fi
	if [ -d "./gauss_js_build" ]; then rm -rf ./gauss_js_build; fi
