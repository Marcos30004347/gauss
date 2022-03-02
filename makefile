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
