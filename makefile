all: environment binaries

build_type ?= Debug

ifeq ($(OS),Windows_NT)
	OPERATIONAL_SYSTEM := WINDOWS
	ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
		PROCESSOR := amd64
	endif
	ifeq ($(PROCESSOR_ARCHITECTURE),x86)
		PROCESSOR := x86
	endif
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		OPERATIONAL_SYSTEM := LINUX
	endif

	ifeq ($(UNAME_S),Darwin)
		OPERATIONAL_SYSTEM := MACOS_DARWIN
	endif

	UNAME_P := $(shell uname -p)

	ifeq ($(UNAME_P),x86_64)
		PROCESSOR := x86_64
	endif

	ifneq ($(filter %86,$(UNAME_P)),)
		PROCESSOR := x86
	endif

	ifneq ($(filter arm%,$(UNAME_P)),)
		PROCESSOR := arm
	endif
endif



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

	cp ./build-wasm/gauss.cjs      ./releases/wasm-js/gauss.cjs
	cp ./build-wasm/gauss.wasm     ./releases/wasm-js/gauss.wasm
	cp ./gaussjs/package.json      ./releases/wasm-js/package.json
	cp ./gaussjs/package-lock.json ./releases/wasm-js/package-lock.json
	cp ./gaussjs/README.md         ./releases/wasm-js/README.md
	cp ./gaussjs/gauss.js          ./releases/wasm-js/gauss.js
	cp -r ./gaussjs/docs              ./releases/wasm-js/docs
	cp ./gaussjs/README.md         ./releases/wasm-js/README.md
	cp ./gaussjs/jsconfig.json     ./releases/wasm-js/jsconfig.json
	cp ./gaussjs/jsdoc.json				 ./releases/wasm-js/jsdoc.json
	cp ./gaussjs/makefile				   ./releases/wasm-js/makefile

# rm -rf build-wasm


release-local: environment binaries
ifeq ($(OPERATIONAL_SYSTEM), LINUX)
	if [ ! -d "./releases/linux-$(PROCESSOR)" ]; then mkdir -p releases/linux-$(PROCESSOR); fi

	cp -r ./build/libgauss.a  ./releases/linux-$(PROCESSOR)/libgauss.a
	cp -r ./build/include     ./releases/linux-$(PROCESSOR)/include

	rm -rf build
else ifeq ($(OPERATIONAL_SYSTEM), MACOS_DARWIN)
	if [ ! -d "./releases/macos-darwin-$(PROCESSOR)" ]; then mkdir -p releases/macos-darwin-$(PROCESSOR); fi

	cp -r ./build/libgauss.a  ./releases/macos-darwin-$(PROCESSOR)/libgauss.a
	cp -r ./build/include     ./releases/macos-darwin-$(PROCESSOR)/include

	rm -rf build
else ifeq ($(OPERATIONAL_SYSTEM), WINDOWS)
	if [ ! -d "./releases/windows-$(PROCESSOR)" ]; then mkdir -p releases/windows-$(PROCESSOR); fi

	cp -r ./build/include   ./releases/windows-$(PROCESSOR)/include

	if [ -d "./build/x86" ]; then \
		cp -r ./build/x86  		./releases/windows-$(PROCESSOR)/x86; \
	fi

	if [ -d "./build/x64" ]; then \
		cp -r ./build/x64  		./releases/windows-$(PROCESSOR)/x64; \
	fi

	if [ -d "./build/Debug" ]; then \
		cp -r ./build/Debug  		./releases/windows-$(PROCESSOR)/Debug; \
	fi

	if [ -d "./build/Release" ]; then \
		cp -r ./build/Release  		./releases/windows-$(PROCESSOR)/Release; \
	fi

	rm -rf build
else
	$(error Handler for release not registered for OS $(OPERATIONAL_SYSTEM))
endif


# releases: release-linux-x86 release-wasm
