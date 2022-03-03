cmake_minimum_required(VERSION 3.10)

project (math)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


option(BUILD_WASM "build wasm binaries" OFF)

if(BUILD_WASM)
	message("BUILDING WASM")
	set( CMAKE_VERBOSE_MAKEFILE on )
	set(CMAKE_SYSTEM_NAME WASI)
	set(CMAKE_SYSTEM_VERSION 1)
	set(CMAKE_SYSTEM_PROCESSOR wasm32)

	set(triple wasm32-unknown-unknown)

	if(WIN32)
		set(WASI_HOST_EXE_SUFFIX ".exe")
	else()
		set(WASI_HOST_EXE_SUFFIX "")
	endif()

	message("WASI SDK PATH = ${WASI_SDK_PATH}")

	set(CMAKE_C_COMPILER ${WASI_SDK_PATH}/bin/clang${WASI_HOST_EXE_SUFFIX})
	set(CMAKE_CXX_COMPILER ${WASI_SDK_PATH}/bin/clang++${WASI_HOST_EXE_SUFFIX})

	set(CMAKE_SYSROOT ${WASI_SDK_PATH}/share/wasi-sysroot)

	set(CMAKE_AR ${WASI_SDK_PATH}/bin/llvm-ar${WASI_HOST_EXE_SUFFIX})
	set(CMAKE_RANLIB ${WASI_SDK_PATH}/bin/llvm-ranlib${WASI_HOST_EXE_SUFFIX})

	set(CMAKE_C_COMPILER_TARGET ${triple})
	set(CMAKE_CXX_COMPILER_TARGET ${triple})

	# Don't look in the sysroot for executables to run during the build
	set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
	# Only look in the sysroot (not in the host paths) for the rest
	set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
	set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
	set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

	set(CMAKE_CXX_FLAGS_DEBUG "-g")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")


	set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++0x -nostartfiles \
  -fno-exceptions \
  -Wl,--strip-all \
  -Wl,--export-dynamic \
  -Wl,--export-all \
  -Wl,--import-memory \
	-Wl,--allow-undefined \
  -fvisibility=hidden")


	set(WASM_MAIN_FILE "src/WebAssembly/main.cpp")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++0x")
	set(CMAKE_CXX_FLAGS_DEBUG "-g")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++0x")
	set(CMAKE_CXX_FLAGS_DEBUG "-g")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

if(NOT BUILD_WASM)
	set(WASM_MAIN_FILE "")
endif()
message("USING CXX compiler = ${CMAKE_CXX_COMPILER}")

add_library(
  math
  src/MathSystem/Algebra/String.cpp
  src/MathSystem/Algebra/Matrix.cpp
  src/MathSystem/Algebra/Expression.cpp
  src/MathSystem/Algebra/Integer.cpp
  src/MathSystem/Algebra/Int.cpp
  # src/MathSystem/Trigonometry/Trigonometry.cpp
  # src/MathSystem/Exponential/Exponential.cpp
  src/MathSystem/Calculus/Derivative.cpp
    # src/MathSystem/Calculus/Integral.cpp
    # src/MathSystem/Rational/Rational.cpp
    src/MathSystem/Primes/Primes.cpp
    # src/MathSystem/Expand/Expand.cpp
    # src/MathSystem/Expand/Multiplication.cpp
    # src/MathSystem/Expand/Polynomial.cpp
    # src/MathSystem/Expand/Division.cpp
    # src/MathSystem/Simplification/Rationals.cpp
    # src/MathSystem/Simplification/Power.cpp
    # src/MathSystem/Simplification/Multiplication.cpp
    # src/MathSystem/Simplification/Addition.cpp
    # src/MathSystem/Simplification/Subtraction.cpp
    # src/MathSystem/Simplification/Division.cpp
    # src/MathSystem/Simplification/Simplification.cpp
    # src/MathSystem/Simplification/Factorial.cpp
    # src/MathSystem/Simplification/Trigonometry.cpp
    # src/MathSystem/Simplification/Exponential.cpp

    src/MathSystem/Factorization/Utils.cpp
    src/MathSystem/Factorization/Hensel.cpp
		src/MathSystem/Factorization/Berlekamp.cpp
		src/MathSystem/Factorization/SquareFree.cpp
    src/MathSystem/Factorization/Zassenhaus.cpp
    src/MathSystem/Factorization/Wang.cpp

		# #src/MathSystem/Polynomial/Factorization.cpp

    src/MathSystem/GaloisField/GaloisField.cpp
		src/MathSystem/Polynomial/Resultant.cpp
		src/MathSystem/Polynomial/Polynomial.cpp

		# src/MathSystem/Polynomial/Algebraic.cpp

		# #src/MathSystem/Polynomial/Hensel.cpp
		src/TypeSystem/Utils.cpp
		src/TypeSystem/AST.cpp
		src/TypeSystem/SymbolTable.cpp
		${WASM_MAIN_FILE}
	)

if(BUILD_WASM)
	target_compile_definitions(math PRIVATE WASM_BUILD=1)

	set_target_properties(
		  math
			PROPERTIES
				  SUFFIX ".wasm"
					)

	# project(math_wasm_lib)

	# add_library(math_wasm_lib src/WebAssembly.main.cpp)

	#target_compile_definitions(math_wasm_lib PRIVATE WASM_BUILD=1)
	#target_link_libraries(math_wasm_lib PUBLIC math)

endif()

target_include_directories(
    math
    PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

if(NOT BUILD_WASM)
	enable_testing()
	add_subdirectory(tests)
endif()
