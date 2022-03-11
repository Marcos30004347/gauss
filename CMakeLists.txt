cmake_minimum_required(VERSION 3.10)

project (gauss)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


option(BUILD_WASM "build wasm binaries" OFF)
option(BUILD_TESTS "build tests" OFF)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

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
  gauss STATIC
  gauss/Algebra/Matrix.cpp
  gauss/Algebra/Integer.cpp
  gauss/Algebra/Expression.cpp
  gauss/Algebra/Trigonometry.cpp
  gauss/Algebra/Lexer.cpp
  gauss/Algebra/Parser.cpp
  gauss/Algebra/Int.cpp
  gauss/Calculus/Derivative.cpp
  gauss/Primes/Primes.cpp
  gauss/Factorization/Utils.cpp
  gauss/Factorization/Hensel.cpp
  gauss/Factorization/Berlekamp.cpp
  gauss/Factorization/SquareFree.cpp
  gauss/Factorization/Zassenhaus.cpp
  gauss/Factorization/Wang.cpp
  gauss/GaloisField/GaloisField.cpp
  gauss/Polynomial/Resultant.cpp
  gauss/Polynomial/Polynomial.cpp
  gauss/SVD/SVD.cpp
  gauss/Gauss.cpp
)

target_include_directories(gauss PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

if(BUILD_WASM)
	project(gaussjs)

	add_executable(gaussjs gaussjs/index.cpp)

  set_target_properties(gaussjs PROPERTIES LINK_FLAGS "--emrun \
	--bind \
	-s WASM=1 \
	-s ALLOW_MEMORY_GROWTH=1\
	-s MALLOC=emmalloc \
	-s MODULARIZE=1")

	target_link_libraries(gaussjs gauss)
	target_compile_definitions(gaussjs PRIVATE WASM_BUILD=1)
endif()

if(BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()
