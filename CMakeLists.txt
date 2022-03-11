cmake_minimum_required(VERSION 3.10)

project (math)

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
  math STATIC
  src/MathSystem/Algebra/Matrix.cpp
  src/MathSystem/Algebra/Expression.cpp
  src/MathSystem/Algebra/Integer.cpp
  src/MathSystem/Algebra/Lexer.cpp
  src/MathSystem/Algebra/Parser.cpp
  src/MathSystem/Algebra/Int.cpp
  src/MathSystem/Calculus/Derivative.cpp
  src/MathSystem/Primes/Primes.cpp
  src/MathSystem/Factorization/Utils.cpp
  src/MathSystem/Factorization/Hensel.cpp
	src/MathSystem/Factorization/Berlekamp.cpp
	src/MathSystem/Factorization/SquareFree.cpp
  src/MathSystem/Factorization/Zassenhaus.cpp
  src/MathSystem/Factorization/Wang.cpp
	src/MathSystem/GaloisField/GaloisField.cpp
	src/MathSystem/Polynomial/Resultant.cpp
	src/MathSystem/Polynomial/Polynomial.cpp
	src/MathSystem/SVD/SVD.cpp
	src/MathSystem/Gauss.cpp
)

target_include_directories(math PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

if(BUILD_WASM)
	project(gaussjs)

	add_executable(gaussjs gaussjs/index.cpp)

  set_target_properties(gaussjs PROPERTIES LINK_FLAGS "--emrun \
	--bind \
	-s WASM=1 \
	-s ALLOW_MEMORY_GROWTH=1\
	-s MALLOC=emmalloc \
	-s MODULARIZE=1")

	target_link_libraries(gaussjs math)
	target_compile_definitions(gaussjs PRIVATE WASM_BUILD=1)
endif()

if(BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()
