cmake_minimum_required(VERSION 3.10)

project (gauss)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


option(BUILD_WASM "build wasm binaries" OFF)
option(BUILD_TESTS "build tests" OFF)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

	set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -std=c++0x")
	set(CMAKE_CXX_FLAGS_DEBUG "-g")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

  set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -std=c++0x")
  set(CMAKE_CXX_FLAGS_DEBUG "-g")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /std:c11 /wd4290 /wd4244")

endif()

set(GAUSS_SOURCES
	gauss/Error/error.cpp
	gauss/Algebra/Matrix.cpp
  gauss/Algebra/Integer.cpp
  gauss/Algebra/Expression.cpp
  gauss/Algebra/Utils.cpp
  gauss/Algebra/Reduction.cpp
  gauss/Algebra/Sorting.cpp
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
  gauss/Polynomial/Roots.cpp
  gauss/SVD/SVD.cpp
  gauss/Gauss.cpp
)

set(GAUSS_INCLUDES
	gauss/Error/error.hpp
	gauss/Algebra/Matrix.hpp
  gauss/Algebra/Integer.hpp
  gauss/Algebra/Expression.hpp
  gauss/Algebra/Utils.hpp
  gauss/Algebra/Reduction.hpp
  gauss/Algebra/Sorting.hpp
  gauss/Algebra/Trigonometry.hpp
  gauss/Algebra/Lexer.hpp
  gauss/Algebra/Parser.hpp
  gauss/Algebra/Int.hpp
  gauss/Calculus/Derivative.hpp
  gauss/Primes/Primes.hpp
  gauss/Factorization/Utils.hpp
  gauss/Factorization/Hensel.hpp
  gauss/Factorization/Berlekamp.hpp
  gauss/Factorization/SquareFree.hpp
  gauss/Factorization/Zassenhaus.hpp
  gauss/Factorization/Wang.hpp
  gauss/GaloisField/GaloisField.hpp
  gauss/Polynomial/Resultant.hpp
  gauss/Polynomial/Polynomial.hpp
  gauss/Polynomial/Roots.hpp
  gauss/SVD/SVD.hpp
  gauss/Gauss.hpp
)

foreach(HEADER ${GAUSS_INCLUDES})
  string(REGEX REPLACE "[^\/]+$" "" PATH ${HEADER})

  message("Copying ${HEADER} to ${CMAKE_CURRENT_BINARY_DIR}/include/${PATH}")

  file(COPY ${HEADER} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/${PATH})
endforeach()

add_library(gauss STATIC ${GAUSS_SOURCES})

target_include_directories(gauss PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/include")

if(BUILD_WASM)
	project(gaussjs)

	add_executable(gaussjs gaussjs/bindings.cpp)

  set_target_properties(gaussjs PROPERTIES PREFIX "")
  set_target_properties(gaussjs PROPERTIES OUTPUT_NAME "gauss-bindings")
  set_target_properties(gaussjs PROPERTIES SUFFIX ".js")

	#TODO: enable fsanitize only on debug builds

  set_target_properties(gaussjs PROPERTIES LINK_FLAGS " \
	  --bind \
    -s WASM=1 \
    -s TOTAL_MEMORY=1024MB \
    -s MALLOC=emmalloc \
	  -s MODULARIZE=1 \
    -s SINGLE_FILE=1 \
    -O3 -g0 \
  ")
  # -s ALLOW_MEMORY_GROWTH=1 \
  # -fsanitize=address \

	target_link_libraries(gaussjs gauss)
endif()

if(BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()
