cmake_minimum_required(VERSION 3.10)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the project name
project(Hospital
  VERSION 1.0
  DESCRIPTION "A brief experimentation for the Hospital"
  LANGUAGES CXX C
  )

# Add the executable
add_executable(Project code.cpp)

# Set the name of executable
set_target_properties(Project PROPERTIES OUTPUT_NAME "${MY_COMPILE_NAME}")

# Add the dir of the headers
include_directories(${CMAKE_SOURCE_DIR}/inc)

# Add the sources
target_sources(Project PRIVATE
  src/constants.cpp
  src/dynamics.cpp
  src/other_functions.cpp
  src/Random64.cpp
  src/reaction.cpp
  src/test.cpp
  src/workers.cpp
  src/trace.cpp
  src/other_CSV.cpp
  )

# Add the flags
target_compile_options(Project PRIVATE
  -Wall  -Werror  -g  -O3
  )

# Add the flags for the sanitizer
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=leak")
