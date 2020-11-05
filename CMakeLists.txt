cmake_minimum_required(VERSION 3.10)

#specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#set the project name
project(Hospital
  VERSION 1.0
  DESCRIPTION "A brief experimentation for the Hospital"
  LANGUAGES CXX C
  )

#add the executable
add_executable(Project code.cpp)

#set the name of executable
set_target_properties(Project PROPERTIES OUTPUT_NAME "x.code")

#add the sources
target_sources(Project PRIVATE
  dynamics.cpp
  other_functions.cpp
  Random64.cpp
  reaction.cpp
  test.cpp
  trabajadores.cpp
  trace.cpp
  )

#add the flags
target_compile_options(Project PRIVATE
  -Wall  -Werror  -g  -O3
  )

#add the sanitizers
#target_link_libraries(Project
#  -fsanitize=address
#  -fsanitize=leak
#  )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=leak")

