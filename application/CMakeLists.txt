cmake_minimum_required(VERSION 3.31)
project(power CXX)

# user options
option(DISABLE_CUDA "Disable CUDA")
option(ENABLE_WARNINGS "Enable strict warnings")

# by default build optimized code with debug symbols
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()
if(NOT CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES native)
endif()

set(CMAKE_ALLOWED_BUILD_TYPES DEBUG RELEASE RELWITHDEBINFO)
if(NOT CMAKE_BUILD_TYPE IN_LIST CMAKE_ALLOWED_BUILD_TYPES)
  message(FATAL_ERROR "Invalid build type ${CMAKE_BUILD_TYPE}. Use one of: ${CMAKE_ALLOWED_BUILD_TYPES}.")
endif()

# target setup
add_executable(project
  main.cpp
  generator.cpp
  detector.cpp
  reconstruction.cpp
  reader.cpp
  cmodules/qconfig.cpp
  cmodules/timer.cpp
)

# warn about undefined symbols already for libraries
if(ENABLE_WARNINGS)
  SET_TARGET_PROPERTIES(project PROPERTIES LINK_FLAGS "-Wl,--as-needed")
  SET_TARGET_PROPERTIES(project PROPERTIES LINK_FLAGS "-Wl,--no-undefined")
  set(CXX_FLAGS_WARNINGS "${CXX_FLAGS_WARNINGS} -Wall")
  set(CXX_FLAGS_WARNINGS "${CXX_FLAGS_WARNINGS} -Wextra")
  set(CXX_FLAGS_WARNINGS "${CXX_FLAGS_WARNINGS} -Werror")
endif()
set(CXX_FLAGS_WARNINGS "${CXX_FLAGS_WARNINGS} -Wno-unused-function -Wno-unused-parameter")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -DDEBUG -D__DEBUG ${CXX_FLAGS_WARNINGS}")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -DNODEBUG ${CXX_FLAGS_WARNINGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -ggdb")
message(STATUS "CXX Flags: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")


if(DISABLE_CUDA STREQUAL "OFF")
  if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    set(CMAKE_CUDA_FLAGS "-Xptxas -O0 -Xcompiler -O0")
  else()
    set(CMAKE_CUDA_FLAGS "-Xptxas -O4 -Xcompiler -O4")
  endif()
  enable_language(CUDA)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --compiler-options='${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}} -DUSE_GPU'")
  set_property(TARGET project PROPERTY CUDA_SEPARABLE_COMPILATION OFF)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "Using C++${CMAKE_CXX_STANDARD}")

target_include_directories(project PRIVATE cmodules)
include_directories(.)

get_property(LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)
if (CUDA IN_LIST LANGUAGES)
  message(STATUS "Using CUDA: ${CUDA_VERSION}")
  set_property(TARGET project APPEND PROPERTY COMPILE_DEFINITIONS USE_GPU)
  target_sources(project PRIVATE gpu.cu)
  set_property(TARGET project PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()

INSTALL(FILES defaultSim.txt DESTINATION .)
INSTALL(TARGETS project DESTINATION .)
