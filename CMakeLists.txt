if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()
#
# ComPWA - A Common framework for Partial Wave Analysis
#
project( COMPWA 
  VERSION 1.0.0
  DESCRIPTION "The common Partial Wave Analysis framework"
  LANGUAGES CXX )

#
# CMake configuration
#
cmake_minimum_required (VERSION 3.4 FATAL_ERROR)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Which target should be build by default?
set(DEFAULT_BUILD_TYPE "Release")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#
# Config options
#
option(SAFE_LINKING "Enable to avoid conflicts between different versions \
of libraries stdc++, boost, ...)" OFF)

# Setting policy CMP0060 to the new version enforces that libraries are linked
# using their full path. That should help in case that multiple versions of a
# library (e.g. boost) are installed on a system
cmake_policy(SET CMP0060 NEW)

# Configure RPATH
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
    STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Additonal compile flags for various compilers
# ${CMAKE_CXX_COMPILER_ID} can be one of {GNU Clang AppleClang Intel MSVC}
# - verbose output on loop vectorization
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang") # using Clang or AppleClang
  #set(CMAKE_CXX_FLAGS
  #    "${CMAKE_CXX_FLAGS} -Rpass-analysis=loop-vectorize -Rpass=loop-vectorize")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") # using GCC
  # gcc 4.8/4.9 have for example no full regex implementation 
  # and will cause runtime errors
  # full c++11 support is only given in gcc 5.1
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.1")
    message(FATAL_ERROR "Version of gcc is too low, and does not have \
    full c++11 support. Please install gcc 5.1 or higher.")
  endif()
  #set(CMAKE_CXX_FLAGS
  #    "${CMAKE_CXX_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=1")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") # using Intel C++
  # nothing here
endif()

# Enable all warnings
# This is not good code, in the sense that its not portable (across compilers)
# The portable way using target_compile_options is just so inconvienient, that
# we use this for now. 
list(APPEND CMAKE_CXX_FLAGS "-Wall")

# Enable -fPIC flag
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

message( STATUS "Global CXX compiler flags: " ${CMAKE_CXX_FLAGS} )

# Setting name prefix for libraries
set( CMAKE_SHARED_LIBRARY_PREFIX "libComPWA_" )


# Code Coverage Configuration
option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()


#
# External dependencies
#
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if(SAFE_LINKING)
  if(NOT APPLE)
    # Setting RPATH instead of RUNPATH
    # RPATH is search before LD_LIBRARY_PATH at runtime
    set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--disable-new-dtags")
  endif()
  if(NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang") # not Clang or AppleClang
    message(STATUS "Option SAFE_LINKING includes static linking of stdc++ "
                   "GCC.")
    # Link libstdc++ statically
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc")
    set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
  endif()
endif()

#
# Boost unit_test_framework and serialization are used throughout the software
# so that they are a requirement
#
if(SAFE_LINKING)
  message(STATUS "Option SAFE_LINKING is enabled. Make sure static boost "
                 "libraries are available and are compiled with -fPIC.")
  # Link boost statically
  set(Boost_USE_STATIC_LIBS ON)
else()
  add_definitions(-DBOOST_TEST_DYN_LINK=1)
endif()
find_package( Boost COMPONENTS
  unit_test_framework
  serialization
  program_options # Examples, Geneva
  filesystem # Geneva
  REQUIRED 
)

# Find python3
# unfortunately there is a change in the interface of CMake at version 12
# in which the PythonInterp gets deprecated and FindPython3 should be used
set(Python3_FOUND FALSE)
if (${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} LESS 12)
  find_package(PythonInterp)
  if(${PYTHONINTERP_FOUND} AND ${PYTHON_VERSION_MAJOR} GREATER 2)
    set(Python3_FOUND TRUE)
    set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
  endif()
else()
  set (CMAKE_FIND_FRAMEWORK NEVER)
  find_package(Python3)
endif()

find_package(ROOT COMPONENTS Minuit2)

find_package(GSL COMPONENTS gsl gslcblas REQUIRED)

# Geneva minimizer module (optional)
find_package(Geneva QUIET)

# Third party libraries included in the repository
# - easyloggingpp, pybind11, qft++, TBB, parallelSTL, EvtGen
add_subdirectory(ThirdParty)

# Setting ComPWA source dir as include directory
include_directories(${COMPWA_SOURCE_DIR})

#
# Enable target 'test'
#
enable_testing()
set(CTEST_OUTPUT_ON_FAILURE TRUE)

#
# Submodules
#
add_subdirectory(Core)
add_subdirectory(Tools)
add_subdirectory(Data)
add_subdirectory(Estimator)
add_subdirectory(Optimizer)
add_subdirectory(Physics)
add_subdirectory(Examples)
add_subdirectory(pycompwa)
