cmake_minimum_required( VERSION 3.10 )
project( oktal LANGUAGES CXX )

# FIX 1: Find HDF5 package and define the HDF5::HDF5 target
find_package(HDF5 REQUIRED)

# *****************************************************************
# FIX 2: Enable policy CMP0079 (NEW) to allow subdirectories to link to
# targets (like 'oktal') defined in parent directories (the root).
cmake_policy(SET CMP0079 NEW)
# *****************************************************************

set( CMAKE_CXX_STANDARD 23 )
set( CMAKE_EXPORT_COMPILE_COMMANDS ON)
set( CMAKE_CXX_STANDARD_REQUIRED )
set( CMAKE_COMPILE_WARNING_AS_ERROR ON )
add_compile_options(-Wall -Wextra -Wconversion)

# Pull Helper Library first (needed for mdspan)

message( STATUS "Downloading AdvPT helper libraries..." )
include( FetchContent )
FetchContent_Declare(
    advpt-helpers
    URL https://teaching.pages.i10git.cs.fau.de/advpt-exercise-ws2526/assets/advpt-helpers.tar.gz
)
FetchContent_MakeAvailable( advpt-helpers )

# Create libraries

add_library( oktal )

target_include_directories( oktal PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${advpt-helpers_SOURCE_DIR}/include
    ${advpt-helpers_SOURCE_DIR}/extern/mdspan/include
)

target_link_libraries( oktal PUBLIC 
    HDF5::HDF5
    advpt::mdspan
)

add_subdirectory( src )
add_subdirectory( apps )

# Add test suite

option( OKTAL_BUILD_TESTSUITE "Build the Oktal test suite" ON )

if( OKTAL_BUILD_TESTSUITE )
    include( CTest )
    add_subdirectory( tests )
endif()