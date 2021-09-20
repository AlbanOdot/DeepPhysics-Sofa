cmake_minimum_required(VERSION 3.12)
project(DeepPhysicsSofa VERSION 1.0.0)

# Policies
cmake_policy(SET CMP0072 NEW)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0076 NEW)
cmake_policy(SET CMP0094 NEW)

set(DPS_VERSION ${PROJECT_VERSION})

## Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

## Change default install prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()

# Include CMake specific tools
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_LIST_DIR}/cmake")
list(INSERT CMAKE_PREFIX_PATH 0 "${CMAKE_BINARY_DIR}/cmake")
include(CheckCXXCompilerFlag)
include(CMakePackageConfigHelpers)
include(CMakeDependentOption)
include(DPSTools)
set(CMAKE_REQUIRED_QUIET 1)


CMAKE_DEPENDENT_OPTION(DPS_OPTIMIZE_FOR_NATIVE "Build with -march=native to allow msse support" ON "COMPILER_SUPPORTS_MARCH_NATIVE" OFF)

# SOFA options
find_package(SOFA COMPONENTS SofaFramework QUIET)
get_filename_component(SOFA_ROOT_DIR "${SofaFramework_DIR}/../../.." ABSOLUTE)
message(STATUS
        "DPS compiling with \n"
        "    SOFA version: ${SofaFramework_VERSION}\n"
        "    SOFA location: ${SOFA_ROOT_DIR}"
        )
list(APPEND CMAKE_PREFIX_PATH ${SOFA_MODULE_DIR})
list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)


#SP3 options
find_package(SofaPython3 QUIET)
if (SofaPython3_FOUND)
    set(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE})
    set(PYBIND11_PYTHON_VERSION ${PYTHON_VERSION})
    get_filename_component(SP3_ROOT_DIR "${SofaPython3_DIR}/../../.." ABSOLUTE)
    find_package(Python ${PYBIND11_PYTHON_VERSION} COMPONENTS Interpreter Development QUIET)
    find_package(pybind11 CONFIG QUIET)
    message(STATUS "DPS using python and pybind \n:
        Version: ${Python_VERSION}
        Executable: ${Python_EXECUTABLE}
        Headers: ${Python_INCLUDE_DIRS}
        Libraries: ${Python_LIBRARIES}
        pybind11 version: ${pybind11_VERSION}
        pybind11 config: ${pybind11_CONFIG}"
            )
else()
    message(STATUS "DPS without SofaPython (not found)")
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Compiler options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED 17)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
set(PYBIND11_FINDPYTHON ON)
set(PYBIND11_PYTHON_VERSION 3)
set(Python_FIND_FRAMEWORK LAST)
if(NOT WIN32)
    if ("${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}" MATCHES ".*[cC][lL][aA][nN][gG]")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
endif()

add_subdirectory(src)

install(
        FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/DPSTools.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindSOFA.cmake"
        DESTINATION lib/cmake/DPS
        COMPONENT headers
)
