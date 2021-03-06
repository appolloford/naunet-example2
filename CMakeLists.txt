# To use CUDA by enabling language
cmake_minimum_required(VERSION 3.18)
project(NAUNET VERSION 0.1.0 LANGUAGES CUDA CXX)

message(STATUS "CMAKE_SYSTEM_NAME          = ${CMAKE_SYSTEM_NAME}"         )
message(STATUS "CMAKE_CXX_COMPILER_ID      = ${CMAKE_CXX_COMPILER_ID}"     )
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_COMPILER_IS_GNUCXX   = ${CMAKE_COMPILER_IS_GNUCXX}"  )

include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_definitions(-DNAUNET_DEBUG)
endif()


set(SUNDIALS_MIN_VERSION "5.6.1")
set(SUNDIALS_DIR "" CACHE PATH "Path to SUNDIALSConfig.cmake")
find_package(SUNDIALS ${SUNDIALS_MIN_VERSION} REQUIRED CONFIG)

# find_package(CUDA REQUIRED)
# find_package(CUDAToolkit REQUIRED)

# if(CMAKE_BUILD_TYPE STREQUAL Debug)
#     cmake_print_properties(TARGETS SUNDIALS::cvode_static PROPERTIES
#                            LOCATION INTERFACE_INCLUDE_DIRECTORIES)
#     cmake_print_properties(TARGETS SUNDIALS::nvecserial_static PROPERTIES
#                            LOCATION INTERFACE_INCLUDE_DIRECTORIES)
# endif()


add_subdirectory(src)
add_subdirectory(test)
