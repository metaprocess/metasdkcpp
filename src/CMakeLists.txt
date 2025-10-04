
cmake_minimum_required(VERSION 3.16)

set(CMAKE_PROJECT_INCLUDE_BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/cmake/EnableCMP0048.cmake")

project(metasdkcpp)

set(PROJECT_VERSION "${sdk_cpp_VERSION}")
set(PROJECT_VERSION_MAJOR "${sdk_cpp_VERSION_MAJOR}")
set(PROJECT_VERSION_MINOR "${sdk_cpp_VERSION_MINOR}")
set(PROJECT_VERSION_PATCH "${sdk_cpp_VERSION_PATCH}")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Include OptimizeFlags.cmake only if not in Debug mode
string(TOLOWER "${CMAKE_BUILD_TYPE}" lowercase_cmake_build_type)
if(NOT lowercase_cmake_build_type STREQUAL "debug")
  include(cmake/OptimizeFlags.cmake)
else()
  include(cmake/DebugFlags.cmake)
endif()

set(CMAKE_VERBOSE_MAKEFILE OFF)

get_property(LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS)

set(CMAKE_INSTALL_PREFIX "/usr" )

if (${LIB64} STREQUAL "TRUE")
	set(LIBSUFFIX 64)
else()
	set(LIBSUFFIX "")
endif()
set(INSTALL_LIB_DIR  lib${LIBSUFFIX} CACHE PATH "Installation directory for libraries")



set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
  if(GIT_FOUND)
  message(STATUS "gitexecuteable=${GIT_EXECUTABLE}")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE "PRJ_BUILD_VERSION"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    message( STATUS "Git version: ${PRJ_BUILD_VERSION}" )
  else()
    set(PRJ_BUILD_VERSION 0)
  endif(GIT_FOUND)
endif()


set(CMAKE_PACKAGE_VERSION_MAJOR "1" CACHE STRING "major version")
set(CMAKE_PACKAGE_VERSION_MINOR "1" CACHE STRING "minor version")
set(CMAKE_PACKAGE_VERSION_PATCH "5" CACHE STRING "patch version")

set(CMAKE_PACKAGE_VERSION "${CMAKE_PACKAGE_VERSION_MAJOR}.${CMAKE_PACKAGE_VERSION_MINOR}.${CMAKE_PACKAGE_VERSION_PATCH}")

set(PROJECT_VERSION_MAJOR "${CMAKE_PACKAGE_VERSION_MAJOR}")
set(PROJECT_VERSION_MINOR "${CMAKE_PACKAGE_VERSION_MINOR}")
set(PROJECT_VERSION_PATCH "${CMAKE_PACKAGE_VERSION_PATCH}")
# set(PROJECT_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")
set(PROJECT_VERSION 1.0.0)

message(STATUS "PACKAGE VERSION:: ${CMAKE_PACKAGE_VERSION}")

add_definitions(-DPRJ_BUILD_VER=\"${PRJ_BUILD_VERSION}\")
add_definitions(-DPRJ_PACKAGE_VER=\"${CMAKE_PACKAGE_VERSION}\")

option(CMAKE_RPM_MODE "When the app is installed in system" OFF)
get_property(LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS)

if (${LIB64} STREQUAL "TRUE")
    set(LIBSUFFIX 64)
else()
    set(LIBSUFFIX "")
endif()
set(INSTALL_LIB_DIR  lib${LIBSUFFIX} CACHE PATH "Installation directory for libraries")

# Include the modular resource embedding script
include(cmake/EmbedResources.cmake)

add_subdirectory(Utils)
add_subdirectory(core)


option(ENABLE_TESTING "Build tests and add resources for testing purposes" ON)

enable_testing()
add_subdirectory(tests)
# if(BUILD_TESTING)
# endif()