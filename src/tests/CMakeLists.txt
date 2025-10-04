cmake_minimum_required(VERSION 3.16)
project(tests)
# Method 2: Submodule (alternative)
add_subdirectory(catch2)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    catch2/src
    ../core/include
    ../core/lee_filter
    ../Utils/eigen-3.4.0
    ../Utils/include
)

file(GLOB SRC_LIST "src/*/.c" "src/*/*.cpp")

add_executable(tests 
    ${SRC_LIST}
)

file(GLOB RCS_LIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/*/resources/*")

embed_resources(${PROJECT_NAME} ${RCS_LIST})

target_link_libraries(tests PRIVATE 
    Catch2::Catch2WithMain
    metasdkcpp_core
)

# Enable testing and add Catch2's test discovery
include(${CMAKE_CURRENT_SOURCE_DIR}/catch2/extras/Catch.cmake)
catch_discover_tests(tests)
