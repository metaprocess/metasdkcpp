project(metasdkcpp_core)

cmake_minimum_required(VERSION 3.16)

set(PROJECT_TARGET_EXE "${PROJECT_NAME}_exe")
set(PROJECT_TARGET_LIB "${PROJECT_NAME}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories(
	include/
	../Utils/eigen-3.4.0/
	../Utils/include/
)

file(GLOB_RECURSE MOCS_HDRS "*.h" "*.hpp")
file(GLOB_RECURSE SRC_LIST "*.c" "*.cpp")
file(GLOB RCS_LIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*/resources/*")
# file(GLOB RCS_LIST "*/resources/*")

message(STATUS "list resources files: ${RCS_LIST}")

# Define a macro with the project name
add_definitions(-DPROJECT_NAME="${PROJECT_TARGET_LIB}")

set(SRCS
	${SRC_LIST}
)

message(STATUS "core source files:\"${SRCS}\"")
add_executable(
    ${PROJECT_TARGET_EXE}
    ${SRCS}
)

add_library(
    ${PROJECT_TARGET_LIB} SHARED
    ${SRCS}
)

embed_resources(${PROJECT_TARGET_EXE} ${RCS_LIST})
embed_resources(${PROJECT_TARGET_LIB} ${RCS_LIST})

set(LIB_DEPS
	metasdkcpp_utilities
)

target_link_libraries(${PROJECT_TARGET_EXE}
	${LIB_DEPS}
	)
	
target_link_libraries(${PROJECT_TARGET_LIB}
	${LIB_DEPS}
)

install(TARGETS ${PROJECT_TARGET_EXE} 
	RUNTIME DESTINATION bin
)

install(TARGETS ${PROJECT_TARGET_LIB} 
	ARCHIVE DESTINATION ./
	LIBRARY DESTINATION ${INSTALL_LIB_DIR}
)
