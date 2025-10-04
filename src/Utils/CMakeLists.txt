project(metasdkcpp_utilities)

cmake_minimum_required(VERSION 3.16)

# set(LIBRARY_OUTPUT_PATH ../)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories(
	include/
	eigen-3.4.0/
)

file(GLOB MOCS_HDRS  "include/*.h" "include/*.hpp")
file(GLOB SRC_LIST  "src/*.c" "src/*.cpp")

set(SRCS
	${SRC_LIST}
)


add_library(
    ${PROJECT_NAME} SHARED 
    ${SRCS}
)

target_link_libraries(${PROJECT_NAME}
	z
	pthread
)

install(TARGETS ${PROJECT_NAME} 
	ARCHIVE DESTINATION ./
	LIBRARY DESTINATION ${INSTALL_LIB_DIR}
)

