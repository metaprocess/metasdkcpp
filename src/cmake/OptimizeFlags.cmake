# Check if CMAKE_CXX_FLAGS_RELEASE is defined
if(NOT DEFINED CMAKE_CXX_FLAGS_RELEASE)
  set(CMAKE_CXX_FLAGS_RELEASE "")
endif()

include(cmake/OptimizeMethods.cmake)
# Function to replace or add flags

# Apply the function to CMAKE_CXX_FLAGS_RELEASE
replace_or_add_flag("${CMAKE_CXX_FLAGS_RELEASE}" "-O2" "-O3" CMAKE_CXX_FLAGS_RELEASE)
add_string_if_not_exists(CMAKE_CXX_FLAGS_RELEASE "-mtune=native")
add_string_if_not_exists(CMAKE_CXX_FLAGS_RELEASE "-march=native")
add_string_if_not_exists(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")

replace_or_add_flag("${CMAKE_C_FLAGS_RELEASE}" "-O2" "-O3" CMAKE_C_FLAGS_RELEASE)
add_string_if_not_exists(CMAKE_C_FLAGS_RELEASE "-mtune=native")
add_string_if_not_exists(CMAKE_C_FLAGS_RELEASE "-march=native")
add_string_if_not_exists(CMAKE_C_FLAGS_RELEASE "-DNDEBUG")
# Output the result (optional, for debugging)
message(STATUS "Modified CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "Modified CMAKE_C_FLAGS_RELEASE: ${CMAKE_C_FLAGS_RELEASE}")
