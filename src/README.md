# MetaSDK C++ Development Guide

This SDK provides a robust C++ framework for scientific computing with efficient resource embedding, parallel processing capabilities, and a modular architecture.

## Table of Contents
- [Compilation and Testing](#compilation-and-testing)
- [Docker Development Environment](#docker-development-environment)
- [Resource Embedding System](#resource-embedding-system)
- [Adding New Process Modules](#adding-new-process-modules)
- [Thread Pool and Parallel Processing](#thread-pool-and-parallel-processing)

## Compilation and Testing

### Prerequisites
- CMake 3.16 or higher
- C++17 compatible compiler (g++/clang++)
- Git (optional, for version information)

### Building the Project
```bash
# Create build directory
mkdir build && cd build

# Configure with CMake
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build all targets
cmake --build .

# Build specific target (e.g., core executable)
cmake --build . --target metasdkcpp_core_exe
```

### Running Tests with CTest
The project includes comprehensive testing using Catch2 framework:

```bash
# Run all tests
ctest

# Run tests with verbose output
ctest -V

# Run specific test (e.g., lee_filter tests)
ctest -R lee_filter

# Build and run tests manually
./build/tests/metasdkcpp_tests
```

Test resources are automatically embedded into the test executable using the resource embedding system (see below).

## Docker Development Environment

### Dockerfile
The provided `Dockerfile` creates a lightweight Alpine Linux container with all necessary dependencies:

- **Base Image**: `alpine:3.20`
- **Build Tools**: cmake, g++, make, ninja, git
- **System Libraries**: zlib-dev, binutils, linux-headers, musl-dev
- **Build Process**: Automatically compiles the project during container build

### Docker Compose
The `docker-compose.yml` file provides easy container orchestration:

```bash
# Build and run the core executable
docker-compose up --build

# Run interactively for development
docker-compose run --rm metasdkcpp-core /bin/sh

# Run specific command
docker-compose run --rm metasdkcpp-core ./build/tests/metasdkcpp_tests
```

### Development Workflow
1. **Build Container**: `docker-compose build`
2. **Run Application**: `docker-compose up`
3. **Development**: Mount source directory and rebuild inside container
4. **Testing**: Run ctest inside the container environment

## Resource Embedding System

The SDK features a lightweight and efficient resource embedding system that compiles binary resources directly into executables, eliminating external file dependencies.

### How It Works
1. **Resource Declaration**: Resources are placed in `resources/` subdirectories
2. **Automatic Embedding**: The `embed_resources()` CMake function automatically processes all files in resource directories
3. **Runtime Access**: Resources are accessible via the `get_resource()` function with string keys

### Usage Example
```cpp
// Load embedded resource
auto resource = get_resource("src/lee_filter/resources/lee_filter_input_image.bin");
auto matrix = MatrixUtils::load_from_resource<Eigen::MatrixXcd>(resource);

// Load from file path (automatically detects embedded vs. filesystem)
auto matrix2 = MatrixUtils::load_from_file(":/src/lee_filter/resources/lee_filter_input_image.bin");
```

### Benefits
- **Zero Runtime Dependencies**: All resources are compiled into the binary
- **Performance**: Direct memory access without file I/O overhead
- **Portability**: Single executable deployment
- **Security**: Resources are embedded and cannot be tampered with externally

## Adding New Process Modules

The SDK follows a modular architecture where new processing modules can be easily added.

### Step 1: Create Module Directory Structure
```
core/
└── your_module/
    ├── your_module.cpp
    ├── your_module.hpp
    └── resources/          # Optional: embedded test data
```

### Step 2: Implement Module Header (`your_module.hpp`)
```cpp
#include <Eigen/Dense>
#include <complex>
#include "MatrixUtils.h"

// Your function declarations
Eigen::MatrixXcd YourFunction(const Eigen::MatrixXcd& input, int parameter);
```

### Step 3: Implement Module Source (`your_module.cpp`)
```cpp
#include "your_module.hpp"
#include "Utils.h"
#include "MatrixUtils.h"
#include <ElapsedTimer.h>
#include RESOURCES_HEADER  // Required for resource access

Eigen::MatrixXcd YourFunction(const Eigen::MatrixXcd& input, int parameter) {
    // Your implementation here
    // Use MatrixUtils for matrix operations
    // Access embedded resources if needed
    return processed_result;
}
```

### Step 4: Create Tests
Create corresponding test directory:
```
tests/src/
└── your_module/
    ├── test.cpp
    └── resources/          # Test input/output data
```

### Step 5: Implement Tests (`test.cpp`)
```cpp
#include <catch2/catch_test_macros.hpp>
#include "your_module.hpp"
#include "CommonConstants.h"
#include <ElapsedTimer.h>

MatrixCompareResult test_your_function() {
    // Load test input from embedded resources
    auto input = std::get<Eigen::MatrixXcd>(
        MatrixUtils::load_from_file(":/src/your_module/resources/input.bin")
    );
    
    // Execute function
    auto result = YourFunction(input, test_parameter);
    
    // Load expected output
    auto expected = std::get<Eigen::MatrixXcd>(
        MatrixUtils::load_from_file(":/src/your_module/resources/expected.bin")
    );
    
    // Compare results
    return MatrixUtils::compare_matrices(result, expected);
}

TEST_CASE("Your Function Test", "[your_module]") {
    auto result = test_your_function();
    REQUIRE(result.max_diff() <= CONST_ACCEPTABLE_MAX_DIFF_ERROR);
    REQUIRE(result.sum_diff() <= CONST_ACCEPTABLE_SUM_DIFF_ERROR);
}
```

### Step 6: CMake Integration
The build system automatically discovers new modules through `file(GLOB_RECURSE)` patterns in the CMakeLists.txt files, so no manual CMake modifications are typically needed.

## Thread Pool and Parallel Processing

The SDK includes a powerful `StdThreadPool` class for easy parallelization of computationally intensive tasks.

### Key Features
- **Automatic Thread Count**: Uses `std::thread::hardware_concurrency()` by default
- **Task Queue**: Thread-safe task submission and execution
- **Parallelize Helper**: Simplifies parallel execution of array operations

### Basic Usage
```cpp
#include "std_thread_pool.h"

// Create thread pool
StdThreadPool pool;

// Submit individual tasks
pool.push_task([]() {
    // Your task code here
});

// Wait for all tasks to complete
pool.wait_for_tasks();
```

### Parallelize Method
The `parallelize()` method is particularly convenient for parallelizing loops or array operations:

```cpp
// Parallelize work across multiple threads
pool.parallelize(
    num_elements,  // Total number of work items
    [](const TaskArgsGeneral& args) {
        // args.index_lower: start index for this thread
        // args.index_upper: end index for this thread  
        // args.index_job: thread job index
        // args.ptr_any: optional user data pointer
        
        for (int i = args.index_lower; i < args.index_upper; ++i) {
            // Process element i
            process_element(i);
        }
    },
    
);
```

### Benefits of Parallelize
- **Automatic Load Balancing**: Work is evenly distributed across available threads
- **Simple API**: No need to manage thread synchronization manually
- **Flexible**: Supports any work that can be divided into index ranges
- **Efficient**: Minimizes thread creation overhead through thread reuse

### Example: Parallel Matrix Processing
```cpp
StdThreadPool pool;
Eigen::MatrixXd matrix(rows, cols);

pool.parallelize(rows, [&](const TaskArgsGeneral& args) {
    for (int row = args.index_lower; row < args.index_upper; ++row) {
        for (int col = 0; col < cols; ++col) {
            matrix(row, col) = expensive_computation(row, col);
        }
    }
});
```

This approach makes it trivial to parallelize computationally intensive operations while maintaining clean, readable code.

