FROM alpine:3.20

# Install build dependencies
RUN apk add --no-cache \
    cmake \
    git \
    g++ \
    make \
    ninja \
    binutils \
    zlib-dev \
    linux-headers \
    musl-dev \
    bash

# Set working directory
WORKDIR /app

# Copy source code
# COPY . .

# # Build the project
# RUN mkdir -p build && \
#     cd build && \
#     cmake -DCMAKE_BUILD_TYPE=Release -G Ninja .. && \
#     cmake --build . --target metasdkcpp_core_exe

# Create a minimal runtime image (multi-stage build optional)
# For now, we'll keep it in the same image for simplicity
# The executable will be in the root directory (CMAKE_RUNTIME_OUTPUT_DIRECTORY ../)

# CMD ["./metasdkcpp_core_exe"]