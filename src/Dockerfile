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
    bash \
    hdf5-dev \
    sudo \
    fftw-dev \
    python3 \
    python3-dev \
    py3-pip \
    gdb

# Create user with UID 1000
RUN addgroup -g 1000 devuser && \
    adduser -D -u 1000 -G devuser devuser && \
    echo "devuser ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/devuser && \
    chmod 0440 /etc/sudoers.d/devuser

# Set working directory
WORKDIR /app

# Change ownership of /app to devuser
RUN chown -R devuser:devuser /app

# Switch to non-root user
USER devuser

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
