# Stage 1: Build
FROM debian:bookworm AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    gettext \
    zlib1g-dev \
    libbz2-dev \
    libomp-dev \
    libpthread-stubs0-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Clone your repo
RUN git clone --recurse-submodules https://github.com/RNABioInfo/mcaat.git .

RUN mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j$(nproc)

# Stage 2: Runtime
FROM debian:bookworm-slim

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libomp5 \
    zlib1g \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Copy binary from builder
COPY --from=builder /app/build/mcaat /usr/local/bin/mcaat

# Set default command
ENTRYPOINT ["mcaat"]