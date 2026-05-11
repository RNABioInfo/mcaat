#!/bin/bash

set -e

INSTALL=false
CLEAN=false

# Parse arguments
for arg in "$@"; do
    case $arg in
        --install)
            INSTALL=true
            shift
            ;;
        --clean)
            CLEAN=true
            shift
            ;;
        *)
            echo "Unknown option: $arg"
            echo "Usage: ./install.sh [--install] [--clean]"
            exit 1
            ;;
    esac
done

if $CLEAN; then
    echo "Cleaning build directory..."
    rm -rf build
    echo "Clean complete."
    exit 0
fi

echo "Installing system dependencies..."
sudo apt-get update
sudo apt-get install -y \
  build-essential \
  cmake \
  git \
  gettext \
  zlib1g-dev \
  libomp-dev \
  libpthread-stubs0-dev

echo "Initializing submodules..."
git submodule update --init --recursive

echo "Creating build directory..."
mkdir -p build
cd build

echo "Running CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

echo "Building..."
make -j$(nproc)

echo "Build complete. Binary: build/mcaat"

if $INSTALL; then
    echo "Installing binary to /usr/local/bin..."
    sudo cp mcaat /usr/local/bin/
    echo "Installed. Run 'mcaat' from anywhere."
fi
