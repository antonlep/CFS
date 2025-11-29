#!/bin/bash
# Portable script to build C++ solver and copy to GUI folder
# Requires: CMake, pybind11, Eigen installed in WSL

# Directories
SOLVER_DIR="solver"
GUI_DIR="gui"
BUILD_DIR="$SOLVER_DIR/build"

echo "=== Building C++ solver with CMake ==="

# Configure + build
cmake -S "$SOLVER_DIR" -B "$BUILD_DIR" || { echo "CMake configure failed"; exit 1; }
cmake --build "$BUILD_DIR" || { echo "Build failed"; exit 1; }

# Find generated .so module
SO_FILE=$(find "$BUILD_DIR" -name "*.so" | head -n 1)

if [ -z "$SO_FILE" ]; then
    echo "❌ No .so file found!"
    exit 1
fi

echo "Found solver module: $SO_FILE"

# Copy it into Windows Python GUI folder
if [ -z "$GUI_DIR" ]; then
    echo "❌ GUI_DIR is not set in config.sh!"
    exit 1
fi

echo "=== Copying .so to GUI folder ==="
cp "$SO_FILE" "$GUI_DIR/solver.so" || { echo "❌ Copy failed"; exit 1; }

echo "✔ Successfully copied solver.so → $GUI_DIR"