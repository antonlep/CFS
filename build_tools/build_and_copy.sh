#!/bin/bash
# Portable script to build C++ solver and copy to Windows GUI folder
# Requires: CMake, pybind11, Eigen installed in WSL

# Load user config
CONFIG_FILE="build_tools/config.sh"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ Config file $CONFIG_FILE not found!"
    echo "Please copy config.sh.example → config.sh and edit the WINDOWS_GUI_DIR path"
    exit 1
fi
source "$CONFIG_FILE"

# Directories
SOLVER_DIR="solver"
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
if [ -z "$WINDOWS_GUI_DIR" ]; then
    echo "❌ WINDOWS_GUI_DIR is not set in config.sh!"
    exit 1
fi

echo "=== Copying .so to Windows GUI folder ==="
cp "$SO_FILE" "$WINDOWS_GUI_DIR/solver.so" || { echo "❌ Copy failed"; exit 1; }

echo "✔ Successfully copied solver.so → $WINDOWS_GUI_DIR"