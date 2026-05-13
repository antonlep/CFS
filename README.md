# CFSolve

## Development Setup

> **All commands must be run from the x64 Native Tools Command Prompt.**

### Prerequisites

- Visual Studio 2022+ with C++ workload
- Python 3.12+
- [Ninja](https://ninja-build.org/)
- [LLVM/clang-cl](https://releases.llvm.org/)

### Initial Setup (once)

```cmd
python -m venv .venv
.\.venv\Scripts\activate
python -m pip install --upgrade pip setuptools wheel
pip install -e .
cmake -S . -B build -G Ninja ^
    -DCMAKE_CXX_COMPILER=clang-cl ^
    -DCMAKE_C_COMPILER=clang-cl ^
    -DCMAKE_BUILD_TYPE=Debug ^
    -DCFS_BUILD_TESTS=ON ^
    -DCMAKE_EXPORT_COMPILER_COMMANDS=ON
```

### GUI

```cmd
python src\gui\main.py
```

### CLI

```cmd
ninja -C build solver_cli
build\src\solver\solver_cli.exe examples\fem.inp
```

### Testing

```cmd
ninja -C build
ctest --test-dir build --output-on-failure
pytest tests -v
```

:: ── C++ formatting ──
clang-format --dry-run --Werror -i src/solver/*.cpp src/solver/*.hpp

:: ── C++ static analysis ──
clang-tidy src/solver/*.cpp -p build

:: ── Python formatting ──
black --check src/gui

:: ── Python linting ──
flake8 src/gui