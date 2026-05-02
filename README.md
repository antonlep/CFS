# ══════════════════════════════════════════
# ALL from x64 Native Tools Command Prompt
# ══════════════════════════════════════════

# ── Setup (once) ──
python -m venv .venv
.\.venv\Scripts\activate
pip install --upgrade pip setuptools wheel
pip install -e .

# ── GUI ──
python src\gui\main.py

# ── Dev build (MSVC, fast iteration) ──
cmake -S . -B build-dev -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCFS_BUILD_TESTS=ON
cmake --build build-dev
ctest --test-dir build-dev --output-on-failure

# ── CLI ──
ninja -C build-dev solver_cli
build-dev\src\solver\solver_cli.exe examples\fem.inp

# ── CI-parity build (before pushing) ──
cmake -S . -B build-ci -G Ninja ^
    -DCMAKE_CXX_COMPILER=clang-cl ^
    -DCMAKE_C_COMPILER=clang-cl ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DCFS_BUILD_TESTS=ON
cmake --build build-ci
ctest --test-dir build-ci --output-on-failure

# ── Python tests (uses pip install -e . build) ──
pytest tests -v
