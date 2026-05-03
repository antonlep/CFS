# ══════════════════════════════════════════
# ALL from x64 Native Tools Command Prompt
# ══════════════════════════════════════════

# ── Setup (once) ──
python -m venv .venv
.\.venv\Scripts\activate
pip install --upgrade pip setuptools wheel
pip install -e .
cmake -S . -B build -G Ninja ^
    -DCMAKE_CXX_COMPILER=clang-cl ^
    -DCMAKE_C_COMPILER=clang-cl ^
    -DCMAKE_BUILD_TYPE=Debug ^
    -DCFS_BUILD_TESTS=ON

# ── GUI ──
python src\gui\main.py

# ── Build & run CLI ──
ninja -C build solver_cli
build\src\solver\solver_cli.exe examples\fem.inp

# ── Build & run tests ──
ninja -C build
ctest --test-dir build --output-on-failure
pytest tests -v