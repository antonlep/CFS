```
python -m venv .venv
.\.venv\Scripts\activate
python -m pip install --upgrade pip setuptools wheel
pip install -e .
python .\src\gui\main.py

x64 Native Tools Command Prompt for VS 2022
cmake -S . -B build-dev -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCFS_BUILD_TESTS=ON
cmake --build build-dev
ctest --test-dir build-dev --output-on-failure
pytest tests -v

local
ninja -C build-dev solver_cli
build-dev\src\solver\solver_cli.exe examples\fem.inp

CI
cmake -S . -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release -DCFS_BUILD_TESTS=ON
cmake --build build-release
ctest --test-dir build-release --output-on-failure
```
