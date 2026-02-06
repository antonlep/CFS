```
python -m venv .venv
.\.venv\Scripts\activate
python -m pip install --upgrade pip setuptools wheel
pip install -e .
python .\src\gui\main.py

x64 Native Tools Command Prompt for VS 2022
cmake -S . -B build-dev -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCFS_BUILD_TESTS=ON
cmake --build build-dev
ctest --test-dir build-dev --verbose
pytest tests -v

CI
pip install -e . --config-settings=cmake.define.CFS_BUILD_TESTS=ON
ctest --test-dir _skbuild/*/build --verbose
pytest tests -v
```
