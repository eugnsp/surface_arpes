#!/bin/sh

find src -type f \( -iname "*.hpp" -or -iname "*.cpp" \) | xargs clang-format -style=file -i
