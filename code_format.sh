#!/bin/sh

CLANG_FORMAT=clang-format-9

find . -type f \( -iname "*.hpp" -or -iname "*.cpp" \) | xargs $CLANG_FORMAT -style=file -i
