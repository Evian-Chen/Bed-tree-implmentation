#!/usr/bin/env bash
# 直接用 g++ 編譯所有工具 (不需要 make / cmake)。輸出在 build/ 。
# 用法: bash scripts/build.sh        (在專案根目錄；Git Bash / MSYS / Linux 皆可)
set -e
cd "$(dirname "$0")/.."
mkdir -p build
CXX=${CXX:-g++}
FLAGS="-std=c++14 -O2 -Wall -Wextra -Iinclude"
# MinGW 靜態連結，避免執行時找不到 libgcc/libstdc++ dll
case "$(uname -s)" in MINGW*|MSYS*|CYGWIN*) FLAGS="$FLAGS -static";; esac

TOOLS="ed_tool gram_tool zorder_tool bptree_tool bedtree_cli gen_data bench"
for t in $TOOLS; do
  echo "building $t"
  $CXX $FLAGS src/tools/$t.cpp src/StringMapper.cpp -o build/$t
done
echo "building test_all"
$CXX $FLAGS tests/test_all.cpp src/StringMapper.cpp -o build/test_all
echo "done -> build/"
