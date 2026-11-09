cmake_minimum_required(VERSION 3.10)

# 專案名稱
project(BedTreeProject)

# 指定使用 C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 告訴編譯器去哪裡找標頭檔 (.h)
include_directories(include)

# 將 src 底下的 cpp 檔編譯成一個叫做 bed_tree_run 的執行檔
add_executable(bed_tree_run 
    src/main.cpp 
    src/StringMapper.cpp
)