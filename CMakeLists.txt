cmake_minimum_required(VERSION 3.10)
project(BedTreeProject CXX)

# 舊版 MinGW (g++ 6.3) 只完整支援 C++14，因此統一用 C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
include_directories(include)

add_library(bedcore STATIC src/StringMapper.cpp)

foreach(tool ed_tool gram_tool zorder_tool bptree_tool bedtree_cli gen_data bench)
  add_executable(${tool} src/tools/${tool}.cpp)
  target_link_libraries(${tool} bedcore)
endforeach()

add_executable(test_all tests/test_all.cpp)
target_link_libraries(test_all bedcore)
