# ==========================================
#  SSE for the 1D FM-TFIM (OBC)
#  Author: Yi-Ming Ding
#  Email: dingyiming@westlake.edu.cn
#  Last updated: April 22, 2024
# ==========================================
cmake_minimum_required(VERSION 3.17)
project(TFIM_1D)
set(CMAKE_CXX_STANDARD 17)
add_compile_options(-O3)
set(THE_SOURCE "./src/main.cpp")
set(OUTPUT_EXE "tfim")
add_executable(${OUTPUT_EXE} ${THE_SOURCE})