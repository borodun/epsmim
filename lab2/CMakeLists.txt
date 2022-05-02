cmake_minimum_required(VERSION 3.10)
project(lab2)

set(CMAKE_CXX_STANDARD 17)

function(setup name parts_str)
    set(dependencies "")
    string(REPLACE " " ";" parts ${parts_str})
    foreach (part ${parts})
        list(APPEND dependencies ${name}/${part})
    endforeach ()

    add_executable(${name} ${dependencies})
endfunction(setup)

setup(1-compiler-vectorization "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")
set_target_properties(1-compiler-vectorization PROPERTIES COMPILE_FLAGS "-O3 -g -ftree-vectorize -fopt-info-vec-optimized -fopt-info-vec-missed")

setup(2-ivdep "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")
set_target_properties(2-ivdep PROPERTIES COMPILE_FLAGS "-O3 -g -ftree-vectorize -fopt-info-vec-optimized -fopt-info-vec-missed")

setup(3-omp "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")
set_target_properties(3-omp PROPERTIES COMPILE_FLAGS "-O3 -g -ftree-vectorize -fopt-info-vec-optimized -fopt-info-vec-missed")

setup(4-avx2 "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")
set_target_properties(4-avx2 PROPERTIES COMPILE_FLAGS "-O3 -g -mavx2")

setup(5-avx512 "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")
set_target_properties(5-avx512 PROPERTIES COMPILE_FLAGS "-O3 -g -mavx512f")