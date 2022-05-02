cmake_minimum_required(VERSION 3.10)
project(lab3)

set(CMAKE_CXX_STANDARD 17)

function(setup name parts_str)
    set(dependencies "")
    string(REPLACE " " ";" parts ${parts_str})
    foreach (part ${parts})
        list(APPEND dependencies ${name}/${part})
    endforeach ()

    add_executable(${name} ${dependencies})
endfunction(setup)

setup(1-avx2 "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")
set_target_properties(1-avx2 PROPERTIES COMPILE_FLAGS "-O3 -g -mavx2")

setup(2-avx512 "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")
set_target_properties(2-avx512 PROPERTIES COMPILE_FLAGS "-O3 -g -mavx512f")