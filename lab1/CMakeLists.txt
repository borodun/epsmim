cmake_minimum_required(VERSION 3.10)
project(lab1)

set(CMAKE_CXX_STANDARD 17)

function(setup name parts_str)
    set(dependencies "")
    string(REPLACE " " ";" parts ${parts_str})
    foreach (part ${parts})
        list(APPEND dependencies ${name}/${part})
    endforeach ()

    add_executable(${name} ${dependencies})
endfunction(setup)

setup(1-no_optimization "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")

setup(2-O1 "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")
set_target_properties(2-O1 PROPERTIES COMPILE_FLAGS "-O1 -g")

setup(3-O2 "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")
set_target_properties(3-O2 PROPERTIES COMPILE_FLAGS "-O2 -g")

setup(4-O3 "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")
set_target_properties(4-O3 PROPERTIES COMPILE_FLAGS "-O3 -g")

setup(5-operation-reorganize "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")
set_target_properties(5-operation-reorganize PROPERTIES COMPILE_FLAGS "-O3 -g")

setup(6-index-optimization "main.cpp waveparams.cpp waveparams.h wavesim.cpp wavesim.h exceptions.h")
set_target_properties(6-index-optimization PROPERTIES COMPILE_FLAGS "-O3 -g")