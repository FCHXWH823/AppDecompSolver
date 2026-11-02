cmake_minimum_required(VERSION 3.16)
project(ising_decomposition_solver LANGUAGES CXX)

option(BUILD_DALTA_BASELINE
       "Build the vendored DALTA comparison executable (requires Gurobi and Boost)"
       OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(ising-decomp src/main.cpp)
target_compile_features(ising-decomp PRIVATE cxx_std_20)
target_compile_options(ising-decomp PRIVATE -O3 -Wall -Wextra)
target_link_libraries(ising-decomp PRIVATE Eigen3::Eigen)

set(GUROBI_HOME "$ENV{GUROBI_HOME}" CACHE PATH "Optional Gurobi installation root")
if(NOT GUROBI_HOME AND APPLE AND EXISTS "/Library/gurobi1202/macos_universal2")
    set(GUROBI_HOME "/Library/gurobi1202/macos_universal2")
endif()
find_path(GUROBI_INCLUDE_DIR gurobi_c++.h HINTS "${GUROBI_HOME}/include")
find_library(GUROBI_CXX_LIBRARY NAMES gurobi_c++ HINTS "${GUROBI_HOME}/lib")
find_library(GUROBI_LIBRARY NAMES gurobi130 gurobi120 gurobi110 gurobi100 gurobi95 gurobi91
             HINTS "${GUROBI_HOME}/lib")
if(GUROBI_INCLUDE_DIR AND GUROBI_CXX_LIBRARY AND GUROBI_LIBRARY)
    target_compile_definitions(ising-decomp PRIVATE HAVE_GUROBI=1)
    target_include_directories(ising-decomp PRIVATE ${GUROBI_INCLUDE_DIR})
    target_link_libraries(ising-decomp PRIVATE ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
endif()
if(BUILD_DALTA_BASELINE)
    if(NOT GUROBI_INCLUDE_DIR OR NOT GUROBI_CXX_LIBRARY OR NOT GUROBI_LIBRARY)
        message(FATAL_ERROR
                "BUILD_DALTA_BASELINE=ON requires Gurobi; set GUROBI_HOME")
    endif()

    find_package(Boost REQUIRED)
    find_package(Threads REQUIRED)
    add_executable(dalta-baseline third_party/dalta-baseline/src/main.cc)
    target_compile_features(dalta-baseline PRIVATE cxx_std_20)
    target_compile_options(dalta-baseline PRIVATE -O3 -Wall -Wextra)
    target_include_directories(dalta-baseline PRIVATE
        ${GUROBI_INCLUDE_DIR}
        ${Boost_INCLUDE_DIRS}
        third_party/dalta-baseline/src)
    target_link_libraries(dalta-baseline PRIVATE
        ${GUROBI_CXX_LIBRARY}
        ${GUROBI_LIBRARY}
        Threads::Threads
        m)
endif()
