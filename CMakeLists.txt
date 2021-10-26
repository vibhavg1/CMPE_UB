cmake_minimum_required(VERSION 3.20)

project(CMPE_UB)
include_directories(/Users/vgogate/software/ortools/include include)
link_directories(/Users/vgogate/software/ortools/lib)

add_definitions(-DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Ofast")
set(CMAKE_CXX_STANDARD 14)


add_executable(CMPE_UB target/main.cpp src/MN.cpp src/myRandom.cpp
        src/CG.cpp src/MN_MCKP_Bound.cpp src/MCKP.cpp)
target_link_libraries(CMPE_UB protobuf glog gflags CbcSolver Cbc OsiCbc Cgl ClpSolver Clp OsiClp Osi CoinUtils ortools)