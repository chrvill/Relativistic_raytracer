cmake_minimum_required(VERSION 3.10)
set(TARGET_NAME Relativistic-renderer)

project(${TARGET_NAME} VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

set(EIGEN3_INCLUDE_DIR "../../../Eigen")

find_package(OpenMP REQUIRED)

include_directories(include)
include_directories(${EIGEN3_INCLUDE_DIR})

add_definitions(-DEIGEN_NO_DEBUG)

file(GLOB SOURCES "src/*.cpp")
add_executable(${TARGET_NAME} ${SOURCES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(${TARGET_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

target_compile_options(${TARGET_NAME} PRIVATE -O3 -march=native)

include(CTest)
enable_testing()

