cmake_minimum_required(VERSION 3.29)
project(diff_evo)

set(CMAKE_CXX_STANDARD 23
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO)

set(MATPLOTLIBCPP_DIR lib/matplotlibcpp)
include_directories(${MATPLOTLIBCPP_DIR})

add_executable(diff_evo
        diff_evo.cpp
        main.cpp
        test_functions.cpp
        diff_evo.h
        structures.h
        test_functions.h
        diff_evo.tpp
)

find_package(Python3 COMPONENTS Development NumPy REQUIRED)

target_include_directories(${PROJECT_NAME} PRIVATE
        python_module
        ${Python3_INCLUDE_DIRS}
        ${Python3_NumPy_INCLUDE_DIRS}
)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE DEBUG)
else ()
    set( CMAKE_CXX_FLAGS  "-O3")
endif()

target_link_libraries(${PROJECT_NAME} PRIVATE ${Python3_LIBRARIES})
