cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(
    NHFINT
    VERSION 0.1
    LANGUAGES C CXX
)

# if (NOT CMAKE_BUILD_TYPE)
#     set(CMAKE_BUILD_TYPE Release)
# endif()

if (WIN32)
    add_definitions(-DNOMINMAX)
endif()

add_subdirectory(extern/nhfmath)
add_subdirectory(src)

add_subdirectory(extern/GoogleTest)
add_subdirectory(test)

