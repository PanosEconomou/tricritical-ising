cmake_minimum_required(VERSION 3.23)

set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(fuse C)
add_executable(fuse)

# Link external libraries
find_package(PkgConfig REQUIRED)
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint)

target_link_libraries(fuse 
    PRIVATE
        PkgConfig::FLINT
)

# Source files
target_sources(fuse
    PRIVATE
        src/main.c
        src/interface.c
    PRIVATE
        FILE_SET HEADERS
        BASE_DIRS
            include
        FILES
            include/interface.h
)
