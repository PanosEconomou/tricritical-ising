cmake_minimum_required(VERSION 3.23)

set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(fuse C CXX)
add_executable(fuse)

# Link external libraries
set_target_properties(fuse PROPERTIES LINKER_LANGUAGE CXX)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint)


find_path(FTI2_INCLUDE_DIR NAMES 4ti2/4ti2.h REQUIRED)
find_library(FTI2_GMP      NAMES 4ti2gmp    REQUIRED)
find_library(FTI2_COMMON   NAMES 4ti2common REQUIRED)
find_library(FTI2_UTIL     NAMES 4ti2util   REQUIRED)
find_library(FTI2_ZSOLVE   NAMES zsolve     REQUIRED)

add_library(4ti2 INTERFACE IMPORTED)
target_include_directories(4ti2 INTERFACE ${FTI2_INCLUDE_DIR})
find_library(GMPXX NAMES gmpxx REQUIRED)
target_link_libraries(4ti2 INTERFACE
    ${FTI2_GMP} ${FTI2_ZSOLVE} ${FTI2_COMMON} ${FTI2_UTIL}
    ${GMPXX}
)

target_link_libraries(fuse 
    PRIVATE
        PkgConfig::FLINT
        4ti2
)

# Source files
target_sources(fuse
    PRIVATE
        src/main.c
        src/interface.c
        src/solver.c
    PRIVATE
        FILE_SET HEADERS
        BASE_DIRS
            include
        FILES
            include/interface.h
            include/solver.h
)
