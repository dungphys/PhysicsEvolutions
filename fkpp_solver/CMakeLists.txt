cmake_minimum_required(VERSION 3.10)
project(FKPP_Solver CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O2")

# Source files
set(SOURCES
    src/main.cpp
    src/solver.cpp
    src/init.cpp
    src/utils.cpp
)

# Header files
set(HEADERS
    include/solver.h
    include/init.h
    include/utils.h
)

# Create executable
add_executable(fkpp_solver ${SOURCES} ${HEADERS})

# Include directories
target_include_directories(fkpp_solver PRIVATE include)

# Optional: Enable testing
enable_testing()