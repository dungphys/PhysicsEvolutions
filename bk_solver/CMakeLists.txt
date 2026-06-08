cmake_minimum_required(VERSION 3.10)
project(bk_solver 
    VERSION 1.0.0
    LANGUAGES CXX
)

# ===== C++ STANDARD =====
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ===== COMPILER FLAGS =====
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -D_DEBUG")

# ===== GENERATE COMPILE COMMANDS (for VS Code IntelliSense) =====
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ===== Find GSL =====
find_package(GSL REQUIRED)
if(NOT GSL_FOUND)
    message(FATAL_ERROR "GSL not found! Install with: apt-get install libgsl-dev")
endif()
message(STATUS "GSL found: ${GSL_INCLUDE_DIRS}")

# ===== SOURCE FILES =====
set(SOURCES
    src/utils.cpp
    src/init_condition.cpp
    src/interp.cpp
    src/solver.cpp
    src/main.cpp
)

set(HEADERS
    include/utils.h
    include/init_condition.h
    include/interp.h
    include/solver.h
)

# ===== CREATE EXECUTABLE ===== 
add_executable(bk_solver ${SOURCES} ${HEADERS})
target_include_directories(bk_solver PRIVATE include ${GSL_INCLUDE_DIRS})
target_link_libraries(bk_solver PRIVATE GSL::gsl GSL::gslcblas m)

# ===== PRINT BUILD INFO =====
message(STATUS "")
message(STATUS "===== Build Configuration =====")
message(STATUS "Project: ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "CXX Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "CXX Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Source files: ${SOURCES}")
message(STATUS "Output: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bk_kl_solver")
message(STATUS "")