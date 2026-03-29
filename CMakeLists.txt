cmake_minimum_required(VERSION 3.15)
project(VOR-rPPG VERSION 1.0.0 LANGUAGES CXX)

# Standard C++20 for modern DSP
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library definition
add_library(vor_core
    src/VOREngine.cpp
    src/vor_math.cpp
    src/vor_extractors.cpp
)

# Include directories
target_include_directories(vor_core PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Compiler flags for optimization
if(MSVC)
    target_compile_options(vor_core PRIVATE /W4 /O2)
else()
    target_compile_options(vor_core PRIVATE -Wall -Wextra -O3 -ffast-math)
endif()

# Installation rules (for library users)
install(TARGETS vor_core
    EXPORT vor_core_targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(DIRECTORY include/ DESTINATION include)

# Export targets
install(EXPORT vor_core_targets
    FILE vor_core_targets.cmake
    NAMESPACE VOR::
    DESTINATION lib/cmake/vor_core
)

# Example executable (optional validation)
add_executable(vor_benchmark examples/benchmark.cpp)
target_link_libraries(vor_benchmark PRIVATE vor_core)
