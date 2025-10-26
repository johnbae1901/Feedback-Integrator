cmake_minimum_required(VERSION 3.10)
project(euler_feedback VERSION 1.0.0)

# Require a C++17 compiler
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3 (make sure it's installed on your system)
find_package(Eigen3 REQUIRED)

# Manually configure OpenMP for macOS
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/Cellar/libomp/19.1.6/include")
set(OpenMP_EXE_LINKER_FLAGS "-L/opt/homebrew/Cellar/libomp/19.1.6/lib -lomp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# Confirm OpenMP configuration
message(STATUS "OpenMP CXX Flags: ${OpenMP_CXX_FLAGS}")
message(STATUS "OpenMP Linker Flags: ${OpenMP_EXE_LINKER_FLAGS}")

add_executable(main main.cpp
                            include/basic_operations.cpp
                            include/euler_feedback.cpp
                            include/euler_vanilla.cpp
                            include/euler_feedback_adaptive.cpp
                            include/stormer_verlet_B.cpp
                            include/euler_feedback_adaptive_light.cpp
                            include/getError.cpp
                            include/linspace.cpp
                            include/estimateL.cpp
                            include/csv_io.cpp)

target_link_libraries(main PRIVATE Eigen3::Eigen)

# Link libraries
target_link_libraries(main PRIVATE 
    ${OpenMP_EXE_LINKER_FLAGS} 
)
