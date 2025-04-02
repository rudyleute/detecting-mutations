cmake_minimum_required(VERSION 3.10)
project(DetectingMutations LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

# Add BamTools include directory
include_directories("/usr/local/include/bamtools")

# Add the zlib include directory (CMake will find zlib if it's installed)
find_package(ZLIB REQUIRED)

# Add the directory where BamTools and zlib libraries are located
link_directories("/usr/local/lib")

# Add your source files to create the executable
add_executable(DetectingMutations main.cpp SamReader.cpp)

# Link BamTools and zlib libraries
target_link_libraries(DetectingMutations bamtools ${ZLIB_LIBRARIES})

# Include zlib headers
include_directories(${ZLIB_INCLUDE_DIRS})
