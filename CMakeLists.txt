cmake_minimum_required(VERSION 3.10)
project(DetectingMutations LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add BamTools include directory
include_directories("/usr/local/include/bamtools")

# Find and Include Boost (header-only for Interval Tree)
find_package(Boost REQUIRED)

# Find zlib (required by BamTools)
find_package(ZLIB REQUIRED)

# Add the directory where BamTools and zlib libraries are located
link_directories("/usr/local/lib")

# Add your source files to create the executable
add_executable(DetectingMutations main.cpp FilesReader.cpp)

# Link BamTools and zlib libraries
target_link_libraries(DetectingMutations bamtools ${ZLIB_LIBRARIES})

# Include zlib headers
include_directories(${ZLIB_INCLUDE_DIRS})

# Include Boost headers (no linking needed for header-only libraries like Boost.Interval)
target_include_directories(DetectingMutations PRIVATE ${Boost_INCLUDE_DIRS})
