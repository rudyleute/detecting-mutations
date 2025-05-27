cmake_minimum_required(VERSION 3.10)
project(DetectingMutations LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

# Add your source files to create the executable
add_executable(DetectingMutations main.cpp FilesManipulator.cpp Comparator.cpp)

# Include directories
include_directories(${HTSLIB_INCLUDE_DIRS})
include_directories("/usr/local/include/bamtools")
include_directories(${ZLIB_INCLUDE_DIRS})
target_include_directories(DetectingMutations PRIVATE ${Boost_INCLUDE_DIRS})

# Link directories (usually not necessary if system libs are found properly)
link_directories(${HTSLIB_LIBRARY_DIRS})
link_directories("/usr/local/lib")

# Find and Include Boost (header-only for Interval Tree)
find_package(Boost REQUIRED)

# Find zlib (required by BamTools)
find_package(ZLIB REQUIRED)

# Link libraries after executable is defined
target_link_libraries(DetectingMutations
        ${HTSLIB_LIBRARIES}
        bamtools
        ${ZLIB_LIBRARIES}
)
