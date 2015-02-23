project (bam_markdup)
cmake_minimum_required (VERSION 2.8)

set (CMAKE_BUILD_TYPE RelWithDebInfo)

# Make the SeqAn CMake modules available.
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                       ${CMAKE_SOURCE_DIR}/ext/seqan/util/cmake)
# Tell CMake to look for headers in ext/seqan/include.
set (CMAKE_INCLUDE_PATH ${CMAKE_LIBRARY_PATH}
                        ${CMAKE_SOURCE_DIR}/ext/seqan/include)

# Set default build type.
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

add_subdirectory (src)
