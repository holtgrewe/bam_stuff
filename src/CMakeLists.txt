# enable C++11 support
find_package (CXX11)

# write resulting binaries into "bin" directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# search SeqAn library
set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

# enable SeqAn dependencies
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS}")
add_definitions (${SEQAN_DEFINITIONS} -DSEQAN_ENABLE_DEBUG=0)
include_directories (${SEQAN_INCLUDE_DIRS})

#if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#    add_definitions("-DSEQAN_ENABLE_DEBUG=1")
#else ()
#    add_definitions("-DSEQAN_ENABLE_DEBUG=0")
#endif ()

# register our target
add_executable (bam_markdup
                bam_markdup.cpp
                bam_markdup_app.cpp
                bam_markdup_app.h
                bam_markdup_options.cpp
                bam_markdup_options.h
                progress_indicator.cpp
                progress_indicator.h)

target_link_libraries (bam_markdup ${SEQAN_LIBRARIES})
