cmake_minimum_required(VERSION 3.6)
project(inprot)

# Enable C++ 14 and define Cmake modules path
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

# Check the required compiler version is met
include(CheckCompilerVersion)

# Check the available CPU's ISAs and configure compiler flags
include(CheckCPUFeatures)
set(CMAKE_VERBOSE_MAKEFILE 1)

#
# Source Files
#
set(SOURCE_FILES
        # Main, fasta & k-mers
        main.cpp
        globals.h
        fasta_seq.h
        fasta_utils.h
        kmer_offset.h
        kmers_manager.h
        group_koff.h

        # Molecular descriptors
        md.h
        charge_scale.h
        hydrophilicity_scale.h
        hydrophobicity_scale.h
        pka_scale.h ras.h
        reduced_alphabets.h

        # CLI, svm
        CLI11.hpp
        cli.h
        rang.hpp
        svm_scaling.h)


# Set to "TRUE" if using LIBSVM => FALSE by default
set(USE_LIBSVM FALSE)

# Check if using the original libsvm or my own's implementation
if (${USE_LIBSVM})
    message(STATUS "Using LIBSVM")
    add_definitions(-DUSE_LIBSVM)
    set(SOURCE_FILES ${SOURCE_FILES} svm.cpp svm.h)
else()
    message(STATUS "Using own's libsvm")
    set(SOURCE_FILES ${SOURCE_FILES} libsvm.h svm_model.h)
endif()

add_executable(inprot ${SOURCE_FILES})


# If you have another Intel TBB installation, set the following three variables:
#set(TBB_ROOT "/usr/")
#set(TBB_LIBRARY "/usr/lib64")
#set(TBBINCLUDE_DIR "/usr/include")

# Find Intel TBB library
find_package(TBB)

# Should check if TBB found
message(STATUS "Found Intel TBB version: ${TBB_VERSION}")

# Only for debug purpose
if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    message(STATUS "Debug mode enabled")
    add_definitions(-DTBB_USE_DEBUG_BUILD)
    add_definitions(-DTBB_USE_DEBUG=1)
    add_definitions(-DDEBUG)
    add_definitions(-DTBB_USE_THREADING_TOOLS)
else()
    message(STATUS "Release mode enabled")
endif()

target_link_libraries(inprot tbb)
