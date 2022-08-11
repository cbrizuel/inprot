#
# Created by germelcar on 10/3/17.
#

# -----------------------------------------------------------------------------------------
# --- Cmake module for checking the required compiler version for fully C++ 14 support. ---
# -----------------------------------------------------------------------------------------
# ** Only tested for GCC 7.2.0 (GNU, on GNU/Linux systems), Clang 5.0 (GNU and Mac systems)
#    and Intel 18.0 (GNU/Linux systems) compilers.
#
#
# Required versions:
# g++ >= 5.0
# clang >= 3.4
# intel c++ >= 17.0
#
#
# NOTE: Change the required version of the compiler depending on the requirements of your source code
#
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

    if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 5.0)
        message(FATAL_ERROR "GCC VERSION SHOULD BE >= 5.0")
    else()
        message(STATUS "Found GCC COMPILER ${CMAKE_CXX_COMPILER_VERSION}")
    endif()

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")

    if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 17.0)
        message(FATAL_ERROR "INTEL COMPILER SHOULD BE >= 17.0")
    else()
        message(STATUS "Found INTEL C++ COMPILER ${CMAKE_CXX_COMPILER_VERSION}")
    endif()

elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")

    if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 3.4)
        message(FATAL_ERROR "CLANG COMPILER SHOULD BE >= 3.4")
    else()
        message(STATUS "Found CLANG COMPILER ${CMAKE_CXX_COMPILER_VERSION}")
    endif()

endif()