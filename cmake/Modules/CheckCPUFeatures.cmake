#
# Created by germelcar on 10/3/17.
#

# ---------------------------------------------------------------------------------------------
# --- Cmake module for checking CPU ISAs available and passing C++ compiler flags according ---
# --- to the ISAs detecteds.                                                                ---
# ---------------------------------------------------------------------------------------------
# ** Only tested for GCC 7.2.0 (GNU, on GNU/Linux systems).
#
# The source code for testing the ISA availables in the CPU was obtained from MMseqs2:
# https://github.com/soedinglab/MMseqs2
#
# with the exception of the code for AVX-512, which was obtained from Intel's webpage:
# https://software.intel.com/articles/how-to-detect-knl-instruction-support
#
#
set(CXX_OLD_FLAGS "${CMAKE_CXX_FLAGS}")
set(C_OLD_FLAGS "${CMAKE_C_FLAGS}")

# ---------------
# --- AVX 512 ---
# ---------------
macro(CHECK_AVX512)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "-mavx512f")
    elseif(CMAKE_CXX_CMPILER_ID STREQUAL "Intel")
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            set(CMAKE_CXX_FLAGS "/QaxCOMMON-AVX512")
        else()
            set(CMAKE_CXX_FLAGS "-axCOMMON-AVX512")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        set(CMAKE_CXX_FLAGS "/arch:AVX512")
    endif()

    try_run(RUNRES_AVX512 COMPRES_AVX512
            "${CMAKE_BINARY_DIR}/check_avx512"
            "${CMAKE_SOURCE_DIR}/cmake/Modules/check_avx512.cpp"
            RUN_OUTPUT_VARIABLE RUNOUT_AVX512
            COMPILE_OUTPUT_VARIABLE COMPOUT_AVX512)

    if(RUNOUT_AVX512 STREQUAL "AVX512")
        set(HAVE_AVX512 TRUE)
    endif()
endmacro(CHECK_AVX512)


# ------------
# --- AVX2 ---
# ------------
macro(CHECK_AVX2)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "-mavx2")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            set(CMAKE_CXX_FLAGS "/QaxCORE-AVX2")
        else()
            set(CMAKE_CXX_FLAGS "-axCORE-AVX2")
        endif()
    elseif(CMAKE_CXX_COMPILER_DIR STREQUAL "MSVC")
        set(CMAKE_CXX_FLAGS "/arch:AVX2")
    endif()

    try_run(RUNRES_AVX2 COMPRES_AVX2
            "${CMAKE_BINARY_DIR}/check_avx2"
            "${CMAKE_SOURCE_DIR}/cmake/Modules/check_avx2.cpp"
            RUN_OUTPUT_VARIABLE RUNOUT_AVX2
            COMPILE_OUTPUT_VARIABLE COMPOUT_AVX2)

    if(RUNOUT_AVX2 STREQUAL "AVX2")
        set(HAVE_AVX2 TRUE)
    endif()
endmacro(CHECK_AVX2)


# -----------
# --- AVX ---
# -----------
macro(CHECK_AVX)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "-mavx")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            set(CMAKE_CXX_FLAGS "/QaxCORE-AVX")
        else()
            set(CMAKE_CXX_FLAGS "-axCORE-AVX")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        set(CMAKE_CXX_FLAGS "/arch:AVX")
    endif()

    try_run(RUNRES_AVX COMPRES_AVX
            "${CMAKE_BINARY_DIR}/check_avx"
            "${CMAKE_SOURCE_DIR}/cmake/Modules/check_avx.cpp"
            RUN_OUTPUT_VARIABLE RUNOUT_AVX
            COMPILE_OUTPUT_VARIABLE COMPOUT_AVX)

    if(RUNOUT_AVX STREQUAL "AVX")
        set(HAVE_AVX TRUE)
    endif()
endmacro(CHECK_AVX)


# ---------------
# --- SSE 4.2 ---
# ---------------
macro(CHECK_SSE42)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "-msse4.2")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            set(CMAKE_CXX_FLAGS "/QaxSSE4.2")
        else()
            set(CMAKE_CXX_FLAGS "-axSSE4.2")
        endif()
    endif()

    try_run(RUNRES_SSE42 COMPRES_SSE42
            "${CMAKE_BINARY_DIR}/check_sse42"
            "${CMAKE_SOURCE_DIR}/cmake/Modules/check_sse42.cpp"
            RUN_OUTPUT_VARIABLE RUNOUT_SSE42
            COMPILE_OUTPUT_VARIABLE COMPOUT_SSE42)

    if(RUNOUT_SSE42 STREQUAL "SSE42")
        set(HAVE_SSE42 TRUE)
    endif()
endmacro(CHECK_SSE42)


# ---------------
# --- SSE 4.1 ---
# ---------------
macro(CHECK_SSE41)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "-msse4.1")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            set(CMAKE_CXX_FLAGS "/QaxSSE4.1")
        else()
            set(CMAKE_CXX_FLAGS "-axSSE4.1")
        endif()
    endif()

    try_run(RUNRES_SSE41 COMPRES_SSE41
            "${CMAKE_BINARY_DIR}/check_sse41"
            "${CMAKE_SOURCE_DIR}/cmake/Modules/check_sse41.cpp"
            RUN_OUTPUT_VARIABLE RUNOUT_SSE41
            COMPILE_OUTPUT_VARIABLE COMPOUT_SSE41)

    if(RUNOUT_SSE41 STREQUAL "SSE41")
        set(HAVE_SSE41 TRUE)
    endif()
endmacro(CHECK_SSE41)


# ---------------
# --- SSE 3 ---
# ---------------
macro(CHECK_SSE3)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "-msse3")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            set(CMAKE_CXX_FLAGS "/QaxSSE3")
        else()
            set(CMAKE_CXX_FLAGS "-axSSE3")
        endif()
    endif()

    try_run(RUNRES_SSE3 COMPRES_SSE3
            "${CMAKE_BINARY_DIR}/check_sse3"
            "${CMAKE_SOURCE_DIR}/cmake/Modules/check_sse3.cpp"
            RUN_OUTPUT_VARIABLE RUNOUT_SSE3
            COMPILE_OUTPUT_VARIABLE COMPOUT_SSE3)

    if(RUNOUT_SSE3 STREQUAL "SSE3")
        set(HAVE_SSE3 TRUE)
    endif()
endmacro(CHECK_SSE3)


# -------------
# --- SSE 2 ---
# -------------
macro(CHECK_SSE2)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "-msse2")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            set(CMAKE_CXX_FLAGS "/QaxSSE2")
        else()
            set(CMAKE_CXX_FLAGS "-axSSE2")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        set(CMAKE_CXX_FLAGS "/arch:SSE2")
    endif()

    try_run(RUNRES_SSE2 COMPRES_SSE2
            "${CMAKE_BINARY_DIR}/check_sse2"
            "${CMAKE_SOURCE_DIR}/cmake/Modules/check_sse2.cpp"
            RUN_OUTPUT_VARIABLE RUNOUT_SSE2
            COMPILE_OUTPUT_VARIABLE COMPOUT_SSE2)

    if(RUNOUT_SSE2 STREQUAL "SSE2")
        set(HAVE_SSE2 TRUE)
    endif()
endmacro(CHECK_SSE2)


# -----------
# --- SSE ---
# -----------
macro(CHECK_SSE)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "-msse")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            set(CMAKE_CXX_FLAGS "/QaxSSE")
        else()
            set(CMAKE_CXX_FLAGS "-axSSE")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        set(CMAKE_CXX_FLAGS "/arch:SSE")
    endif()

    try_run(RUNRES_SSE COMPRES_SSE
            "${CMAKE_BINARY_DIR}/check_sse"
            "${CMAKE_SOURCE_DIR}/cmake/Modules/check_sse.cpp"
            RUN_OUTPUT_VARIABLE RUNOUT_SSE
            COMPILE_OUTPUT_VARIABLE COMPOUT_SSE)

    if(RUNOUT_SSE STREQUAL "SSE")
        set(HAVE_SSE TRUE)
    endif()
endmacro(CHECK_SSE)

# Run/execute the macros
CHECK_AVX512()
CHECK_AVX2()
CHECK_AVX()
CHECK_SSE42()
CHECK_SSE41()
CHECK_SSE3()
CHECK_SSE2()
CHECK_SSE()

set(CMAKE_CXX_FLAGS "${CXX_OLD_FLAGS}")
set(CMAKE_C_FLAGS "${C_OLD_FLAGS}")

# NOTE:
#
# In the case of Clang compiler, for Mac systems (mac sierra tested), if clang was installed
# using brew, the "CMAKE_CXX_COMPILER_ID" variable is equals to "AppleClang", that is the
# reason why the variable is "compared" with the keyword "MATCHES" instead of "STREQUAL".
#

if(HAVE_AVX512)
    message(STATUS "Found AVX512 extensions.")
elseif(HAVE_AVX2)
    message(STATUS "Found AVX2 extensions.")
elseif(HAVE_AVX)
    message(STATUS "Found AVX extensions.")
elseif(HAVE_SSE42)
    message(STATUS "Found SSE 4.2 extensions.")
elseif(HAVE_SSE41)
    message(STATUS "Found SSE 4.1 extensions.")
elseif(HAVE_SSE3)
    message(STATUS "Found SSE 3 extensions.")
elseif(HAVE_SSE2)
    message(STATUS "Found SSE 2 extensions.")
elseif(HAVE_SSE)
    message(STATUS "Found SSE extensions.")
endif()

#
# GNU (g++) compiler and Clang++ compiler
#
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")

    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        add_compile_options(-Wall -ftree-vectorize -march=native -mfpmath=sse -funroll-loops
                -fopt-info-vec-optimized -fdiagnostics-color=auto -mtune=generic)

    else()
        add_compile_options(-Weverything -fvectorize -fslp-vectorize -march=native -funroll-loops
                -Rpass=loop-vectorize -fcolor-diagnostics)
    endif()

    if(HAVE_AVX512)
        add_compile_options(-mavx512f)
    elseif(HAVE_AVX2)
        add_compile_options(-mavx2)
    elseif(HAVE_AVX)
        add_compile_options(-mavx)
    elseif(HAVE_SSE42)
        add_compile_options(-msse4.2)
    elseif(HAVE_SSE41)
        add_compile_options(-msse4.1)
    elseif(HAVE_SSE3)
        add_compile_options(-msse3)
    elseif(HAVE_SSE2)
        add_compile_options(-msse2)
    elseif(HAVE_SSE)
        add_compile_options(-msse)
    endif()

#
# Intel compiler
#
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        add_compile_options(/xHOST /Qipo /Wall)
    else()
        add_compile_options(-xHOST -ipo -Wall)
    endif()

    if(HAVE_AVX512)
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            add_compile_options(/QaxCOMMON-AVX512)
        else()
            add_compile_options(-axCOMMON-AVX512)
        endif()

    elseif(HAVE_AVX2)
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            add_compile_options(/QaxCORE-AVX2 /Qtune:core-avx2)
        else()
            add_compile_options(-axCORE-AVX2 -mtune=core-avx2 -march=core-avx2)
        endif()

    elseif(HAVE_AVX)
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /QaxCORE-AVX")
        else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -axCORE-AVX")
        endif()

    elseif(HAVE_SSE42)
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            add_compile_options(/QaxSSE4.2)
        else()
            add_compile_options(-axSSE4.2)
        endif()

    elseif(HAVE_SSE41)
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            add_compile_options(/QaxSSE4.1)
        else()
            add_compile_options(-axSSE4.1)
        endif()

    elseif(HAVE_SSE3)
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            add_compile_options(/QaxSSE3)
        else()
            add_compile_options(-axSSE3)
        endif()

    elseif(HAVE_SSE2)
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            add_compile_options(/QaxSSE2)
        else()
            add_compile_options(-axSSE2)
        endif()

    elseif(HAVE_SSE)
        if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            add_compile_options(/QaxSSE)
        else()
            add_compile_options(-axSSE)
        endif()

    endif()

#
# MSVC compiler
#
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    add_compile_options(/Wall /Qvec-report:1 /GL)

    if(HAVE_AVX512)
        add_compile_options(/arch:AVX512)
    elseif(HAVE_AVX2)
        add_compile_options(/arch:AVX2)
    elseif(HAVE_AVX)
        add_compile_options(/arch:AVX)
    elseif(HAVE_SSE2)
        add_compile_options(/arch:SSE2)
    elseif(HAVE_SSE)
        add_compile_options(/arch:SSE)
    endif()
endif()
