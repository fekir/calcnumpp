cmake_minimum_required(VERSION 3.3)

if( CMAKE_SYSTEM_NAME STREQUAL "Windows" )
        include(windows_flags)
endif()

if     ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        include(msvc_flags)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        include(gcc_flags)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        include(clang_flags)
endif()
