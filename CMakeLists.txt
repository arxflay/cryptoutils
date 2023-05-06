cmake_minimum_required(VERSION 3.25)

project (cryptoutils VERSION 1.0 LANGUAGES CXX C)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindFmt.cmake)

if (NOT ${LIBFMT_FOUND})
    message(STATUS "Using LIBFMT submodule instead")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/third_party/fmt)
endif()

list(APPEND SOURCES
    src/algos.cpp
    src/generators.cpp
    src/ec.cpp
    src/rsa.cpp
    src/shamir.cpp
    src/elgamal.cpp
    src/handlers.cpp
    src/main.cpp
)

set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")

add_executable(cryptoutils ${SOURCES})
target_include_directories(cryptoutils PUBLIC include)

if (WIN32)
    target_link_libraries(cryptoutils PRIVATE fmt::fmt)
else()
    target_link_libraries(cryptoutils PRIVATE m fmt::fmt)
endif()

target_compile_options(cryptoutils PRIVATE "$<${gcc_like_cxx}:-Wall;-Wextra;-Wconversion;-std=c++17>"  
                                                        "$<${msvc_cxx}:/W3;/std:c++17>")
