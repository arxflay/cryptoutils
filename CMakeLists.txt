cmake_minimum_required(VERSION 3.25)
include(cmake/ProjectVersion.cmake)
project (cryptoutils VERSION ${VER_MAJOR}.${VER_MINOR} LANGUAGES CXX C)
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
include(libs/ProjectLibs.cmake)

if (WIN32)
    target_link_libraries(cryptoutils PRIVATE fmt::fmt)
else()
    target_link_libraries(cryptoutils PRIVATE m fmt::fmt)
endif()

target_compile_options(cryptoutils PRIVATE "$<${gcc_like_cxx}:-Wall;-Wextra;-Wconversion;-std=c++17>"  
                                                        "$<${msvc_cxx}:/W3;/std:c++17>")
install(TARGETS cryptoutils RUNTIME DESTINATION bin)
include(cmake/ProjectInstaller.cmake)
