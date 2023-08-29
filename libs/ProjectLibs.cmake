include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindFmt.cmake)

if (NOT ${LIBFMT_FOUND})
    set(CPM_LOCAL_PACKAGES_ONLY False)
    set(CPM_SOURCE_CACHE ${CMAKE_CURRENT_LIST_DIR}/third_party/CPM)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FetchCPM.cmake)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CPM/CPM.cmake)
    CPMAddPackage("gh:fmtlib/fmt#7.1.3")
endif()
