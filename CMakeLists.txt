cmake_minimum_required(VERSION 3.12)
project(mcaat)

# Enable verbose makefile output
set(CMAKE_VERBOSE_MAKEFILE OFF)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Add optimization flags for performance
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

find_package(ZLIB)
if(NOT ZLIB_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        zlib
        GIT_REPOSITORY https://github.com/madler/zlib.git
        GIT_TAG v1.3.1
    )
    FetchContent_MakeAvailable(zlib)
    set(ZLIB_LIBRARIES zlibstatic)
    set(ZLIB_INCLUDE_DIRS ${zlib_SOURCE_DIR} ${zlib_BINARY_DIR})
endif()

# Add parallel hashmap for better performance
include(FetchContent)
FetchContent_Declare(
    phmap
    GIT_REPOSITORY https://github.com/greg7mdp/parallel-hashmap.git
    GIT_TAG v1.3.11
)
FetchContent_MakeAvailable(phmap)

find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(BZip2 REQUIRED)

# Options
option(COVERAGE "Generate coverage report" ON)
option(STATIC_BUILD "Build static execution" OFF)
option(SANITIZER "Enable sanitizers" OFF)
option(TSAN "Enable thread sanitizers" OFF)
option(DEBUG "Enable development features (JSON, debug, etc)" OFF)

# Only define DEBUG if ENABLE_DEBUG is ON - for release, no macro (uses #else in main.cpp)
option(ENABLE_DEBUG "Enable DEBUG macro for main function" OFF)
if(ENABLE_DEBUG)
    add_definitions(-DDEBUG)
endif()

# Include directories
include_directories(include)
include_directories(tests)

add_subdirectory(libs/cft)
include_directories(libs/cft/src)

add_subdirectory(libs/googletest)
include_directories(libs/googletest/googletest/include)

add_subdirectory(libs/kseqpp)
include_directories(libs/kseqpp/include)

include_directories(libs/megahit/src)

add_subdirectory(libs/rapidfuzz-cpp)
include_directories(libs/rapidfuzz-cpp)
include_directories(${phmap_SOURCE_DIR})

include_directories(libs/spoa)
include_directories(libs/spoa/include)

# Configure spoa options (disable exe and tests to build only the library)
set(spoa_build_exe OFF CACHE BOOL "Disable spoa executable build" FORCE)
set(spoa_build_tests OFF CACHE BOOL "Disable spoa tests build" FORCE)
set(spoa_use_cereal OFF CACHE BOOL "Disable cereal dependency" FORCE)  # Optional: Skip serialization if not needed
set(spoa_use_simde OFF CACHE BOOL "Disable SIMDe for simplicity" FORCE)  # Optional: Skip SIMD if not needed

# Add spoa subdirectory (builds the library)
add_subdirectory(libs/spoa)

# Include fmt
include(FetchContent)
FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 10.2.1
)
FetchContent_MakeAvailable(fmt)

# Source files
file(GLOB_RECURSE MAIN_SOURCE "src/*.cpp")
file(GLOB_RECURSE SOURCES_WITHOUT_MAIN "src/*.cpp")
list(REMOVE_ITEM SOURCES_WITHOUT_MAIN "${CMAKE_SOURCE_DIR}/src/main.cpp")
file(GLOB_RECURSE INCLUDES "include/*.h")
file(GLOB_RECURSE TESTS_SOURCE "tests/*.cpp")
list(APPEND MAIN_SOURCE src/isolate_protospacers.cpp)
file(GLOB_RECURSE ASMBL_SOURCE "libs/megahit/src/assembly/*.cpp")
file(GLOB_RECURSE LCASM_SOURCE "libs/megahit/src/localasm/*.cpp")
file(GLOB_RECURSE IDBA_SOURCE "libs/megahit/src/idba/*.cpp")
file(GLOB_RECURSE SDBG_SOURCE "libs/megahit/src/sdbg/*.cpp")
file(GLOB_RECURSE CX1_SOURCE "libs/megahit/src/sorting/*.cpp")
file(GLOB_RECURSE SEQ_SOURCE "libs/megahit/src/sequence/*.cpp")
file(GLOB_RECURSE TOOLKIT_SOURCE "libs/megahit/src/tools/*.cpp")
file(GLOB_RECURSE PROGRESSBAR_SOURCE "libs/progressbar/include/*.hpp")


# Additional source files
list(APPEND OTHER_SOURCE
    libs/megahit/src/main_assemble.cpp
    libs/megahit/src/main_buildlib.cpp
    libs/megahit/src/main_iterate.cpp
    libs/megahit/src/main_local_assemble.cpp
    libs/megahit/src/main_sdbg_build.cpp
    libs/megahit/src/utils/options_description.cpp
)

list(APPEND SOME_OTHERS libs/megahit/src/utils/options_description.cpp)

# Static build configuration
if (STATIC_BUILD)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif (STATIC_BUILD)

# Testing
enable_testing()
set(TEST_IMPL_SOURCES
    ${SOURCES_WITHOUT_MAIN}
    ${ASMBL_SOURCE}
    ${IDBA_SOURCE}
    ${SDBG_SOURCE}
    ${LCASM_SOURCE}
    ${SEQ_SOURCE}
    ${CX1_SOURCE}
    ${TOOLKIT_SOURCE}
    ${PROGRESSBAR_SOURCE}
    ${SOME_OTHERS}
)

# Compiler and linker flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXXH_INLINE_ALL -ftemplate-depth=3000")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprefetch-loop-arrays -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__XROOT__='\"${CMAKE_SOURCE_DIR}/src\"'")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ZLIB_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# Release and Debug flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O1 -D_LIBCPP_DEBUG -D_GLIBCXX_DEBUG")

# Add executable
add_executable(mcaat
    ${MAIN_SOURCE}
    ${ASMBL_SOURCE}
    ${IDBA_SOURCE}
    ${SDBG_SOURCE}
    ${LCASM_SOURCE}
    ${SEQ_SOURCE}
    ${CX1_SOURCE}
    ${TOOLKIT_SOURCE}
    ${SOME_OTHERS}
)

# Static build linker flags
if (STATIC_BUILD)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -static")
    set_target_properties(mcaat PROPERTIES LINK_SEARCH_START_STATIC ON CXX_STANDARD 17)
endif (STATIC_BUILD)

# Link libraries
target_link_libraries(mcaat PRIVATE ${ZLIB_LIBRARIES} rapidfuzz::rapidfuzz spoa::spoa fmt::fmt)


# Threading configuration
if(THREADS_HAVE_PTHREAD_ARG)
    target_compile_options(mcaat PUBLIC "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
    target_link_libraries(mcaat "${CMAKE_THREAD_LIBS_INIT}")
endif()

if(DEBUG)
    include(FetchContent)
    FetchContent_Declare(
        json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
    )
    FetchContent_MakeAvailable(json)
    target_link_libraries(mcaat PRIVATE nlohmann_json::nlohmann_json)
    target_compile_definitions(mcaat PRIVATE DEBUG)
endif()

# Doxygen Documenation for cpp configurations
find_package(Doxygen)
if (DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/include/Doxyfile)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/include/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} COPYONLY)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
else()
    message(STATUS "Doxygen not found, 'doc' target will not be available.")
endif()

# Testing
# if (UNIT_TESTS)
add_executable(runTests ${TESTS_SOURCE} ${TEST_IMPL_SOURCES})

target_link_libraries(runTests
    gtest
    gtest_main
    ${ZLIB_LIBRARIES}
    pthread
    fmt::fmt
    spoa::spoa
)

include(GoogleTest)
gtest_discover_tests(runTests)
# endif()
