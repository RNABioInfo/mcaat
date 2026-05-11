cmake_minimum_required(VERSION 3.12)
project(mcaat)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Dependencies ──────────────────────────────────────────────────────────────

find_package(ZLIB)
if(NOT ZLIB_FOUND)
    include(FetchContent)
    FetchContent_Declare(zlib
        GIT_REPOSITORY https://github.com/madler/zlib.git
        GIT_TAG v1.3.1)
    FetchContent_MakeAvailable(zlib)
    set(ZLIB_LIBRARIES zlibstatic)
    set(ZLIB_INCLUDE_DIRS ${zlib_SOURCE_DIR} ${zlib_BINARY_DIR})
endif()

include(FetchContent)

FetchContent_Declare(phmap
    GIT_REPOSITORY https://github.com/greg7mdp/parallel-hashmap.git
    GIT_TAG v1.3.11)
FetchContent_MakeAvailable(phmap)

FetchContent_Declare(fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 10.2.1)
FetchContent_MakeAvailable(fmt)

find_package(OpenMP REQUIRED)

# ── Options ───────────────────────────────────────────────────────────────────

option(STATIC_BUILD "Build static executable"      OFF)
option(SANITIZER    "Enable address/UB sanitizers" OFF)
option(ENABLE_DEBUG "Enable DEBUG macro (dev main)" OFF)

if(ENABLE_DEBUG)
    add_definitions(-DDEBUG)
    include(FetchContent)
    FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3)
    FetchContent_MakeAvailable(json)
endif()

# ── Compiler flags ────────────────────────────────────────────────────────────

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXXH_INLINE_ALL -ftemplate-depth=3000")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprefetch-loop-arrays -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -ggdb -O1 -D_GLIBCXX_DEBUG")

if(SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=undefined")
endif()

if(STATIC_BUILD)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

# ── Sub-libraries ─────────────────────────────────────────────────────────────

set(MEGAHIT_SRC "${CMAKE_SOURCE_DIR}/libs/megahit/src")

add_subdirectory(libs/cft)
add_subdirectory(libs/kseqpp)

set(spoa_build_exe   OFF CACHE BOOL "" FORCE)
set(spoa_build_tests OFF CACHE BOOL "" FORCE)
set(spoa_use_cereal  OFF CACHE BOOL "" FORCE)
set(spoa_use_simde   OFF CACHE BOOL "" FORCE)
add_subdirectory(libs/spoa)

# ── Include paths ─────────────────────────────────────────────────────────────

include_directories(include)
include_directories(${MEGAHIT_SRC})
include_directories(libs/cft/src)
include_directories(libs/kseqpp/include)
include_directories(libs/spoa/include)
include_directories(${phmap_SOURCE_DIR})

# ── Source globs ──────────────────────────────────────────────────────────────

file(GLOB_RECURSE MAIN_SOURCE          "src/*.cpp")
file(GLOB_RECURSE SOURCES_WITHOUT_MAIN "src/*.cpp")
list(REMOVE_ITEM SOURCES_WITHOUT_MAIN "${CMAKE_SOURCE_DIR}/src/main.cpp")

file(GLOB_RECURSE SDBG_SOURCE    "${MEGAHIT_SRC}/sdbg/*.cpp")
file(GLOB_RECURSE ASMBL_SOURCE   "${MEGAHIT_SRC}/assembly/*.cpp")
file(GLOB_RECURSE IDBA_SOURCE    "${MEGAHIT_SRC}/idba/*.cpp")
file(GLOB_RECURSE SEQ_SOURCE     "${MEGAHIT_SRC}/sequence/*.cpp")
file(GLOB_RECURSE CX1_SOURCE     "${MEGAHIT_SRC}/sorting/*.cpp")
file(GLOB_RECURSE TOOLKIT_SOURCE "${MEGAHIT_SRC}/tools/*.cpp")

set(MEGAHIT_EXTRAS
    ${MEGAHIT_SRC}/main_assemble.cpp
    ${MEGAHIT_SRC}/main_buildlib.cpp
    ${MEGAHIT_SRC}/main_iterate.cpp
    ${MEGAHIT_SRC}/main_sdbg_build.cpp
    ${MEGAHIT_SRC}/utils/options_description.cpp
)

set(MEGAHIT_SOURCES
    ${SDBG_SOURCE}
    ${ASMBL_SOURCE}
    ${IDBA_SOURCE}
    ${SEQ_SOURCE}
    ${CX1_SOURCE}
    ${TOOLKIT_SOURCE}
    ${MEGAHIT_EXTRAS}
)

# ── Executable: mcaat ─────────────────────────────────────────────────────────

add_executable(mcaat ${MAIN_SOURCE} ${MEGAHIT_SOURCES})

target_link_libraries(mcaat PRIVATE
    ${ZLIB_LIBRARIES}
    spoa::spoa
    fmt::fmt
    OpenMP::OpenMP_CXX
)

if(ENABLE_DEBUG)
    target_link_libraries(mcaat PRIVATE nlohmann_json::nlohmann_json)
    target_compile_definitions(mcaat PRIVATE DEBUG)
endif()

if(STATIC_BUILD)
    set_target_properties(mcaat PROPERTIES LINK_SEARCH_START_STATIC ON)
    target_link_options(mcaat PRIVATE -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -static)
endif()

