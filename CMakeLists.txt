cmake_minimum_required(VERSION 3.12)
project(megahit)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

find_package(TBB REQUIRED)


option(COVERAGE "Generate coverage report" OFF)
option(STATIC_BUILD "Build static executation" OFF)
option(SANITIZER "Enable sanitizers" OFF)
option(TSAN "Enable thread sanitizers" OFF)

include_directories(libs/megahit/src)
FILE(GLOB_RECURSE ASMBL_SOURCE "libs/megahit/src/assembly/*.cpp")
FILE(GLOB_RECURSE LCASM_SOURCE "libs/megahit/src/localasm/*.cpp")
FILE(GLOB_RECURSE IDBA_SOURCE "libs/megahit/src/idba/*.cpp")
FILE(GLOB_RECURSE SDBG_SOURCE "libs/megahit/src/sdbg/*.cpp")
FILE(GLOB_RECURSE CX1_SOURCE "libs/megahit/src/sorting/*.cpp")
FILE(GLOB_RECURSE SEQ_SOURCE "libs/megahit/src/sequence/*.cpp")
FILE(GLOB_RECURSE TOOLKIT_SOURCE "libs/megahit/src/tools/*.cpp")
FILE(GLOB_RECURSE PROGRESSBAR_SOURCE "libs/progressbar/include/*.hpp")
FILE(GLOB_RECURSE MAIN_SOURCE "src/*.cpp")
LIST(APPEND OTHER_SOURCE
        libs/megahit/src/main.cpp
        libs/megahit/src/main_assemble.cpp
        libs/megahit/src/main_buildlib.cpp
        libs/megahit/src/main_iterate.cpp
        libs/megahit/src/main_local_assemble.cpp
        libs/megahit/src/main_sdbg_build.cpp
        libs/megahit/src/utils/options_description.cpp
        )
LIST(APPEND SOME_OTHERS libs/megahit/src/utils/options_description.cpp)
if (STATIC_BUILD)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif (STATIC_BUILD)
add_subdirectory(libs/vendor/edlib EXCLUDE_FROM_ALL)

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)
find_package (SeqAn REQUIRED)
include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
# -fsanitize=address  -fsanitize=leak
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXXH_INLINE_ALL -ftemplate-depth=3000")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprefetch-loop-arrays -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__XROOT__='\"${CMAKE_SOURCE_DIR}/src\"'")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ZLIB_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# -fsanitize=address -fsanitize=leak
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O1 -D_LIBCPP_DEBUG -D_GLIBCXX_DEBUG")

add_executable(mcaat src/main.cpp ${ASMBL_SOURCE} ${IDBA_SOURCE} ${SDBG_SOURCE} ${LCASM_SOURCE}
        ${SEQ_SOURCE} ${CX1_SOURCE} ${TOOLKIT_SOURCE} ${SOME_OTHERS} ${TBB_LIBRARIES} ${PROGRESSBAR_SOURCE})

if (STATIC_BUILD)
    # TODO dirty
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -static")
    set_target_properties(mcaat PROPERTIES LINK_SEARCH_START_STATIC ON CXX_STANDARD 17)
endif (STATIC_BUILD)

target_link_libraries(mcaat ${ZLIB_LIBRARIES} ${TBB_LIBRARIES} ${SEQAN_LIBRARIES} edlib)

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
  target_compile_options(mcaat PUBLIC "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
  target_link_libraries(mcaat "${CMAKE_THREAD_LIBS_INIT}")
endif()