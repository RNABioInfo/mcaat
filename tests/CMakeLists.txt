cmake_minimum_required(VERSION 3.16)
project(mcaat_search_checks LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
enable_testing()
get_filename_component(MCAAT_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(MEGAHIT_SRC "${MCAAT_ROOT}/libs/megahit/src")
add_executable(test_unitig_cycle_search test_unitig_cycle_search.cpp)
target_include_directories(test_unitig_cycle_search PRIVATE "${MCAAT_ROOT}/include")
target_include_directories(test_unitig_cycle_search SYSTEM PRIVATE "${MEGAHIT_SRC}")
add_test(NAME unitig_cycle_search COMMAND test_unitig_cycle_search)
option(BUILD_SDBG_BENCHMARK "Build the actual-SDBG benchmark" ON)
option(SEARCH_SANITIZERS "Check search tests with ASan/UBSan" OFF)
if(SEARCH_SANITIZERS)
    target_compile_options(test_unitig_cycle_search PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(test_unitig_cycle_search PRIVATE -fsanitize=address,undefined)
endif()
if(BUILD_SDBG_BENCHMARK)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    find_package(OpenMP REQUIRED)
    add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/original_precheck.h"
        COMMAND "${Python3_EXECUTABLE}" "${MCAAT_ROOT}/bench/extract_precheck.py"
            "${MCAAT_ROOT}/src/crispr_array/cycle_finder.cpp"
            "${CMAKE_CURRENT_BINARY_DIR}/original_precheck.h"
        DEPENDS "${MCAAT_ROOT}/bench/extract_precheck.py"
                "${MCAAT_ROOT}/src/crispr_array/cycle_finder.cpp")
    add_executable(benchmark_unitig_search
        "${MCAAT_ROOT}/bench/benchmark_unitig_search.cpp"
        "${CMAKE_CURRENT_BINARY_DIR}/original_precheck.h"
        "${MEGAHIT_SRC}/sdbg/sdbg_raw_content.cpp"
        "${MEGAHIT_SRC}/sdbg/sdbg_meta.cpp")
    target_include_directories(benchmark_unitig_search PRIVATE
        "${MCAAT_ROOT}/include" "${CMAKE_CURRENT_BINARY_DIR}")
    target_include_directories(benchmark_unitig_search SYSTEM PRIVATE "${MEGAHIT_SRC}")
    target_link_libraries(benchmark_unitig_search PRIVATE OpenMP::OpenMP_CXX)
endif()

add_executable(test_graph_repeat_filter test_graph_repeat_filter.cpp)
target_include_directories(test_graph_repeat_filter PRIVATE "${MCAAT_ROOT}/include")
add_test(NAME graph_repeat_filter COMMAND test_graph_repeat_filter)
add_executable(test_candidate_quality test_candidate_quality.cpp)
target_include_directories(test_candidate_quality PRIVATE "${MCAAT_ROOT}/include")
add_test(NAME candidate_quality COMMAND test_candidate_quality)
add_executable(benchmark_repeat_filter "${MCAAT_ROOT}/bench/benchmark_repeat_filter.cpp")
target_include_directories(benchmark_repeat_filter PRIVATE "${MCAAT_ROOT}/include")
if(SEARCH_SANITIZERS)
    target_compile_options(test_candidate_quality PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(test_candidate_quality PRIVATE -fsanitize=address,undefined)
    target_compile_options(test_graph_repeat_filter PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(test_graph_repeat_filter PRIVATE -fsanitize=address,undefined)
endif()

if(BUILD_SDBG_BENCHMARK)
    add_executable(benchmark_postprocess "${MCAAT_ROOT}/bench/benchmark_postprocess.cpp"
        "${MEGAHIT_SRC}/sdbg/sdbg_raw_content.cpp"
        "${MEGAHIT_SRC}/sdbg/sdbg_meta.cpp")
    target_include_directories(benchmark_postprocess PRIVATE "${MCAAT_ROOT}/include")
    target_include_directories(benchmark_postprocess SYSTEM PRIVATE "${MEGAHIT_SRC}")
    target_link_libraries(benchmark_postprocess PRIVATE OpenMP::OpenMP_CXX)
    add_executable(test_graph_path_io test_graph_path_io.cpp
        "${MEGAHIT_SRC}/sdbg/sdbg_raw_content.cpp"
        "${MEGAHIT_SRC}/sdbg/sdbg_meta.cpp")
    target_include_directories(test_graph_path_io PRIVATE "${MCAAT_ROOT}/include")
    target_include_directories(test_graph_path_io SYSTEM PRIVATE "${MEGAHIT_SRC}")
    target_link_libraries(test_graph_path_io PRIVATE OpenMP::OpenMP_CXX)
endif()

add_executable(replay_array_candidates "${MCAAT_ROOT}/bench/replay_array_candidates.cpp")
target_include_directories(replay_array_candidates PRIVATE "${MCAAT_ROOT}/include")
