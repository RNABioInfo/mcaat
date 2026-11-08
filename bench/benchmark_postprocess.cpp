#include "crispr_array/post_processing.h"
#include <chrono>
#include <iostream>
#include <omp.h>

int main(int argc,char** argv) {
    try {
        if (argc!=5) throw std::runtime_error("benchmark_postprocess GRAPH CYCLES OUTPUT SUPPORT");
        omp_set_num_threads(1);
        auto timer=std::chrono::steady_clock::now();
        SDBG graph;graph.LoadFromFile(argv[1]);
        std::cout<<"graph_load_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-timer).count()<<std::endl;
        Settings settings;settings.sdbg=&graph;settings.output_folder=argv[3];
        settings.dna_sequence_settings.repeat_support=mcaat::arrays::parse_repeat_support(argv[4]);
        timer=std::chrono::steady_clock::now();
        PostProcessor processor(settings);processor.run_analysis(argv[2]);
        std::cout<<"postprocess_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-timer).count()<<std::endl;
    } catch (const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
