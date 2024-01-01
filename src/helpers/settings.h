#include <string>

/* the following are settings
input folder: in/001/
output graph: graphs/001/
output cycles: 18102023/cycles/
host_memory: 128000000000
mem: 1
k: 23
num_cpu_threads: 20
*/

struct settings{
    std::string input_folder;
    std::string output_graph;
    std::string output_cycle;
    std::string host_memory;
    std::string mem;
    std::string threads;
    std::string kmer_size;
};