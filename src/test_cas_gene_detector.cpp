#include <iostream>
#include "sdbg/sdbg.h"
#include "cas_workflow.h"

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: test_cas_gene_detector <graph> <profiles_dir> <repeat_kmer>" << std::endl;
        return 1;
    }

    std::string graph = argv[1];
    std::string profiles = argv[2];
    std::string repeat = argv[3];

    SDBG sdbg;
    sdbg.LoadFromFile(graph.c_str());

    // find repeat node
    uint64_t repeat_node = SDBG::kNullID;
    if (repeat.length() == (size_t)sdbg.k()) {
        uint8_t* enc = new uint8_t[sdbg.k()];
        for (uint32_t i = 0; i < sdbg.k(); ++i) {
            char c = repeat[i];
            enc[i] = (c == 'A' || c == 'a') ? 1 : (c=='C'||c=='c')?2:(c=='G'||c=='g')?3:4;
        }
        repeat_node = sdbg.IndexBinarySearch(enc);
        delete[] enc;
    }

    if (repeat_node == SDBG::kNullID) {
        std::cerr << "Repeat kmer not found" << std::endl;
        return 1;
    }

    CasWorkflow workflow(sdbg, profiles);
    workflow.SetParameters(50, 5000, 100, 41591, 10, 5000);

    auto res = workflow.DetectCasOperon(repeat_node);
    std::cout << "Detected " << res.genes.size() << " genes." << std::endl;

    for (size_t i = 0; i < res.genes.size(); ++i) {
        auto &g = res.genes[i];
        std::cout << "Gene #" << (i+1) << ": " << g.gene_name << " start=" << g.start_node << " dist=" << g.distance_from_repeat << " len=" << g.gene_length << " score=" << g.score << std::endl;
    }
    return 0;
}
