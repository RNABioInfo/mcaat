#ifndef GRAPH_WRITER_H
#define GRAPH_WRITER_H

#include <../src/sdbg/sdbg.h>
#include <vector>
#include <map>
#include <set>
#include <array>
#include "helpers/generic_functions.h"
using namespace std;

class GraphWriter {

    private:
        SDBG& succinct_de_bruijn_graph;
        string output_file_name;

    public:
        GraphWriter(SDBG& sdbg, string output_file_name = "data/graph/graph") : 
        succinct_de_bruijn_graph(sdbg), output_file_name(output_file_name) {};

        /// @brief Generate dot file from SDBG
        /// @return void
        void GenerateDotFile(){
            ofstream dot_file;
            dot_file.open(this->output_file_name+".dot");
            dot_file << "digraph G {" << endl;
            for (uint64_t i = 0; i < this->succinct_de_bruijn_graph.size(); i++)
            {
                int outdegree = this->succinct_de_bruijn_graph.EdgeOutdegree(i);
                if (outdegree > 0)
                {
                    uint64_t *outgoings = new uint64_t[outdegree];
                    this->succinct_de_bruijn_graph.OutgoingEdges(i,outgoings);
                    for (int j = 0; j < outdegree; j++)
                    {
                        dot_file << i << " -> " << outgoings[j] << endl;
                    }
                }
            }
            dot_file << "}" << endl;
            dot_file.close();
            //convert to svg
            string command = "dot -Tsvg " + this->output_file_name + " -o " + this->output_file_name + ".svg";
            system(command.c_str());
        };

        /// @brief Write nodes with SDBG labels to file
        /// @return void

        void WriteNodes(){
            ofstream node_file;
            node_file.open(this->output_file_name+".nodes");
            for (uint64_t node = 0; node < this->succinct_de_bruijn_graph.size(); node++)
            {
                node_file << node << " : " << FetchNodeLabel(this->succinct_de_bruijn_graph,node) <<endl;
            }
            node_file.close();
        };
        ~GraphWriter(){};
};
#endif