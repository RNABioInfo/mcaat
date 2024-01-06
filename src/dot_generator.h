#ifndef DOT_GENERATOR_H
#define DOT_GENERATOR_H

#include <../src/sdbg/sdbg.h>
#include <vector>
#include <map>
#include <set>
#include <array>

using namespace std;

class DotGenerator {

    private:
        SDBG& succinct_de_bruijn_graph;
        string genome_name;
        string output_file;
        void _Generate_Dot(){
            // generate graphviz dot file
            ofstream dot_file;
            dot_file.open(this->output_file);
            dot_file << "digraph G {\n";
            for (uint64_t i = 0; i < this->succinct_de_bruijn_graph.size(); i++)
            {
                int outdegree = this->succinct_de_bruijn_graph.EdgeOutdegree(i);
                if (outdegree > 0)
                {
                    uint64_t *outgoings = new uint64_t[outdegree];
                    this->succinct_de_bruijn_graph.OutgoingEdges(i,outgoings);
                    for (int j = 0; j < outdegree; j++)
                    {
                        dot_file << i << " -> " << outgoings[j] << ";\n";
                    }
                }

            }
            dot_file << "}";
            dot_file.close();
            //convert dot file to svg
            string command = "dot -Tsvg " + this->output_file + " -o " + this->output_file + ".svg";
            system(command.c_str());

        }
    public:
        DotGenerator(SDBG& succinct_de_bruijn_graph) :
            succinct_de_bruijn_graph(succinct_de_bruijn_graph),
            output_file("data/out/graph.dot")
        {
            this->_Generate_Dot();
        }


        ~DotGenerator() {
        }
};

#endif