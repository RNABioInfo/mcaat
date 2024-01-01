#ifndef CYCLE_FINDER_H
#define CYCLE_FINDER_H

#include <../src/sdbg/sdbg.h>
#include <vector>
#include <map>
#include <set>
#include <array>
#include "helpers.h"

using namespace std;

class CycleFinder {

    private:
        vector<uint64_t> path;
        map<uint64_t, int> lock;
        vector<set<uint64_t>> stack;
        vector<int> backtrack_lengths;
        int maximal_length;
        int minimal_length;
        SDBG& succinct_de_bruijn_graph;
        string genome_name;
        
        /// @brief Find all nodes with indegree > 1, used for setting all the start nodes 
        /// @return uint64_t vector of all nodes with indegree > 1
        vector<uint64_t> _IndegreeMoreThanOne(){
            vector<uint64_t> nodes_more_than_one;
            for (uint64_t i = 0; i < this->succinct_de_bruijn_graph.size(); i++)
            {
                int indegree = this->succinct_de_bruijn_graph.EdgeIndegree(i);
                int multiplicity = this->succinct_de_bruijn_graph.EdgeMultiplicity(i);
                if (indegree >1)
                {
                    nodes_more_than_one.push_back(i);
                }
            }
            return nodes_more_than_one;
        };

        /// @brief Helper Function gets all outgoing edges of a node
        /// @param uint64_t node
        /// @return Outgoing edges of a node
        set<uint64_t> _GetOutgoings(uint64_t node){
            
            if (this->succinct_de_bruijn_graph.EdgeOutdegree(node) == 0)
                return {};
            
            uint64_t *outgoings = new uint64_t[this->succinct_de_bruijn_graph.EdgeOutdegree(node)];
            this->succinct_de_bruijn_graph.OutgoingEdges(node,outgoings);
            set<uint64_t> outgoings_set(outgoings, outgoings + this->succinct_de_bruijn_graph.EdgeOutdegree(node));
            return outgoings_set;
        };
        
        /// @brief Helper Function gets all incoming edges of a node
        /// @param uint64_t node
        /// @return Incoming edges of a node
        set<uint64_t> _GetIncomings(uint64_t node){
            
            if (this->succinct_de_bruijn_graph.EdgeIndegree(node) == 0)
                return {};
            
            uint64_t *incomings = new uint64_t[this->succinct_de_bruijn_graph.EdgeIndegree(node)];
            this->succinct_de_bruijn_graph.IncomingEdges(node,incomings);
            set<uint64_t> incomings_set(incomings, incomings + succinct_de_bruijn_graph.EdgeIndegree(node));
            return incomings_set;
        };

        /// @brief Set initial values for the cycle finding algorithm
        /// @param uint64_t start_node 
        void _InitialValues(uint64_t start_node){
            this->path.push_back(start_node);
            this->lock[start_node] = 0;
            this->stack.push_back(this->_GetOutgoings(start_node));
            this->backtrack_lengths.push_back(this->maximal_length);
        };
        
        /// @brief Mark a node as visited and add it to the path 
        /// @param uint64_t neighbor
        /// @return always false
        bool _MarkNode(uint64_t neighbor){
            this->path.push_back(neighbor);
            this->backtrack_lengths.push_back(maximal_length);
            this->lock[neighbor] = path.size();
            this->stack[this->stack.size() - 1].erase(neighbor);
            this->stack.push_back(_GetOutgoings(neighbor));
            return false;   
        };

        /// @brief Relax the lock of a node
        /// @param int backtrack_length
        void _RelaxLock(int backtrack_length, uint64_t v){
            vector<pair<int, int> > relax_stack;
            relax_stack.push_back(make_pair(backtrack_length, v));
            while (relax_stack.size() > 0) {
                int bl = relax_stack[relax_stack.size() - 1].first;
                int u = relax_stack[relax_stack.size() - 1].second;
                relax_stack.pop_back();
                if (this->lock.try_emplace(u,this->maximal_length).first->second < this->maximal_length - bl + 1) {
                    this->lock[u] = this->maximal_length - bl + 1;
                    for (auto w : this->_GetIncomings(u)) 
                        if (find(path.begin(), path.end(), w) == path.end()) 
                            relax_stack.push_back(make_pair(bl + 1, w));
                    }
                }
            return;
    }

    public:
        /// @brief Constructor for CycleFinder 
        /// @param SDBG& Succinct de Bruijn graph, int lengthBound
        CycleFinder(SDBG& succinct_de_bruijn_graph, int length_bound, int minimal_length, string genome_name) : 
            succinct_de_bruijn_graph(succinct_de_bruijn_graph), 
            maximal_length(length_bound),
            minimal_length(minimal_length),
            genome_name(genome_name) {
                // initialize values
                this->path = {};
                this->lock = {};
                this->stack = {};
                this->backtrack_lengths = {};
                this->FindCycles();
            };

        int FindCycle(uint64_t start_node){  
            // --- 3. Inner Loop: ---
            //--- 3.1. Count cycles and assign current_node to startnode ---
            int counter = 0;
            int current_node=start_node;
            // --- 3.2. Iterate through stack ---
            while(this->stack.size()>0){
                set<uint64_t> neighbors = this->stack.back(); // get neighbors of current node
                bool flag = true; 
                for (auto neighbor : neighbors){
                    current_node = neighbor;
                    if(current_node==start_node){
                        if (this->path.size() > this->minimal_length){
                            Helpers("a", this->succinct_de_bruijn_graph, 
                                    this->path, this->genome_name);
                            counter+=1;
                        }
                        this->backtrack_lengths[this->backtrack_lengths.size() - 1] = 1;
                    }
                    else if (this->path.size() < this->lock.try_emplace(neighbor, this->maximal_length).first->second)
                    {
                        neighbors.erase(neighbor);
                        flag=this->_MarkNode(neighbor);
                        break;
                    }
                }

                if(flag){
                    this->stack.pop_back();
                    uint64_t v = this->path[this->path.size() - 1];
                    this->path.pop_back();
                    int backtrack_length = this->backtrack_lengths[backtrack_lengths.size() - 1];
                    this->backtrack_lengths.pop_back();
                    if (this->backtrack_lengths.size() > 0) {
                        this->backtrack_lengths[this->backtrack_lengths.size() - 1] = 
                        min(this->backtrack_lengths[this->backtrack_lengths.size() - 1], backtrack_length);
                    }
                    if (backtrack_length < this->maximal_length) {
                        _RelaxLock(backtrack_length,v);
                    }
                }
            }
            return counter;
            
        };
        

        /// @brief Cycle Finding Algorithm
        /// based on a paper Finding All Bounded-Length Simple Cycles in a Directed Graph
        /// @return int number of cycles found in the graph        
        int FindCycles(){
            int cumulative = 0;
            // ----------------- 1. Find all nodes with indegree > 1 -----------------
            vector<uint64_t> start_nodes = this->_IndegreeMoreThanOne();
            
            // ------------------2. Outer Loop: ------------------
            while (start_nodes.size() > 0)
	        {
                // ---- 2.1. Pick a possible start_node and remove from stack ---
                uint64_t start_node = start_nodes.back();
		        start_nodes.pop_back();

                // ---- 2.2. Set initial values ----
                this->_InitialValues(start_node); 

                // ---- 2.3. Inner Loop: ----
                cumulative += this->FindCycle(start_node);

                // ---- 2.4. Set node as invalid(remove/subgraph) ----
                this->succinct_de_bruijn_graph.SetInvalidEdge(start_node);
            }

            std::cout << "Number of Cycles: " << cumulative << std::endl;
            return cumulative;
        };

        ~CycleFinder(){};
};  
#endif