#ifndef CYCLE_FINDER_H
#define CYCLE_FINDER_H

#include <../src/sdbg/sdbg.h>
#include <vector>
#include <map>
#include <set>
#include <array>
#include <stack>
#include "path_writer.h"
#include <omp.h>
#include "helpers/node.h"
using namespace std;
class CycleFinder {

    private:
        uint8_t maximal_length;
        uint8_t minimal_length;
        SDBG& succinct_de_bruijn_graph;
        string genome_name;
        uint8_t cluster_bounds;
        bool* visited;
         

        /// @brief Helper Function gets all outgoing edges of a node
        /// @param uint64_t node
        /// @return Outgoing edges of a node
        set<uint64_t> _GetOutgoings(uint64_t node){
            return GetOutgoingNodes(succinct_de_bruijn_graph, node);
        };
        
        /// @brief Helper Function gets all incoming edges of a node
        /// @param uint64_t node
        /// @return Incoming edges of a node
        set<uint64_t> _GetIncomings(uint64_t node){
            return GetIncomingNodes(succinct_de_bruijn_graph, node);
        };
        void WriteInfo(string cycles){
            ofstream info_file;
            string info_filename = "proof_of_concept/data/"+genome_name+"/cycles/info.txt";
            info_file.open(info_filename,std::ios_base::app);
            info_file << "Number of cycles: " << cycles << "\n";
            info_file.close();
        }
        /// @brief Set initial values for the cycle finding algorithm
        /// @param uint64_t start_node 
        
        
        /// @brief Mark a node as visited and add it to the path 
        /// @param uint64_t neighbor
        /// @return always false
        bool _MarkNode(uint64_t neighbor, vector<uint64_t>& path, map<uint64_t, int>& lock, vector<set<uint64_t>>& stack, vector<int>& backtrack_lengths){
            path.push_back(neighbor);
            backtrack_lengths.push_back(maximal_length);
            lock[neighbor] = path.size();
            stack[stack.size() - 1].erase(neighbor);
            stack.push_back(_GetOutgoings(neighbor));
            return false;
        };

        /// @brief Relax the lock of a node
        /// @param int backtrack_length
        
        void _RelaxLock(int backtrack_length, uint64_t v, vector<uint64_t>& path, map<uint64_t, int>& lock, vector<set<uint64_t>>& stack, vector<int>& backtrack_lengths){
            
            vector<pair<int, int> > relax_stack;
            relax_stack.push_back(make_pair(backtrack_length, v));
            
            while (relax_stack.size() > 0) {
                int bl = relax_stack[relax_stack.size() - 1].first;
                int u = relax_stack[relax_stack.size() - 1].second;
                
                relax_stack.pop_back();
                if (lock.try_emplace(u,maximal_length).first->second < maximal_length - bl + 1) {
                    lock[u] = maximal_length - bl + 1;
                    for (auto w : _GetIncomings(u)) 
                        if (find(path.begin(), path.end(), w) == path.end()) 
                            relax_stack.push_back(make_pair(bl + 1, w));
                    }
                }
            return;
    }
    

    vector<node> _SortNodesByMultiplicities(){
        //nodes is array of nodes with multiple edges
        vector<node> nodes;
        for (uint64_t i = 0; i < succinct_de_bruijn_graph.size(); i++){
            if(succinct_de_bruijn_graph.EdgeIndegree(i)>1){
                nodes.push_back(node(i, succinct_de_bruijn_graph.EdgeMultiplicity(i)));
            }
        }
        

        sort(nodes.begin(), nodes.end());
        return nodes;
    }
    vector<uint64_t> _SortNodesByIds(){
        //nodes is array of nodes with multiple edges
        vector<uint64_t> nodes;
        for (uint64_t i = 0; i < succinct_de_bruijn_graph.size(); i++){
            if(succinct_de_bruijn_graph.EdgeIndegree(i)>1){
                nodes.push_back(i);
            }
        }
        

        sort(nodes.begin(), nodes.end());
        return nodes;
    }
    public:
        /// @brief Constructor for CycleFinder 
        /// @param SDBG& Succinct de Bruijn graph, int lengthBound
        CycleFinder(SDBG& succinct_de_bruijn_graph, int length_bound, int minimal_length, string genome_name)
            :   maximal_length(length_bound),
                minimal_length(minimal_length),
                succinct_de_bruijn_graph(succinct_de_bruijn_graph),
                genome_name(genome_name),
                cluster_bounds(300),
                visited(new bool[succinct_de_bruijn_graph.size()])
            {
                FindCycles();
            };

        int FindCycle(uint64_t start_node, vector<uint64_t> path, map<uint64_t, int> lock, vector<set<uint64_t>> stack, vector<int> backtrack_lengths){  
            // --- 3. Inner Loop: ---
            //--- 3.1. Count cycles and assign current_node to startnode ---
            int counter = 0;
            int current_node = start_node;
            bool count_flag = false;
            vector<vector<uint64_t>> cycles;
            // --- 3.2. Iterate through stack ---

            while(stack.size()>0){
                map<int, set<uint64_t>> keep_nodes;
                set<uint64_t> neighbors = stack.back(); // get neighbors of current node
                bool flag = true; 
                // create a dictionary with key as a edgeindegree and value as list of nodes
               
                for (auto neighbor : neighbors){
                    current_node = neighbor;
                    if(current_node==start_node){
                        if (path.size() > minimal_length){
                            
                            counter+=1;
                            cycles.push_back(path);
                            
                            if (counter >= cluster_bounds) 
                                return 0;          
                        }
                        backtrack_lengths[backtrack_lengths.size() - 1] = 1;
                    }
                    else if (path.size() < lock.try_emplace(neighbor, maximal_length).first->second)
                    {
                        neighbors.erase(neighbor);
                        flag=_MarkNode(neighbor, path, lock, stack, backtrack_lengths);
                        break;
                    }
                }         
                
                
                if(flag){
                    stack.pop_back();
                    uint64_t v = path[path.size() - 1];
                    path.pop_back();
                    int backtrack_length = backtrack_lengths[backtrack_lengths.size() - 1];
                    backtrack_lengths.pop_back();
                    
                    if (backtrack_lengths.size() > 0) 
                        backtrack_lengths[backtrack_lengths.size() - 1] = 
                        min(backtrack_lengths[backtrack_lengths.size() - 1], backtrack_length);
                    
                    if (backtrack_length < maximal_length) 
                        _RelaxLock(backtrack_length, v, path, lock, stack, backtrack_lengths);
                }
                
            }
            
            for (auto cycle : cycles)
                for (auto node : cycle){
                    if (node==start_node)
                        continue;
                    if(succinct_de_bruijn_graph.EdgeIndegree(node)>1)
                        visited[node] = true;
                    }
            for (auto cycle : cycles)
                PathWriter("a", succinct_de_bruijn_graph, 
                             cycle, genome_name, "fasta");
            return counter;
            
        };
        

        /// @brief Cycle Finding Algorithm
        /// based on a paper Finding All Bounded-Length Simple Cycles in a Directed Graph
        /// @return int number of cycles found in the graph
        int FindCycles(){
            int cumulative = 0;
            printf("Number of nodes in a graph: %lu\n", succinct_de_bruijn_graph.size());
            int counter = 0;
            vector<uint64_t> path;
            map<uint64_t, int> lock;
            vector<set<uint64_t>> stack;
            vector<int> backtrack_lengths;
            fill_n(visited, succinct_de_bruijn_graph.size(), false);
            //vector<node> start_nodes = _SortNodesByMultiplicities();
            
            
            string mode = "fastq";
            if (mode=="fasta"){
                vector<uint64_t> start_nodes = _SortNodesByIds();
                cout << "Number of nodes with multiple edges: " << start_nodes.size() << endl;
                for (uint64_t start_node_index = start_nodes.size()-1; start_node_index > 0; start_node_index--)
                {   
                    //if(visited[start_nodes[start_node_index]]) continue;
                    uint64_t start_node = start_nodes[start_node_index];
                    counter += 1;
                    path.push_back(start_node);
                    lock[start_node] = 0;
                    stack.push_back(_GetOutgoings(start_node));
                    backtrack_lengths.push_back(maximal_length);
                    cumulative += FindCycle(start_node, path, lock, stack, backtrack_lengths);
                    
                    path.clear();
                    lock.clear();
                    stack.clear();
                    backtrack_lengths.clear();
                    succinct_de_bruijn_graph.SetInvalidEdge(start_node);
                }
            }
            else{
                vector<node> start_nodes = _SortNodesByMultiplicities();
                cout << "Number of nodes with multiple edges: " << start_nodes.size() << endl;
                for (uint64_t start_node_index = start_nodes.size()-1; start_node_index > 0; start_node_index--)
                {   
                    if(visited[start_nodes[start_node_index].id]) continue;
                    uint64_t start_node = start_nodes[start_node_index].id;
                    counter += 1;
                    path.push_back(start_node);
                    lock[start_node] = 0;
                    stack.push_back(_GetOutgoings(start_node));
                    backtrack_lengths.push_back(maximal_length);
                    cumulative += FindCycle(start_node, path, lock, stack, backtrack_lengths);
                    
                    path.clear();
                    lock.clear();
                    stack.clear();
                    backtrack_lengths.clear();
                    succinct_de_bruijn_graph.SetInvalidEdge(start_node);
                }
            }
            //print out number of visited nodes
            int visited_nodes = 0;
            for (uint64_t i = 0; i < succinct_de_bruijn_graph.size(); i++)
                if(visited[i]) visited_nodes++;
            cout << "Number of visited nodes: " << visited_nodes << endl;
            std::cout << "Number of Cycles: " << cumulative << std::endl;
            return cumulative;
        };

        ~CycleFinder(){};
};  
#endif