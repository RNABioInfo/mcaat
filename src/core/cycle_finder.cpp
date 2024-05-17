#include "cycle_finder.h"


set<uint64_t> CycleFinder::_GetOutgoings(uint64_t node){
    return GetOutgoingNodes(this->sdbg, node);
}

set<uint64_t> CycleFinder::_GetIncomings(uint64_t node){
    return GetIncomingNodes(this->sdbg, node);
}

bool CycleFinder::_MarkNode(uint64_t neighbor, vector<uint64_t>& path, map<uint64_t, int>& lock, vector<set<uint64_t>>& stack, vector<int>& backtrack_lengths){
    
    path.push_back(neighbor);
    backtrack_lengths.push_back(this->maximal_length);
    lock[neighbor] = path.size();
    stack[stack.size() - 1].erase(neighbor);
    stack.push_back(_GetOutgoings(neighbor));
    return false;
}

void CycleFinder::_RelaxLock(int backtrack_length, uint64_t v, vector<uint64_t>& path, map<uint64_t, int>& lock, vector<set<uint64_t>>& stack, vector<int>& backtrack_lengths){
    
    vector<pair<int, int> > relax_stack;
    relax_stack.push_back(make_pair(backtrack_length, v));
    
    while (relax_stack.size() > 0) {
        int bl = relax_stack[relax_stack.size() - 1].first;
        int u = relax_stack[relax_stack.size() - 1].second;
        
        relax_stack.pop_back();
        if (lock.try_emplace(u,this->maximal_length).first->second < this->maximal_length - bl + 1) {
            lock[u] = this->maximal_length - bl + 1;
            for (auto w : _GetIncomings(u)) 
                if (find(path.begin(), path.end(), w) == path.end()) 
                    relax_stack.push_back(make_pair(bl + 1, w));
            }
        }
    return;
}

vector<node> CycleFinder::_SortNodesByMultiplicities(){
    vector<node> nodes;
    for (uint64_t i = 0; i < this->sdbg.size(); i++){
        if(this->sdbg.EdgeIndegree(i)>1){
            nodes.push_back(node(i, this->sdbg.EdgeMultiplicity(i)));
        }
    }
    //sort in descending order
    sort(nodes.begin(), nodes.end(), [](node a, node b) { return a.multiplicity > b.multiplicity; });
    return nodes;
}

vector<uint64_t> CycleFinder::_SortNodesByIds(){
    vector<uint64_t> nodes;
    for (uint64_t i = 0; i < this->sdbg.size(); i++){
        if(this->sdbg.EdgeIndegree(i)>1){
            nodes.push_back(i);
        }
    }
    

    sort(nodes.begin(), nodes.end());
    return nodes;
}

CycleFinder::CycleFinder(SDBG& sdbg, int length_bound, int minimal_length, string genome_name)
    :   maximal_length(length_bound),
        minimal_length(minimal_length),
        sdbg(sdbg),
        genome_name(genome_name),
        cluster_bounds(500)
{
    this->FindCycles();
}

vector<vector<uint64_t>> CycleFinder::FindCycle(uint64_t start_node, vector<uint64_t> path, map<uint64_t, int> lock, vector<set<uint64_t>> stack, vector<int> backtrack_lengths){
    int counter = 0;
    int current_node = start_node;
    bool count_flag = false;
    vector<vector<uint64_t>> cycles;
    // --- 3.2. Iterate through stack ---

    while(stack.size()>0){
        set<uint64_t> neighbors = stack.back(); // get neighbors of current node
        bool flag = true; 
        // create a dictionary with key as a edgeindegree and value as list of nodes
        //bool edge_indegree_flag = false;
        for (auto neighbor : neighbors){
            current_node = neighbor;
            if(current_node==start_node){
                backtrack_lengths[backtrack_lengths.size() - 1] = 1;
                if (path.size() > this->minimal_length){
                    counter+=1;
                    cycles.push_back(path);
                    if (counter >= this->cluster_bounds){ 
                        cycles.clear();
                        return cycles;          
                    }
                }
            //edge_indegree_flag = sdbg.EdgeIndegree(current_node) && sdbg.EdgeMultiplicity(current_node)>sdbg.EdgeMultiplicity(start_node) ? true : false;    
            }

            else if (path.size() < lock.try_emplace(neighbor, this->maximal_length).first->second)
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
            
            if (backtrack_length < this->maximal_length) 
                _RelaxLock(backtrack_length, v, path, lock, stack, backtrack_lengths);
        }
        
    }
    

    return cycles;
}

vector<vector<uint64_t>> CycleFinder::FindCycleUtil(uint64_t start_node){
    vector<uint64_t> path;
    map<uint64_t, int> lock;
    vector<set<uint64_t>> stack;
    vector<int> backtrack_lengths;
    path.push_back(start_node);
    lock[start_node] = 0;
    stack.push_back(this->_GetOutgoings(start_node));
    backtrack_lengths.push_back(maximal_length);
   vector<vector<uint64_t>> cycles = FindCycle(start_node, path, lock, stack, backtrack_lengths);
    return cycles;
}
//group the startnodes by multiplicity into map<int, vector<node>>
void CycleFinder::ChunkStartNodes(map<int, vector<uint64_t>, std::greater<int>>& start_nodes_chunked){
    vector<node> nodes = _SortNodesByMultiplicities();
    //initialize the map

    for (auto node : nodes){
        start_nodes_chunked[node.multiplicity].push_back(node.id);
    }
}
 



int CycleFinder::FindCycles(){
    int cumulative = 0;
    printf("Number of nodes in a graph: %lu\n", this->sdbg.size());
    int counter = 1;
    string mode = "fastq";
    this->visited.resize(this->sdbg.size(), false);
    std::cout << "Starting the cycle search:\n\n";
    map<int, vector<uint64_t>,std::greater<int>> start_nodes_chunked;
    ChunkStartNodes(start_nodes_chunked);
    progressbar bar(start_nodes_chunked.size());
    std::cout << "Number of chunks: " << start_nodes_chunked.size() << std::endl;
    bar.set_done_char("▓");
    for (const auto & [key, start_nodes] : start_nodes_chunked){
        bar.update();
        #pragma omp parallel for ordered num_threads(22) shared(visited) reduction(+:cumulative)
        for (uint64_t start_node_index = 0; start_node_index < start_nodes.size(); start_node_index++)
        {   
            uint64_t start_node = start_nodes[start_node_index];
            if(this->visited[start_nodes[start_node_index]]) continue;
            auto cycles = this->FindCycleUtil(start_node);
            cumulative+=cycles.size();
            for (auto cycle : cycles)
                for (auto this_node : cycle){
                    if (this_node==start_node)
                        continue;
                    this->visited[this_node] = true;
                }
            for (auto cycle : cycles)
                PathWriter("a", this->sdbg, cycle, this->genome_name, "fastq");
        }
    }
    
    std::cout<<endl;

    std::cout << "Number of Cycles: " << cumulative << std::endl;
    
    return cumulative;
}

CycleFinder::~CycleFinder(){}
