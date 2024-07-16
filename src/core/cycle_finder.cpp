#include "cycle_finder.h"


void CycleFinder::_GetOutgoings(uint64_t node, std::set<uint64_t>& outgoings_set) {
    // if outdegree is 0, assign an empty set to outgoings_set and return
    int edge_outdegree = sdbg.EdgeOutdegree(node);
    if (edge_outdegree == 0 || node == -1) {
        
        outgoings_set = {};
        return;
    }
    // get the outgoings of the node if the node_multiplicity is greater than 20
    uint64_t outgoings[edge_outdegree];
    sdbg.OutgoingEdges(node, outgoings);
    for (const auto& outgoing : outgoings)
        if (this->sdbg.EdgeMultiplicity(outgoing) > 20)
            outgoings_set.insert(outgoing);
    return;
}

void CycleFinder::_GetIncomings(uint64_t node, std::set<uint64_t>& incomings_set) {
    int edge_indegree = sdbg.EdgeIndegree(node);
    
    if (edge_indegree == 0 || node == -1) {
        incomings_set = {};
        return;
    }
    uint64_t incomings[edge_indegree];
    sdbg.IncomingEdges(node, incomings);
    for (const auto& incoming : incomings) 
        if (this->sdbg.EdgeMultiplicity(incoming) > 20)
            incomings_set.insert(incoming);
    
    return;
}



bool CycleFinder::_MarkNode(uint64_t neighbor, vector<uint64_t>& path, map<uint64_t, int>& lock, vector<set<uint64_t>>& stack, vector<int>& backtrack_lengths){
    path.push_back(neighbor);
    backtrack_lengths.push_back(this->maximal_length);
    lock[neighbor] = path.size();
    stack[stack.size() - 1].erase(neighbor);
    set<uint64_t> outgoings;
    this->_GetOutgoings(neighbor, outgoings);
    stack.push_back(outgoings);
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
            set<uint64_t> incomings;
            this->_GetIncomings(u, incomings);
            for (auto w : incomings)
                if (find(path.begin(), path.end(), w) == path.end())
                    relax_stack.push_back(make_pair(bl + 1, w));
        }
    }
    return;
}


CycleFinder::CycleFinder(SDBG& sdbg, int length_bound, int minimal_length, string genome_name)
        :   maximal_length(length_bound),
            minimal_length(minimal_length),
            sdbg(sdbg),
            genome_name(genome_name),
            cluster_bounds(500)
{
    /*
    vector<uint64_t> tips;
    #pragma omp parallel for num_threads(22)
    for (uint64_t node = 0; node < this->sdbg.size(); node++)
        if (this->sdbg.EdgeOutdegree(node) == 0 && this->sdbg.EdgeIndegree(node) == 1){
            #pragma omp critical
            tips.push_back(node);
        }
    auto count_of_nodes_less_twenty=0;
    auto rest = 0;
    //count the number of nodes with multiplicities greater than 20
    for (uint64_t node = 0; node < this->sdbg.size(); node++)
        if (this->sdbg.EdgeMultiplicity(node) > 20)
            count_of_nodes_less_twenty+=1;
        else
            rest+=1;
    cout << "Number of nodes with multiplicities greater than 20: " << count_of_nodes_less_twenty << endl;
    cout << "Number of nodes with multiplicities less than 20: " << rest << endl;
    cout << "Graph size: " << this->sdbg.size() << endl;
    cout << "Number of tips: " << tips.size() << endl;
    char num_threads = tips.size() < 22 ? tips.size() : 22;
    size_t removed_nodes = 0;
    #pragma omp parallel for num_threads(22)
    for (auto tip : tips)
        removed_nodes += this->RemoveLinearPathsBacktrack(tip);
    cout << "Number of removed nodes: " << removed_nodes << endl;
    */
    this->FindCycles();
}

size_t CycleFinder::RemoveLinearPathsBacktrack(uint64_t node){
    vector<uint64_t> stack;
    stack.push_back(node);
    size_t number_of_nodes_removed = 0;
    while (stack.size() > 0){
        uint64_t current_node = stack[stack.size() - 1];
        stack.pop_back();
        char indegree = this->sdbg.EdgeIndegree(current_node);
        char outdegree = this->sdbg.EdgeOutdegree(current_node);
        //bool first_node = (current_node == node);
        if(indegree <= 1 && outdegree<2){
            uint64_t incomings[indegree];
            this->sdbg.IncomingEdges(current_node, incomings);
            number_of_nodes_removed+=1;
            if (incomings[0] != -1){
                stack.push_back(incomings[0]);
                #pragma omp critical
                this->sdbg.SetInvalidEdge(current_node);
            }
        }
    } 
    return number_of_nodes_removed;
}
vector<vector<uint64_t>> CycleFinder::FindCycle(uint64_t start_node, vector<uint64_t> path, map<uint64_t, int> lock, vector<set<uint64_t>> stack, vector<int> backtrack_lengths){
    int counter = 0;
    int current_node = start_node;
    bool count_flag = false;
    vector<vector<uint64_t>> cycles;
    // --- 3.2. Iterate through stack ---
    int steps_counter = 0;
    size_t current_mem_info;
    while(stack.size()>0){
        
        set<uint64_t> neighbors = stack.back(); // get neighbors of current node
        steps_counter+=1;
        bool flag = true;
        if (steps_counter == 10000000){
            struct mallinfo mem_info = mallinfo();
            current_mem_info = mem_info.uordblks;
            cout << "Memory allocated: " << current_mem_info << endl;
            malloc_trim(0);
            return cycles;
        }
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
            {   neighbors.erase(neighbor);
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
                this->_RelaxLock(backtrack_length, v, path, lock, stack, backtrack_lengths);     
            
        }
    }
    malloc_trim(0);
    return cycles;
}

vector<vector<uint64_t>> CycleFinder::FindCycleUtil(uint64_t start_node){
    vector<uint64_t> path;
    map<uint64_t, int> lock;
    vector<set<uint64_t>> stack;
    vector<int> backtrack_lengths;
    path.push_back(start_node);
    lock[start_node] = 0;
    set<uint64_t> outgoings;
    this->_GetOutgoings(start_node, outgoings);
    stack.push_back(outgoings);
    backtrack_lengths.push_back(maximal_length);
    vector<vector<uint64_t>> cycles = FindCycle(start_node, path, lock, stack, backtrack_lengths);
    return cycles;
}
//group the startnodes by multiplicity into map<int, vector<node>>
void CycleFinder::ChunkStartNodes(map<int, vector<uint64_t>, std::greater<int>>& start_nodes_chunked){
    //vector<node> nodes = _SortNodesByMultiplicities();
    //initialize the map
    //find log2 of the number of nodes
    #pragma omp parallel for num_threads(22)
    for (uint64_t node = 0; node < this->sdbg.size(); node++)
    {
        if (this->sdbg.EdgeIndegree(node) >=2 && this->sdbg.EdgeMultiplicity(node) > 20 && this->sdbg.IsValidEdge(node)){
            double log2_mult = ceil(log2(double(this->sdbg.EdgeMultiplicity(node))));
            #pragma omp critical
            {
                start_nodes_chunked[log2_mult].push_back(node);
            }
            size_t pad = 0; // Example: Release all unused memory
            malloc_trim(pad);
        }
    }

    //print the start_nodes_chunked size
    for (auto [key, value] : start_nodes_chunked)
        std::cout << "log2_mult: " << key << " Number of nodes: " << value.size() << std::endl;
    
}


int CycleFinder::FindCycles(){
    struct mallinfo mem_info = mallinfo();
    size_t graph_mem_info = mem_info.uordblks;
    int cumulative = 0;
    printf("Number of nodes in a graph: %lu\n", this->sdbg.size());
    string mode = "fastq";
    this->visited.resize(this->sdbg.size(), false);
    std::cout << "Starting the cycle search:\n\n";
    map<int, vector<uint64_t>,std::greater<int>> start_nodes_chunked;
    ChunkStartNodes(start_nodes_chunked);
    progressbar bar(start_nodes_chunked.size());
    std::cout << "Number of chunks: " << start_nodes_chunked.size() << std::endl;
    bar.set_done_char("▓");
    int counter = 0;
    std::cin.ignore();
    for (auto nodes_iterator = start_nodes_chunked.begin(); nodes_iterator != start_nodes_chunked.end(); ){
        //bar.update();
        
        auto thread_count = 22;
        if (nodes_iterator->second.size() < thread_count)
            thread_count = nodes_iterator->second.size();
        #pragma omp parallel for num_threads(thread_count) reduction(+:cumulative) shared(nodes_iterator, sdbg, visited)
        for (uint64_t start_node_index = 0; start_node_index < nodes_iterator->second.size(); start_node_index++)
        {
            malloc_trim(0);
            // pick a start_node
            uint64_t start_node = nodes_iterator->second[start_node_index];
            // check if the node has been visited
            if(this->visited[nodes_iterator->second[start_node_index]]) continue;
            // invoke the FindCycleUtil method
            vector<vector<uint64_t>> cycles = this->FindCycleUtil(start_node);
            // update the number of cycles
            cumulative+=cycles.size();
            // write the cycles to a file
            for (auto cycle : cycles){
                for (auto this_node : cycle){
                    visited[this_node] = true;
                }
                PathWriter("a", this->sdbg, cycle, this->genome_name, "fastq");
            }
            malloc_trim(0);
            // Release all unused heap memory
        }
        
        cout << "processed " << nodes_iterator->second.size() << " ";
        cout<<"Press something to continue\n"<<endl;
        //std::cin.ignore();
        nodes_iterator = start_nodes_chunked.erase(nodes_iterator);
        printf("number of Cycles: %d\n", cumulative);
       
    }
    std::cout<<endl;

    std::cout << "Number of Cycles: " << cumulative << std::endl;

    return cumulative;
}

CycleFinder::~CycleFinder(){}
