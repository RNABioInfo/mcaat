#include "cycle_finder.h"


void CycleFinder::_GetOutgoings(uint64_t node, set<uint64_t>& outgoings_set){
    int edge_outdegree = sdbg.EdgeOutdegree(node);
    if (edge_outdegree == 0 || node == -1){
        outgoings_set = {};
        return ;
    }
    uint64_t *outgoings= new uint64_t[edge_outdegree];
    sdbg.OutgoingEdges(node,outgoings);
    for (int i = 0; i < edge_outdegree; i++){ 
        outgoings_set.insert(outgoings[i]);
    }
    delete[] outgoings;
}
void CycleFinder::_GetIncomings(uint64_t node, set<uint64_t>& incomings_set){
    int edge_indegree = sdbg.EdgeIndegree(node);
    if (edge_indegree == 0 || node == -1){
        incomings_set = {};
        return;
    }
    uint64_t *incomings= new uint64_t[edge_indegree];
    sdbg.IncomingEdges(node,incomings);
    for (int i = 0; i < edge_indegree; i++){ 
        incomings_set.insert(incomings[i]);
    }
    delete[] incomings;
}

bool CycleFinder::_MarkNode(uint64_t neighbor, vector<uint64_t>& path, map<uint64_t, int>& lock, vector<set<uint64_t>>& stack, vector<int>& backtrack_lengths){
    
    
    return false;
}

void CycleFinder::_RelaxLock(int backtrack_length, uint64_t v, vector<uint64_t>& path, map<uint64_t, int>& lock, vector<set<uint64_t>>& stack, vector<int>& backtrack_lengths){
    
    
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
    //sort(nodes.begin(), nodes.end(), [](node a, node b) { return a.multiplicity > b.multiplicity; });
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
    int steps_counter = 0;
    while(stack.size()>0){
        set<uint64_t> neighbors = stack.back(); // get neighbors of current node
        steps_counter+=1;
        
        
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
                path.push_back(neighbor);
                backtrack_lengths.push_back(this->maximal_length);
                lock[neighbor] = path.size();
                stack[stack.size() - 1].erase(neighbor);
                set<uint64_t> outgoings;
                this->_GetOutgoings(neighbor, outgoings);
                stack.push_back(outgoings);
                flag=false;
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
            
            if (backtrack_length < this->maximal_length) {
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
            }
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
    set<uint64_t> outgoings;
    this->_GetOutgoings(start_node, outgoings);
    stack.push_back(outgoings);
    backtrack_lengths.push_back(maximal_length);
    
    vector<vector<uint64_t>> cycles = FindCycle(start_node, path, lock, stack, backtrack_lengths);
    lock.clear();
    stack.clear();
    backtrack_lengths.clear();
    path.clear();
    return cycles;
}
//group the startnodes by multiplicity into map<int, vector<node>>
void CycleFinder::ChunkStartNodes(map<int, vector<uint64_t>, std::greater<int>>& start_nodes_chunked){
    //vector<node> nodes = _SortNodesByMultiplicities();
    //initialize the map
    //find log2 of the number of nodes
    std::mutex mtx;

    std::vector<std::thread> threads;
    for (int i = 0; i < 22; ++i) {
        threads.emplace_back([&, i] {
            for (uint64_t node = i; node < this->sdbg.size(); node += 22) {
                if (this->sdbg.EdgeIndegree(node) >= 2 && this->sdbg.EdgeMultiplicity(node) > 20) {
                    double log2_mult = ceil(log2(double(this->sdbg.EdgeMultiplicity(node))));
                    mtx.lock();
                    start_nodes_chunked[log2_mult].push_back(node);
                    mtx.unlock();
                }
            }
        });
    }

    for (auto& thread : threads) {
        thread.join();
    }

    // print the start_nodes_chunked size
    for (auto [key, value] : start_nodes_chunked)
        std::cout << "log2_mult: " << key << " Number of nodes: " << value.size() << std::endl;
}


int CycleFinder::FindCycles(){
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
    /*
    
    for (auto nodes_iterator = start_nodes_chunked.begin(); nodes_iterator != start_nodes_chunked.end(); ){
        //bar.update();
        auto thread_count = 22;
        if (nodes_iterator->second.size() < thread_count)
            thread_count = nodes_iterator->second.size();
        #pragma omp parallel for num_threads(4) shared(sdbg) reduction(+:cumulative)
        for (uint64_t start_node_index = 0; start_node_index < nodes_iterator->second.size(); start_node_index++)
        {   
            
            uint64_t start_node = nodes_iterator->second[start_node_index];
            if(this->visited[nodes_iterator->second[start_node_index]]) continue;
            vector<vector<uint64_t>> cycles = this->FindCycleUtil(start_node);
            cumulative+=cycles.size();
            
            for (auto cycle : cycles){
                    for (auto this_node : cycle){
                        visited[this_node] = true;
                    }
                //PathWriter("a", this->sdbg, cycle, this->genome_name, "fastq");
            }
            
            counter+=1;
            cycles.clear();
        }
        cout << "processed " << nodes_iterator->second.size() << " ";
        counter=0;
        cout<<"Press something to continue\n"<<endl;
        //    std::cin.ignore();
        nodes_iterator = start_nodes_chunked.erase(nodes_iterator);
        printf("number of Cycles: %d\n", cumulative);
    }
    */std::mutex mtx;
    auto nodes_iterator = start_nodes_chunked.begin();
    while (nodes_iterator != start_nodes_chunked.end()) {
        std::vector<std::thread> threads;
        uint64_t cumulative = 0;
        int counter = 0;

        auto thread_count = 22;
        if (nodes_iterator->second.size() < thread_count)
            thread_count = nodes_iterator->second.size();

        for (int i = 0; i < thread_count; ++i) {
            threads.emplace_back([&, i] {  // Capture variables by reference and 'i' by value
                for (uint64_t start_node_index = i; start_node_index < nodes_iterator->second.size(); start_node_index += thread_count) {
                    uint64_t start_node = nodes_iterator->second[start_node_index];
                    mtx.lock();
                    if(this->visited[start_node]) {
                        mtx.unlock();
                        continue;
                    }
                    mtx.unlock();
                    vector<vector<uint64_t>> cycles = this->FindCycleUtil(start_node);
                    mtx.lock();
                    cumulative += cycles.size();
                    mtx.unlock();

                    for (auto cycle : cycles){
                        for (auto this_node : cycle){
                            visited[this_node] = true;
                        }
                        //PathWriter("a", this->sdbg, cycle, this->genome_name, "fastq");
                    }
                    cycles.clear();
                }
            });
        }

        for (auto& thread : threads) {
            thread.join();
        }

        cout << "processed " << nodes_iterator->second.size() << " ";
        counter = 0;
        cout << "Press something to continue\n" << endl;
        // std::cin.ignore();
        nodes_iterator = start_nodes_chunked.erase(nodes_iterator);
        printf("number of Cycles: %d\n", cumulative);
    }
    std::cout<<endl;

    std::cout << "Number of Cycles: " << cumulative << std::endl;
    
    return cumulative;
}

CycleFinder::~CycleFinder(){}
 