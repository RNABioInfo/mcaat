#include <sdbg/sdbg.h>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <string>
#include <stack>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <cstdio>
#include <cstring>
using std::istringstream;
using namespace std;

class Filters {
private:
    SDBG& sdbg;
    unordered_map<uint64_t, vector<vector<uint64_t>>> cycles;

    string _FetchNodeLabel(size_t node);
    char _FetchLastCharacter(size_t node);

public:

    Filters(SDBG& sdbg, std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles);
    bool DFS(uint64_t start, uint64_t goal, const std::unordered_set<uint64_t>& valid_nodes, 
             std::vector<uint64_t>& path);
    std::vector<uint64_t> FindRepeatNodePaths(vector<uint64_t> repeat_nodes,uint64_t start_node);
    void rotateLeft(std::vector<uint64_t>& arr, int k);
    pair<vector<uint64_t>, vector<vector<uint64_t>>> _FindCRISPRArrayNodes(uint64_t start_node);
   unordered_map<string, string>  ListArrays(int& number_of_spacers);
   string findMostFrequentSequence(const std::string& input, int minLength, int maxLength, int minDistance, int maxDistance);
   bool isDistanceValid(const std::vector<int>& positions, int minDistance, int maxDistance);
    int WriteToFile(const string& filename);
    int openFile(const char* filename);
    // Function to get the size of the file
    off_t getFileSize(int fd);
    char* mapFile(int fd, off_t size);
    void searchPattern(char* mapped, off_t size, const char* pattern);
    void cleanUp(char* mapped, off_t size, int fd);
};

