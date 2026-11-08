#!/usr/bin/env python3
"""Exercise the production enumeration body on a tiny branching graph.
The graph backend is a stub; a test watchdog makes the old non-exiting cap fail.
"""
import argparse, pathlib, subprocess, tempfile
p=argparse.ArgumentParser();p.add_argument('--cxx',default='c++');a=p.parse_args()
root=pathlib.Path(__file__).resolve().parents[1];source=(root/'src/crispr_array/cycle_finder.cpp').read_text()
start=source.index('vector<vector<uint64_t>> CycleFinder::FindCycle(');brace=source.index('{',start);depth=1;end=brace+1
while depth:
    depth+=(source[end]=='{')-(source[end]=='}');end+=1
method=source[start:end].replace('CycleFinder::','Harness::')
method=method.replace('while (!stack.empty()) {','while (!stack.empty()) {\n        if (++iterations > 1000) throw std::runtime_error("cycle cap failed to exit traversal");',1)
harness=r'''
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;
namespace phmap {template<class K,class V>using flat_hash_map=unordered_map<K,V>;template<class K>using flat_hash_set=unordered_set<K>;}
int omp_get_thread_num(){return 0;}
struct SmallNodeSet {
 uint64_t nodes[4];uint8_t count=0;
 void insert(uint64_t n){nodes[count++]=n;}
 void erase(uint64_t n){for(uint8_t i=0;i<count;++i)if(nodes[i]==n){nodes[i]=nodes[--count];return;}}
 uint64_t* begin(){return nodes;}uint64_t* end(){return nodes+count;}
};
struct Harness {
 struct Graph {int EdgeMultiplicity(uint64_t){return 10;}} graph;
 struct Settings {struct Limits{int cycle_min_length=1,cycle_max_length=5;}cycle_finder_settings;Graph* sdbg;}settings;
 size_t cluster_bounds=3,iterations=0;
 vector<vector<uint64_t>> per_thread_visited{vector<uint64_t>(8,0)};
 unordered_map<uint64_t,vector<uint64_t>> adjacency{{0,{1}},{1,{2,3,4}},{2,{0}},{3,{0}},{4,{0}}};
 Harness(){settings.sdbg=&graph;}
 void _GetOutgoings(uint64_t n,SmallNodeSet& out,int){for(auto v:adjacency[n])out.insert(v);}
 void _GetIncomings(uint64_t n,SmallNodeSet& out,int){for(auto& [u,vs]:adjacency)if(find(vs.begin(),vs.end(),n)!=vs.end())out.insert(u);}
 vector<vector<uint64_t>> FindCycle(uint64_t,vector<uint64_t>&,phmap::flat_hash_map<uint64_t,int>&,vector<SmallNodeSet>&,vector<int>&);
 vector<vector<uint64_t>> run(){vector<uint64_t> path{0};phmap::flat_hash_map<uint64_t,int> lock{{0,0}};vector<SmallNodeSet> stack(1);_GetOutgoings(0,stack[0],10);vector<int> back{settings.cycle_finder_settings.cycle_max_length};return FindCycle(0,path,lock,stack,back);}
};
'''
main=r'''
int main(){try {
 Harness limited;auto none=limited.run();
 if(!none.empty() || limited.iterations>20 || limited.per_thread_visited[0][0]!=0)throw runtime_error("capped tangle retained or marked");
 Harness allowed;allowed.cluster_bounds=4;auto cycles=allowed.run();
 if(cycles.size()!=3 || allowed.per_thread_visited[0][0]==0)throw runtime_error("below-cap cycles lost");
 for(size_t length:{46,47,48,99,100,101}) {
   Harness ring;ring.cluster_bounds=500;ring.adjacency.clear();
   ring.settings.cycle_finder_settings.cycle_min_length=47;ring.settings.cycle_finder_settings.cycle_max_length=100;
   for(size_t i=0;i<length;++i)ring.adjacency[i]={(i+1)%length};
   auto found=ring.run();bool expected=length>=47 && length<=100;
   if((found.size()==1)!=expected || (!found.empty() && found[0].size()!=length))throw runtime_error("incorrect inclusive cycle bound at "+to_string(length));
 }
 cout<<"Cycle limit checks passed: cap exits after "<<limited.iterations<<" steps; below cap retains 3 cycles; 47..100 endpoints verified\n";
 }catch(const exception& e){cerr<<e.what()<<'\n';return 1;}}
'''
with tempfile.TemporaryDirectory() as temp:
    path=pathlib.Path(temp);(path/'check.cpp').write_text(harness+method+main)
    subprocess.run([a.cxx,'-std=c++17','-O1',str(path/'check.cpp'),'-o',str(path/'check')],check=True)
    subprocess.run([str(path/'check')],check=True,timeout=10)
