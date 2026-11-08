#include "sdbg/sdbg.h"
#include "crispr_array/unitig_cycle_search.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include "original_precheck.h"

using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now()-t).count();
}
// Independent exact, unweighted BFS. All first arrivals have shortest depth.
class ExactBFS {
    const SDBG& g;
    std::vector<uint64_t> bits, dirty;
    struct State { uint64_t node; int depth; };
    std::vector<State> queue;
public:
    explicit ExactBFS(const SDBG& graph):g(graph),bits((g.size()+63)/64) {}
    bool has_cycle(uint64_t start, int limit) {
        for (auto word:dirty) bits[word]=0;
        dirty.clear(); queue.clear();
        auto mark=[&](uint64_t v) {
            auto& word=bits[v>>6];
            if (!word) dirty.push_back(v>>6);
            word |= uint64_t(1)<<(v&63);
        };
        queue.push_back({start,0});mark(start);
        const auto mult=g.EdgeMultiplicity(start);
        for (size_t i=0;i<queue.size();++i) {
            const auto state=queue[i];
            if (state.depth>=limit) continue;
            uint64_t neighbors[4];
            int degree=g.OutgoingEdges(state.node,neighbors);
            for (int j=0;j<degree;++j) {
                auto next=neighbors[j];
                if (!g.IsValidEdge(next)) continue;
                auto m=g.EdgeMultiplicity(next);
                if (!m || mult/m>500) continue;
                if (next==start) return true;
                if (!(bits[next>>6]>>(next&63)&1)) {
                    mark(next);queue.push_back({next,state.depth+1});
                }
            }
        }
        return false;
    }
};

int main(int argc,char**argv) {
    try {
        if (argc<3) {
            std::cerr<<"Usage: benchmark_unitig_search GRAPH_PREFIX RESULTS.csv "
                        "[MAX_CANDIDATES=0] [REPEATS=3] [LIMIT=77] [quick|all]\n";
            return 2;
        }
        size_t cap=argc>3?std::stoull(argv[3]):0;
        int repeats=argc>4?std::stoi(argv[4]):3;
        int limit=argc>5?std::stoi(argv[5]):77;
        const bool quick=argc>6 && std::string(argv[6])=="quick";
        if (repeats<1 || limit<1) throw std::runtime_error("positive repeats/limit required");
        omp_set_num_threads(1);
        if (!std::ifstream(std::string(argv[1])+".sdbg_info"))
            throw std::runtime_error("graph metadata not found");
        auto timer=Clock::now();
        SDBG g;g.LoadFromFile(argv[1]);
        const double load_seconds=elapsed(timer);
        OriginalPrecheck original(g);
        timer=Clock::now();
        // Match FindApproximateCRISPRArrays' exact pruning order and helpers.
        auto tips=original.CollectTips();
        original.InvalidateMultiplicityOneNodes();
        for (auto v:tips) original.RecursiveReduction(v);
        const double prune_seconds=elapsed(timer);
        timer=Clock::now();
        std::vector<uint64_t> all;
        for (uint64_t v=0;v<g.size();++v) {
            if (!g.IsValidEdge(v) || g.EdgeMultiplicity(v)<=20 || g.EdgeIndegree(v)<2)
                continue;
            uint64_t in[4];int n=g.IncomingEdges(v,in);
            if (std::find(in,in+n,v)==in+n) all.push_back(v);
        }
        std::vector<uint64_t> starts;
        if (cap && all.size()>cap) {
            for (size_t i=0;i<cap;++i) starts.push_back(all[i*all.size()/cap]);
        } else starts=all;
        const double scan_seconds=elapsed(timer);
        std::cout<<"Loaded "<<g.size()<<" SDBG entries, k="<<g.k()<<", "
                 <<all.size()<<" candidates, testing "<<starts.size()<<std::endl;
        ExactBFS reference(g);
        std::vector<bool> expected;
        timer=Clock::now();
        size_t positives=0;
        for (auto v:starts) {bool yes=reference.has_cycle(v,limit);expected.push_back(yes);positives+=yes;}
        const double oracle_seconds=elapsed(timer);
        std::vector<bool> original_answers;
        for (auto v:starts) original_answers.push_back(original.BidirectionalBFS(v,limit));
        std::vector<uint8_t> parallel_answers(starts.size());
        #pragma omp parallel num_threads(4)
        {
            mcaat::UnitigCycleSearch<SDBG> worker(g);
            #pragma omp for schedule(dynamic,256)
            for (size_t i=0;i<starts.size();++i)
                parallel_answers[i]=worker.has_cycle(starts[i],mcaat::upstream_precheck_bound(limit));
        }
        for (size_t i=0;i<starts.size();++i)
            if (bool(parallel_answers[i])!=original_answers[i])
                throw std::runtime_error("parallel cached/upstream disagreement");
        std::cout<<"Four-worker candidate comparison: identical\n";
        std::cout<<"Exact reference: "<<positives<<" accepted in "<<oracle_seconds<<" s\n";
        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("cannot write results");
        out<<"graph,entries,k,total_candidates,tested,limit,method,trial,seconds,accepted,"
               "mismatches,legacy_bound_only,cache_hits,cache_misses,built_steps,cached_steps,"
               "cache_entries,load_seconds,prune_seconds,scan_seconds,oracle_seconds\n";
        const char* methods[]={"upstream_bfs","bidir_no_jumps_exact","unitig_uncached_exact",
                               "unitig_cached_exact","unitig_cached_compatible"};
        const std::vector<int> modes=quick?std::vector<int>{0,4}:std::vector<int>{0,1,2,3,4};
        for (int rep=0;rep<repeats;++rep) for (size_t pos=0;pos<modes.size();++pos) {
            int mode=modes[(pos+rep)%modes.size()];
            // A fresh cache each trial: construction cost is included. Cache
            // entries may be reused across candidates, just as in production.
            mcaat::UnitigCycleSearch<SDBG> search(g,mode==2?0:16384);
            size_t accepted=0;
            std::vector<bool> answers;answers.reserve(starts.size());
            timer=Clock::now();
            for (auto v:starts) {
                bool yes=mode==0?original.BidirectionalBFS(v,limit):search.has_cycle(
                    v,mode==4?mcaat::upstream_precheck_bound(limit):limit,mode!=1);
                answers.push_back(yes);accepted+=yes;
            }
            const double seconds=elapsed(timer);
            size_t mismatches=0,bound_only=0;
            for (size_t i=0;i<starts.size();++i) if (answers[i]!=expected[i]) {
                ++mismatches;
                if ((mode==0 || mode==4) && answers[i] && reference.has_cycle(
                        starts[i],mcaat::upstream_precheck_bound(limit)))
                    ++bound_only;
            }
            if (mode==4 && answers!=original_answers)
                throw std::runtime_error("compatible optimized/upstream disagreement");
            if (mode!=0 && mode!=4 && mismatches)
                throw std::runtime_error("optimized/reference disagreement");
            auto s=search.statistics();
            out<<std::setprecision(10)<<'"'<<argv[1]<<"\","<<g.size()<<','<<g.k()<<','
               <<all.size()<<','<<starts.size()<<','<<limit<<','<<methods[mode]<<','<<rep+1<<','
               <<seconds<<','<<accepted<<','<<mismatches<<','<<bound_only<<','<<s.cache_hits<<','
               <<s.cache_misses<<','<<s.built_steps<<','<<s.cached_steps<<','<<search.cache_size()<<','
               <<load_seconds<<','<<prune_seconds<<','<<scan_seconds<<','<<oracle_seconds<<'\n';
            out.flush();
            std::cout<<methods[mode]<<" trial "<<rep+1<<": "<<seconds<<" s, "
                     <<accepted<<" accepted; oracle differences="<<mismatches
                     <<", cache hits="<<s.cache_hits<<std::endl;
        }
    } catch (const std::exception&e) {
        std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;
    }
}
