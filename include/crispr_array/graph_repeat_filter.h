#ifndef MCAAT_GRAPH_REPEAT_FILTER_H
#define MCAAT_GRAPH_REPEAT_FILTER_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace mcaat::arrays {
using Node = uint64_t;
using Path = std::vector<Node>;
constexpr Node missing_node = std::numeric_limits<Node>::max();

inline double parse_repeat_support(const std::string& text) {
    size_t used=0;
    double value;
    try {value=std::stod(text,&used);}
    catch (...) {throw std::invalid_argument("repeat-support must be a finite fraction in (0, 1]");}
    if (used!=text.size() || !std::isfinite(value) || value<=0 || value>1)
        throw std::invalid_argument("repeat-support must be a finite fraction in (0, 1]");
    return value;
}

struct Options {
    size_t k = 23;
    double repeat_support = 0.85;
    size_t repeat_min = 23, repeat_max = 55, spacer_min = 23, spacer_max = 50;
    void validate() const {
        if (!std::isfinite(repeat_support) || repeat_support <= 0 || repeat_support > 1)
            throw std::invalid_argument("repeat-support must be a finite fraction in (0, 1]");
        if (!k || repeat_min > repeat_max || spacer_min > spacer_max || !spacer_min)
            throw std::invalid_argument("invalid graph repeat/spacer length bounds");
    }
};

struct Observation {
    Path cycle;                    // Rotated to the actual repeat start.
    size_t repeat_nodes = 0;
    std::string repeat, spacer;    // Observed strings, never synthetic consensus.
    bool supports_path = false;
};
struct System {
    Path repeat_path;
    std::vector<Path> candidate_cycles;
    std::vector<Observation> observations;
    std::vector<Node> sources;
    size_t support_count = 0;
    size_t alternate_repeat_paths = 0;
};
struct Statistics {
    size_t input_cycles = 0, invalid_cycles = 0, duplicate_cycles = 0;
    size_t input_groups = 0, groups_without_repeat = 0;
    size_t duplicate_systems = 0, merged_systems = 0;
};

inline std::string reverse_complement(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (auto p=sequence.rbegin();p!=sequence.rend();++p) {
        switch (*p) {
            case 'A': result+='T'; break; case 'C': result+='G'; break;
            case 'G': result+='C'; break; case 'T': result+='A'; break;
            default: throw std::invalid_argument("non-ACGT graph label");
        }
    }
    return result;
}
inline Path rotate(Path path, size_t offset) {
    if (!path.empty()) std::rotate(path.begin(),path.begin()+offset,path.end());
    return path;
}
inline Path cycle_key(Path path) {
    if (!path.empty()) path=rotate(path,std::min_element(path.begin(),path.end())-path.begin());
    return path;
}
inline bool contains_path(const Path& cycle, size_t start, const Path& path) {
    if (path.size()>cycle.size()) return false;
    for (size_t i=0;i<path.size();++i)
        if (cycle[(start+i)%cycle.size()]!=path[i]) return false;
    return true;
}

struct Parts {Path repeat,spacer,repeat_to_spacer,spacer_to_repeat;};
inline Parts partition(const Path& cycle,size_t repeat_nodes,size_t k) {
    if (!repeat_nodes || repeat_nodes>cycle.size() || !k || repeat_nodes+k-1>=cycle.size())
        throw std::runtime_error("invalid cycle partition");
    Parts p;
    const size_t r=repeat_nodes+k-1,m=cycle.size();
    for (size_t i=0;i<m;++i) {
        if (i<repeat_nodes) p.repeat.push_back(cycle[i]);
        if (i>=r && i+k<=m) p.spacer.push_back(cycle[i]);
        if (i<r && i+k>r) p.repeat_to_spacer.push_back(cycle[i]);
        if (i+k>m) p.spacer_to_repeat.push_back(cycle[i]);
    }
    return p;
}

// Graph provides label(vertex), adjacent(u,v), reverse_node(vertex). Vertices
// identify complete k-mers, not individual outgoing SDBG edge handles.
template<class Graph>
class RepeatFilter {
public:
    RepeatFilter(Graph& graph, Options options):g_(graph),o_(options) {o_.validate();}
    Statistics stats;

    std::vector<System> run(const std::map<Node,std::vector<Path>>& input) {
        stats={};
        std::vector<System> systems;
        for (const auto& [source,raw]:input) {
            ++stats.input_groups;
            std::set<Path> distinct;
            for (Path cycle:raw) {
                ++stats.input_cycles;
                if (cycle.size()>1 && cycle.front()==cycle.back()) cycle.pop_back();
                if (!valid_cycle(cycle)) {++stats.invalid_cycles;continue;}
                if (!distinct.insert(cycle_key(std::move(cycle))).second) ++stats.duplicate_cycles;
            }
            System system;
            if (infer({distinct.begin(),distinct.end()},system)) {
                system.sources={source};orient(system);systems.push_back(std::move(system));
            } else ++stats.groups_without_repeat;
        }
        deduplicate(systems);
        std::sort(systems.begin(),systems.end(),[](const System& a,const System& b) {
            return a.repeat_path!=b.repeat_path ? a.repeat_path<b.repeat_path
                                               : observed_cycles(a)<observed_cycles(b);
        });
        return systems;
    }

    std::string spell(const Path& path) {
        if (path.empty()) return {};
        std::string sequence=g_.label(path.front());
        if (sequence.size()!=o_.k) throw std::runtime_error("wrong graph k-mer length");
        for (size_t i=1;i<path.size();++i) {
            const auto next=g_.label(path[i]);
            if (next.size()!=o_.k || sequence.compare(sequence.size()-o_.k+1,o_.k-1,next,0,o_.k-1)!=0)
                throw std::runtime_error("non-overlapping graph path");
            sequence+=next.back();
        }
        return sequence;
    }

private:
    Graph& g_;
    Options o_;
    struct Occurrence {size_t cycle,position;};

    bool valid_cycle(const Path& cycle) {
        if (cycle.size()<2 || cycle.size()<o_.k) return false;
        std::set<Node> visited;
        for (size_t i=0;i<cycle.size();++i) {
            if (!visited.insert(cycle[i]).second || !g_.adjacent(cycle[i],cycle[(i+1)%cycle.size()]))
                return false;
            auto a=g_.label(cycle[i]),b=g_.label(cycle[(i+1)%cycle.size()]);
            if (a.size()!=o_.k || b.size()!=o_.k || a.substr(1)!=b.substr(0,o_.k-1)) return false;
        }
        return true;
    }
    static std::set<Path> observed_cycles(const System& system) {
        std::set<Path> cycles;
        for (const auto& row:system.observations) cycles.insert(cycle_key(row.cycle));
        return cycles;
    }

    void place(System& system) {
        system.observations.clear();system.support_count=0;
        const auto& repeat=system.repeat_path;
        for (const auto& cycle:system.candidate_cycles) {
            auto begin=std::find(cycle.begin(),cycle.end(),repeat.front());
            auto end=std::find(cycle.begin(),cycle.end(),repeat.back());
            if (begin==cycle.end() || end==cycle.end()) continue;
            const size_t offset=begin-cycle.begin();
            const size_t nodes=(size_t(end-cycle.begin())+cycle.size()-offset)%cycle.size()+1;
            const size_t repeat_length=o_.k+nodes-1;
            if (repeat_length>=cycle.size() || repeat_length<o_.repeat_min || repeat_length>o_.repeat_max)
                continue;
            const size_t spacer_length=cycle.size()-repeat_length;
            if (spacer_length<o_.spacer_min || spacer_length>o_.spacer_max) continue;
            Path rotated=rotate(cycle,offset);
            const bool supported=contains_path(cycle,offset,repeat);
            const std::string sequence=spell(rotated);
            system.observations.push_back({std::move(rotated),nodes,
                sequence.substr(0,repeat_length),sequence.substr(repeat_length,spacer_length),supported});
            system.support_count+=supported;
        }
        std::sort(system.observations.begin(),system.observations.end(),[](const auto& a,const auto& b) {
            if (a.spacer!=b.spacer) return a.spacer<b.spacer;
            if (a.supports_path!=b.supports_path) return a.supports_path>b.supports_path;
            return a.cycle<b.cycle;
        });
    }

    bool infer(const std::vector<Path>& cycles,System& result) {
        if (cycles.size()<2) return false;
        const size_t required=std::max<size_t>(2,size_t(std::ceil(
            o_.repeat_support*double(cycles.size())-1e-12)));
        std::map<Node,std::vector<Occurrence>> index;
        for (size_t c=0;c<cycles.size();++c)
            for (size_t p=0;p<cycles[c].size();++p) index[cycles[c][p]].push_back({c,p});
        std::set<Path> paths;
        // Enumerate observed frequent substrings, intersecting cycle support at
        // every extension. Per-node majority alone could synthesize a path that
        // no sufficient set of cycles actually contains.
        auto extend=[&](auto&& self,Path& path,const std::vector<Occurrence>& support)->void {
            std::map<Node,std::vector<Occurrence>> next;
            for (auto p:support) {
                const auto& c=cycles[p.cycle];
                if (path.size()<c.size()) next[c[(p.position+path.size())%c.size()]].push_back(p);
            }
            bool extended=false;
            for (const auto& [node,subset]:next) if (subset.size()>=required &&
                    std::find(path.begin(),path.end(),node)==path.end()) {
                extended=true;
                // Once over the maximum, all extensions are invalid. Do not
                // truncate a longer conserved region into an artificial repeat.
                if (o_.k+path.size()<=o_.repeat_max) {
                    path.push_back(node);self(self,path,subset);path.pop_back();
                }
            }
            if (extended) return;
            const size_t length=o_.k+path.size()-1;
            if (length<o_.repeat_min || length>o_.repeat_max) return;
            std::map<Node,size_t> previous;
            for (auto p:support) {
                const auto& c=cycles[p.cycle];
                ++previous[c[(p.position+c.size()-1)%c.size()]];
            }
            for (auto [node,count]:previous)
                if (count>=required && std::find(path.begin(),path.end(),node)==path.end()) return;
            paths.insert(path);
        };
        for (const auto& [node,support]:index) if (support.size()>=required) {
            Path path{node};extend(extend,path,support);
        }
        bool found=false;
        size_t alternatives=0;
        for (const auto& path:paths) {
            System candidate;candidate.repeat_path=path;candidate.candidate_cycles=cycles;
            place(candidate);
            if (candidate.support_count<required) continue;
            std::set<std::string> spacers;
            std::vector<size_t> lengths;
            for (const auto& row:candidate.observations)
                if (spacers.insert(row.spacer).second) lengths.push_back(row.spacer.size());
            if (spacers.size()<2) continue;
            std::sort(lengths.begin(),lengths.end());
            const double median=lengths.size()%2 ? lengths[lengths.size()/2] :
                (double(lengths[lengths.size()/2-1])+lengths[lengths.size()/2])/2;
            if (median/double(o_.k+path.size()-1)<0.5) continue;
            ++alternatives;
            if (!found || path.size()>result.repeat_path.size() ||
                (path.size()==result.repeat_path.size() && candidate.support_count>result.support_count) ||
                (path.size()==result.repeat_path.size() && candidate.support_count==result.support_count && path<result.repeat_path)) {
                result=std::move(candidate);found=true;
            }
        }
        if (found) result.alternate_repeat_paths=alternatives-1;
        return found;
    }

    Path reversed(const Path& path) {
        Path result;
        for (auto i=path.rbegin();i!=path.rend();++i) {
            Node node=g_.reverse_node(*i);
            if (node==missing_node) return {};
            result.push_back(node);
        }
        return result;
    }
    Path strand_key(const Path& path) {
        auto key=cycle_key(path),rc=reversed(path);
        return rc.empty()?key:std::min(key,cycle_key(std::move(rc)));
    }
    bool flip(System& system) {
        auto repeat=reversed(system.repeat_path);
        if (repeat.empty()) return false;
        std::vector<Path> cycles;
        for (const auto& c:system.candidate_cycles) {
            auto r=reversed(c);
            if (r.empty() || !valid_cycle(r)) return false;
            cycles.push_back(cycle_key(std::move(r)));
        }
        std::sort(cycles.begin(),cycles.end());
        system.repeat_path=std::move(repeat);system.candidate_cycles=std::move(cycles);place(system);
        return true;
    }
    void orient(System& system) {
        const auto sequence=spell(system.repeat_path);
        auto repeat=reversed(system.repeat_path);
        if (repeat.empty()) return;
        const auto rc=reverse_complement(sequence);
        if (rc>sequence || (rc==sequence && repeat>=system.repeat_path)) return;
        flip(system);
    }

    std::pair<std::set<Node>,std::set<Node>> roles(const System& s) {
        std::set<Node> repeats,spacers;
        for (const auto& row:s.observations) {
            repeats.insert(row.cycle.begin(),row.cycle.begin()+row.repeat_nodes);
            // K-mers crossing either boundary are deliberately not pure spacers.
            const size_t r=o_.k+row.repeat_nodes-1;
            for (size_t i=r;i+o_.k<=row.cycle.size();++i) spacers.insert(row.cycle[i]);
        }
        return {repeats,spacers};
    }
    static bool intersects(const std::set<Node>& a,const std::set<Node>& b) {
        for (Node node:a) if (b.count(node)) return true;
        return false;
    }
    bool compatible(const System& a,const System& b) {
        auto ar=roles(a),br=roles(b);
        return intersects(ar.first,br.first) && !intersects(ar.first,br.second) &&
               !intersects(br.first,ar.second);
    }
    static bool subset(const std::set<Path>& a,const std::set<Path>& b) {
        return std::includes(b.begin(),b.end(),a.begin(),a.end());
    }
    static void add_sources(System& keep,const System& other) {
        keep.sources.insert(keep.sources.end(),other.sources.begin(),other.sources.end());
        std::sort(keep.sources.begin(),keep.sources.end());
        keep.sources.erase(std::unique(keep.sources.begin(),keep.sources.end()),keep.sources.end());
    }
    void deduplicate(std::vector<System>& systems) {
        // Only compare systems sharing an entire observed cycle. A common
        // repeat alone, or a spacer shared with another repeat, is insufficient.
        bool changed=true;
        while (changed) {
            changed=false;
            std::map<Path,std::vector<size_t>> owners;
            for (size_t i=0;i<systems.size();++i)
                for (const auto& c:observed_cycles(systems[i])) owners[strand_key(c)].push_back(i);
            std::set<std::pair<size_t,size_t>> pairs;
            for (const auto& [cycle,ids]:owners)
                for (size_t i=0;i<ids.size();++i) for (size_t j=i+1;j<ids.size();++j)
                    pairs.emplace(ids[i],ids[j]);
            for (auto [i,j]:pairs) {
                if (i==j) continue;
                auto& a=systems[i];auto b=systems[j];
                // Slightly different repeat boundaries can select opposite
                // canonical strands. Align by shared full cycles before roles.
                const auto ac=observed_cycles(a);
                auto bc=observed_cycles(b);
                bool common=false;
                for (const auto& c:ac) common|=bc.count(c)>0;
                if (!common) {if (!flip(b)) continue;bc=observed_cycles(b);}
                if (!compatible(a,b)) continue;
                if (subset(ac,bc) || subset(bc,ac)) {
                    if (subset(ac,bc) && !subset(bc,ac)) {add_sources(b,a);a=std::move(b);}
                    else add_sources(a,b);
                    orient(a);systems.erase(systems.begin()+j);++stats.duplicate_systems;changed=true;break;
                }
                std::set<Path> united(a.candidate_cycles.begin(),a.candidate_cycles.end());
                united.insert(b.candidate_cycles.begin(),b.candidate_cycles.end());
                System merged;
                if (!infer({united.begin(),united.end()},merged)) continue;
                const auto mc=observed_cycles(merged);
                if (!compatible(merged,a) || !compatible(merged,b) || !subset(ac,mc) || !subset(bc,mc)) continue;
                add_sources(merged,a);add_sources(merged,b);orient(merged);a=std::move(merged);
                systems.erase(systems.begin()+j);++stats.merged_systems;changed=true;break;
            }
        }
    }
};
} // namespace mcaat::arrays
#endif
