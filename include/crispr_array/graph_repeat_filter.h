#ifndef MCAAT_GRAPH_REPEAT_FILTER_H
#define MCAAT_GRAPH_REPEAT_FILTER_H
#include "crispr_array/candidate_quality.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mcaat::arrays {
using Node = uint64_t;
using Path = std::vector<Node>;
using Progress = std::function<void(const char*,size_t,size_t)>;
constexpr Node missing_node = std::numeric_limits<Node>::max();

inline double parse_fraction(const std::string& text,const std::string& name) {
    size_t used=0;
    double value;
    try {value=std::stod(text,&used);}
    catch (...) {throw std::invalid_argument(name+" must be a finite fraction in (0, 1]");}
    if (used!=text.size() || !std::isfinite(value) || value<=0 || value>1)
        throw std::invalid_argument(name+" must be a finite fraction in (0, 1]");
    return value;
}
inline double parse_repeat_support(const std::string& text) {return parse_fraction(text,"repeat-support");}
inline bool parse_filter_mode(const std::string& text) {
    if (text=="conservative") return false;
    if (text=="strict") return true;
    throw std::invalid_argument("array-filter-mode must be conservative or strict");
}
inline size_t parse_max_cycles(const std::string& text) {
    size_t used=0;unsigned long long value=0;
    try {value=std::stoull(text,&used);} catch (...) {used=0;}
    if (text.empty() || text.find_first_not_of("0123456789")!=std::string::npos ||
        used!=text.size() || value<3 || value>std::numeric_limits<int>::max())
        throw std::invalid_argument("max-cycles must be an integer from 3 to INT_MAX");
    return size_t(value);
}

struct Options {
    size_t k = 23;
    double repeat_support = 0.85;
    size_t repeat_min = 23, repeat_max = 55, spacer_min = 23, spacer_max = 50;
    QualityOptions quality;
    size_t cycle_min=1,cycle_max=std::numeric_limits<int>::max();
    void validate() const {
        if (!std::isfinite(repeat_support) || repeat_support <= 0 || repeat_support > 1)
            throw std::invalid_argument("repeat-support must be a finite fraction in (0, 1]");
        if (!k || repeat_min > repeat_max || spacer_min > spacer_max || !spacer_min)
            throw std::invalid_argument("invalid graph repeat/spacer length bounds");
        quality.validate();
        if (!cycle_min || cycle_min>cycle_max) throw std::invalid_argument("invalid cycle length range");
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
    size_t candidate_id = 0, candidate_count = 0;
    Quality quality;
};
struct CandidateDecision {
    System original;
    std::string status;
    std::set<size_t> competing_candidates;
    std::set<Path> retained_paths;
    size_t removed_cycles=0, duplicate_strand_cycles=0, out_of_range_cycles=0;
};
struct Statistics {
    size_t input_cycles = 0, invalid_cycles = 0, duplicate_cycles = 0;
    size_t input_groups = 0, groups_without_repeat = 0;
    size_t duplicate_systems = 0, merged_systems = 0;
    size_t dedup_pair_checks = 0, dedup_index_builds = 0;
    size_t quality_rejected=0, overlap_rejected=0, shared_cycles_removed=0, duplicate_strand_cycles=0;
    size_t out_of_range_cycles=0, curation_out_of_range_cycles=0;
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
inline std::string circular_sequence_key(const std::string& sequence) {
    if (sequence.empty()) return sequence;
    std::string best=sequence;
    for (const auto& direction:{sequence,reverse_complement(sequence)})
        for (size_t i=0;i<direction.size();++i) {
            if (direction[i]>best.front()) continue;
            auto rotated=direction.substr(i)+direction.substr(0,i);
            if (rotated<best) best=std::move(rotated);
        }
    return best;
}
inline bool contains_path(const Path& cycle, size_t start, const Path& path) {
    if (path.size()>cycle.size()) return false;
    for (size_t i=0;i<path.size();++i)
        if (cycle[(start+i)%cycle.size()]!=path[i]) return false;
    return true;
}

inline std::string path_text(const Path& path) {
    if (path.empty()) return "-";
    std::ostringstream out;
    for (size_t i=0;i<path.size();++i) {if(i) out<<',';out<<path[i];}
    return out.str();
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
    std::vector<CandidateDecision> decisions;

    std::vector<System> run(const std::map<Node,std::vector<Path>>& input,const Progress& progress={}) {
        stats={};
        decisions.clear();
        validity_cache_.clear();
        std::vector<System> systems;
        if (progress) progress("Repeat inference",0,input.size());
        for (const auto& [source,raw]:input) {
            ++stats.input_groups;
            std::set<Path> distinct;
            for (Path cycle:raw) {
                ++stats.input_cycles;
                if (cycle.size()>1 && cycle.front()==cycle.back()) cycle.pop_back();
                if (cycle.size()<o_.cycle_min || cycle.size()>o_.cycle_max) {++stats.out_of_range_cycles;continue;}
                cycle=cycle_key(std::move(cycle));
                // Repeated saved records and rotated rediscoveries must not
                // revalidate every graph edge before being recognized as copies.
                if (distinct.count(cycle)) {++stats.duplicate_cycles;continue;}
                if (!valid_cycle(cycle)) {++stats.invalid_cycles;continue;}
                distinct.insert(std::move(cycle));
            }
            System system;
            if (infer({distinct.begin(),distinct.end()},system)) {
                system.sources={source};orient(system);systems.push_back(std::move(system));
            } else ++stats.groups_without_repeat;
            if (progress) progress("Repeat inference",stats.input_groups,input.size());
        }
        deduplicate(systems,progress);
        std::sort(systems.begin(),systems.end(),[](const System& a,const System& b) {
            return a.repeat_path!=b.repeat_path ? a.repeat_path<b.repeat_path
                                               : observed_cycles(a)<observed_cycles(b);
        });
        return curate(std::move(systems),progress);
    }

    // Final quality and evidence-ownership pass. Also usable by report-replay
    // checks: candidate_count must retain the original denominator, including
    // unresolved cycles when only assigned observations are available.
    std::vector<System> curate(std::vector<System> systems,const Progress& progress={}) {
        decisions.clear();
        stats.quality_rejected=stats.overlap_rejected=stats.shared_cycles_removed=stats.duplicate_strand_cycles=0;
        stats.curation_out_of_range_cycles=0;
        std::map<Path,std::string> keys;
        auto sequence_key=[&](const Path& cycle)->const std::string& {
            auto key=cycle_key(cycle);auto it=keys.find(key);
            if (it!=keys.end()) return it->second;
            return keys.emplace(std::move(key),circular_sequence_key(spell(cycle).substr(0,cycle.size()))).first->second;
        };
        auto assess=[&](System& system) {
            std::vector<std::string> spacers;
            for (const auto& row:system.observations) spacers.push_back(row.spacer);
            system.quality=assess_quality(spell(system.repeat_path),std::move(spacers),system.candidate_count,o_.quality);
        };
        if (progress) progress("Candidate quality",0,systems.size());
        for (size_t i=0;i<systems.size();++i) {
            auto& system=systems[i];system.candidate_id=i+1;
            if (!system.candidate_count) system.candidate_count=system.candidate_cycles.size();
            if (system.candidate_count<system.candidate_cycles.size() || system.support_count>system.candidate_count ||
                system.observations.size()>system.candidate_cycles.size())
                throw std::runtime_error("inconsistent candidate cycle counts");
            for (const auto& row:system.observations)
                keys.emplace(cycle_key(row.cycle),circular_sequence_key(row.repeat+row.spacer));
            assess(system);
            if (progress) progress("Candidate quality",i+1,systems.size());
        }
        std::vector<size_t> order(systems.size());std::iota(order.begin(),order.end(),0);
        // Prefer candidates with more distinct spacer families, then better
        // sequence separation. Variant counts are deliberately not abundance.
        if (o_.quality.strict) std::sort(order.begin(),order.end(),[&](size_t i,size_t j) {
            const auto& a=systems[i];const auto& b=systems[j];
            if (a.quality.reason.empty()!=b.quality.reason.empty()) return a.quality.reason.empty();
            if (a.quality.families.size()!=b.quality.families.size()) return a.quality.families.size()>b.quality.families.size();
            if (a.quality.spacer_similarity!=b.quality.spacer_similarity) return a.quality.spacer_similarity<b.quality.spacer_similarity;
            if (a.repeat_path.size()!=b.repeat_path.size()) return a.repeat_path.size()>b.repeat_path.size();
            return i<j;
        });
        std::map<std::string,size_t> owners;
        std::vector<System> accepted;
        decisions.resize(systems.size());
        if (progress) progress("Resolve shared cycles",0,systems.size());
        size_t done=0;
        for (size_t index:order) {
            auto& decision=decisions[index];decision.original=std::move(systems[index]);
            const auto& original=decision.original;
            if (!original.quality.reason.empty()) {
                decision.status=original.quality.reason;++stats.quality_rejected;
            } else {
                System keep=original;
                keep.candidate_cycles.clear();keep.observations.clear();keep.support_count=0;
                std::map<std::string,size_t> remaining,assigned;
                for (const auto& cycle:original.candidate_cycles) {
                    if (cycle.size()<o_.cycle_min || cycle.size()>o_.cycle_max) {++decision.out_of_range_cycles;continue;}
                    const auto& key=sequence_key(cycle);auto owned=owners.find(key);
                    if (!o_.quality.strict || owned==owners.end()) {
                        if (remaining.emplace(key,keep.candidate_cycles.size()).second) keep.candidate_cycles.push_back(cycle);
                        else ++decision.duplicate_strand_cycles;
                    }
                    else {decision.competing_candidates.insert(owned->second);++decision.removed_cycles;}
                }
                keep.candidate_count-=decision.removed_cycles+decision.duplicate_strand_cycles+decision.out_of_range_cycles;
                for (const auto& row:original.observations) {
                    if (row.cycle.size()<o_.cycle_min || row.cycle.size()>o_.cycle_max) continue;
                    const auto& key=sequence_key(row.cycle);
                    if (!o_.quality.strict || !owners.count(key)) {
                        auto [position,fresh]=assigned.emplace(key,keep.observations.size());
                        if (fresh) keep.observations.push_back(row);
                        else if (row.supports_path && !keep.observations[position->second].supports_path)
                            keep.observations[position->second]=row;
                    }
                }
                for (const auto& row:keep.observations) {
                    keep.support_count+=row.supports_path;
                    keep.candidate_cycles.at(remaining.at(sequence_key(row.cycle)))=cycle_key(row.cycle);
                }
                const size_t required=std::max<size_t>(2,size_t(std::ceil(o_.repeat_support*keep.candidate_count-1e-12)));
                if (keep.observations.size()<2) decision.status=decision.out_of_range_cycles?"cycle_length_range":
                    (o_.quality.strict?"shared_cycle_subset":"fewer_than_two_cycles");
                else if (keep.support_count<required) decision.status=o_.quality.strict?"repeat_support_after_overlap":"repeat_support_after_deduplication";
                else {
                    assess(keep);decision.status=keep.quality.reason;
                }
                if (decision.status.empty()) {
                    decision.status=decision.removed_cycles?"accepted_without_shared_cycles":
                        (decision.duplicate_strand_cycles?"accepted_without_strand_duplicates":"accepted");
                    for (const auto& row:keep.observations) {
                        const auto& key=sequence_key(row.cycle);
                        if (!owners.emplace(key,keep.candidate_id).second && o_.quality.strict)
                            throw std::runtime_error("internal error: cycle assigned to multiple accepted arrays");
                        decision.retained_paths.insert(cycle_key(row.cycle));
                    }
                    accepted.push_back(std::move(keep));
                } else ++stats.overlap_rejected;
                stats.shared_cycles_removed+=decision.removed_cycles;
                stats.duplicate_strand_cycles+=decision.duplicate_strand_cycles;
                stats.curation_out_of_range_cycles+=decision.out_of_range_cycles;
            }
            if (progress) progress("Resolve shared cycles",++done,decisions.size());
        }
        std::sort(accepted.begin(),accepted.end(),[](const auto& a,const auto& b) {return a.candidate_id<b.candidate_id;});
        if (!o_.quality.strict) {
            std::map<std::string,std::vector<size_t>> shared;
            for (size_t i=0;i<accepted.size();++i)
                for (const auto& row:accepted[i].observations) shared[sequence_key(row.cycle)].push_back(i);
            for (const auto& [key,indexes]:shared) if (indexes.size()>1) {
                for (size_t i:indexes) {
                    auto& s=accepted[i];auto& d=decisions.at(s.candidate_id-1);
                    if (std::find(s.quality.flags.begin(),s.quality.flags.end(),"shared_cycle_evidence")==s.quality.flags.end())
                        s.quality.flags.push_back("shared_cycle_evidence");
                    d.status="accepted_with_shared_evidence";
                    // Complete membership is serialized once per shared cycle.
                    // Do not materialize a quadratic pairwise candidate clique.
                }
            }
        }
        return accepted;
    }

    std::string spell(const Path& path) {
        if (path.empty()) return {};
        std::string sequence=g_.label(path.front());
        if (sequence.size()!=o_.k) throw std::runtime_error("wrong graph k-mer length");
        for (size_t i=1;i<path.size();++i) {
            const auto& next=g_.label(path[i]);
            if (next.size()!=o_.k || sequence.compare(sequence.size()-o_.k+1,o_.k-1,next,0,o_.k-1)!=0)
                throw std::runtime_error("non-overlapping graph path");
            sequence+=next.back();
        }
        return sequence;
    }

private:
    Graph& g_;
    Options o_;
    std::map<Path,bool> validity_cache_;
    static constexpr size_t validity_cache_limit=65536;
    struct Occurrence {size_t cycle,position;};

    bool valid_cycle(const Path& cycle) {
        auto key=cycle_key(cycle);
        auto cached=validity_cache_.find(key);
        if (cached!=validity_cache_.end()) return cached->second;
        const bool valid=check_cycle(cycle);
        if (validity_cache_.size()>=validity_cache_limit) validity_cache_.clear();
        validity_cache_.emplace(std::move(key),valid);
        return valid;
    }
    bool check_cycle(const Path& cycle) {
        if (cycle.size()<2 || cycle.size()<o_.k) return false;
        std::set<Node> visited;
        for (size_t i=0;i<cycle.size();++i) {
            if (!visited.insert(cycle[i]).second || !g_.adjacent(cycle[i],cycle[(i+1)%cycle.size()]))
                return false;
            const auto& a=g_.label(cycle[i]);
            const auto& b=g_.label(cycle[(i+1)%cycle.size()]);
            if (a.size()!=o_.k || b.size()!=o_.k || a.compare(1,o_.k-1,b,0,o_.k-1)!=0) return false;
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
            System candidate;candidate.repeat_path=path;candidate.candidate_cycles=cycles;candidate.candidate_count=cycles.size();
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
        result.reserve(path.size());
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

    using Roles = std::pair<std::set<Node>,std::set<Node>>;
    Roles roles(const System& s) {
        std::set<Node> repeats,spacers;
        for (const auto& row:s.observations) {
            repeats.insert(row.cycle.begin(),row.cycle.begin()+row.repeat_nodes);
            // K-mers crossing either boundary are deliberately not pure spacers.
            const size_t r=o_.k+row.repeat_nodes-1;
            for (size_t i=r;i+o_.k<=row.cycle.size();++i) spacers.insert(row.cycle[i]);
        }
        return {repeats,spacers};
    }
    template<class T>
    static bool intersects(const std::set<T>& a,const std::set<T>& b) {
        auto i=a.begin(),j=b.begin();
        while (i!=a.end() && j!=b.end()) {
            if (*i<*j) ++i;
            else if (*j<*i) ++j;
            else return true;
        }
        return false;
    }
    static bool compatible(const Roles& ar,const Roles& br) {
        return intersects(ar.first,br.first) && !intersects(ar.first,br.second) &&
               !intersects(br.first,ar.second);
    }
    static bool subset(const std::set<Path>& a,const std::set<Path>& b) {
        return std::includes(b.begin(),b.end(),a.begin(),a.end());
    }
    static void add_sources(System& keep,const System& other) {
        const bool append_sorted=keep.sources.empty() || other.sources.empty() ||
                                 keep.sources.back()<other.sources.front();
        const size_t middle=keep.sources.size();
        keep.sources.insert(keep.sources.end(),other.sources.begin(),other.sources.end());
        if (!append_sorted) {
            std::inplace_merge(keep.sources.begin(),keep.sources.begin()+middle,keep.sources.end());
            keep.sources.erase(std::unique(keep.sources.begin(),keep.sources.end()),keep.sources.end());
        }
    }
    void deduplicate(std::vector<System>& systems,const Progress& progress) {
        // Only compare systems sharing an entire observed cycle. A common
        // repeat alone, or a spacer shared with another repeat, is insufficient.
        // Stable indices preserve the original lexicographic merge order. Keep
        // the index across removals; do not rebuild every system or materialize
        // the quadratic clique of pairs belonging to one common cycle.
        struct Metadata {
            std::set<Path> observed, strands;
            Roles parts;
        };
        std::vector<Metadata> metadata(systems.size());
        std::vector<bool> alive(systems.size(),true);
        std::map<Path,std::set<size_t>> owners;
        std::set<size_t> pending;
        auto index=[&](size_t i) {
            auto& m=metadata[i];
            m.observed=observed_cycles(systems[i]);m.parts=roles(systems[i]);m.strands.clear();
            for (const auto& c:m.observed) m.strands.insert(strand_key(c));
            for (const auto& c:m.strands) owners[c].insert(i);
            ++stats.dedup_index_builds;
        };
        auto unindex=[&](size_t i) {
            for (const auto& c:metadata[i].strands) {
                auto it=owners.find(c);
                it->second.erase(i);
                if (it->second.empty()) owners.erase(it);
            }
        };
        auto retire=[&](size_t j) {
            unindex(j);alive[j]=false;pending.erase(j);
            metadata[j]={};systems[j]={};
        };
        if (progress) progress("Indexing candidate systems",0,systems.size());
        for (size_t i=0;i<systems.size();++i) {
            index(i);pending.insert(i);
            if (progress) progress("Indexing candidate systems",i+1,systems.size());
        }
        if (progress) progress("Duplicate comparisons",0,0);
        while (!pending.empty()) {
            const size_t i=*pending.begin();pending.erase(pending.begin());
            std::set<size_t> neighbors;
            for (const auto& c:metadata[i].strands) {
                const auto& ids=owners.at(c);
                neighbors.insert(ids.upper_bound(i),ids.end());
            }
            for (size_t j:neighbors) {
                if (!alive[j]) continue;
                ++stats.dedup_pair_checks;
                if (progress) progress("Duplicate comparisons",stats.dedup_pair_checks,0);
                auto& a=systems[i];
                const System* b=&systems[j];
                const auto& ac=metadata[i].observed;
                const auto& ar=metadata[i].parts;
                const auto* bc=&metadata[j].observed;
                const auto* br=&metadata[j].parts;
                System flipped;
                std::set<Path> flipped_cycles;
                Roles flipped_parts;
                // Slightly different repeat boundaries can select opposite
                // canonical strands. Align by shared full cycles before roles.
                if (!intersects(ac,*bc)) {
                    flipped=*b;
                    if (!flip(flipped)) continue;
                    b=&flipped;flipped_cycles=observed_cycles(flipped);flipped_parts=roles(flipped);
                    bc=&flipped_cycles;br=&flipped_parts;
                }
                if (!compatible(ar,*br)) continue;
                const bool a_in_b=subset(ac,*bc),b_in_a=subset(*bc,ac);
                if (b_in_a) {
                    // Only provenance changes. All earlier rejected pairs are
                    // still rejected, so continue without rebuilding even i.
                    add_sources(a,*b);retire(j);++stats.duplicate_systems;
                    continue;
                }
                System replacement;
                if (a_in_b) {
                    replacement=*b;add_sources(replacement,a);
                    ++stats.duplicate_systems;
                } else {
                    std::set<Path> united(a.candidate_cycles.begin(),a.candidate_cycles.end());
                    united.insert(b->candidate_cycles.begin(),b->candidate_cycles.end());
                    if (!infer({united.begin(),united.end()},replacement)) continue;
                    const auto mc=observed_cycles(replacement);
                    const auto mr=roles(replacement);
                    if (!compatible(mr,ar) || !compatible(mr,*br) || !subset(ac,mc) || !subset(*bc,mc)) continue;
                    add_sources(replacement,a);add_sources(replacement,*b);++stats.merged_systems;
                }
                orient(replacement);unindex(i);retire(j);a=std::move(replacement);index(i);
                // A changed repeat boundary or cycle set can make an earlier
                // rejected pair eligible. Revisit those earlier owners first.
                pending.insert(i);
                for (const auto& c:metadata[i].strands) {
                    const auto& ids=owners.at(c);
                    pending.insert(ids.begin(),ids.lower_bound(i));
                }
                break;
            }
        }
        size_t kept=0;
        for (size_t i=0;i<systems.size();++i) if (alive[i]) {
            if (kept!=i) systems[kept]=std::move(systems[i]);
            ++kept;
        }
        systems.resize(kept);
    }
};
} // namespace mcaat::arrays
#endif
