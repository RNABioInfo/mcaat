#ifndef MCAAT_PROTO_SPACER_BOUNDARIES_H
#define MCAAT_PROTO_SPACER_BOUNDARIES_H
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace mcaat::proto {
struct BoundaryOptions {
    bool enabled=true;
    double support=.85;
    size_t min_spacers=5,min_bases=6;
    void validate() const {
        if(!std::isfinite(support) || support<=.5 || support>1 || min_spacers<2 || !min_bases)
            throw std::invalid_argument("boundary support must be in (0.5,1], minimum spacers >=2 and minimum bases >0");
    }
};
struct BoundaryInput {std::string group,repeat,spacer;};
struct BoundaryResult {
    std::string repeat,spacer,status="unchanged";
    size_t left=0,right=0;
};
inline std::string boundary_rc(std::string s) {
    std::reverse(s.begin(),s.end());
    for(char& c:s) {
        switch(c){case 'A':c='T';break;case 'T':c='A';break;case 'C':c='G';break;case 'G':c='C';break;
                  default:throw std::invalid_argument("boundary input must be ACGT");}
    }
    return s;
}
inline std::string conserved_terminal(const std::set<std::string>& spacers,bool right,const BoundaryOptions& options) {
    if(spacers.size()<options.min_spacers)return {};
    size_t shortest=spacers.begin()->size();for(const auto& s:spacers)shortest=std::min(shortest,s.size());
    std::string result;
    for(size_t n=1;n<shortest;++n) {
        std::map<std::string,size_t> counts;
        for(const auto& s:spacers)++counts[right?s.substr(s.size()-n):s.substr(0,n)];
        auto best=std::max_element(counts.begin(),counts.end(),[](const auto& a,const auto& b){return a.second<b.second;});
        if(double(best->second)+1e-12<options.support*spacers.size())break;
        result=best->first;
    }
    return result.size()>=options.min_bases?result:std::string{};
}
inline std::vector<BoundaryResult> refine_spacer_boundaries(const std::vector<BoundaryInput>& input,size_t k,const BoundaryOptions& options) {
    options.validate();
    std::map<std::string,std::set<std::string>> groups;
    for(const auto& s:input)groups[s.group].insert(s.spacer);
    std::map<std::string,std::pair<std::string,std::string>> motifs;
    if(options.enabled)for(const auto& [group,spacers]:groups)
        motifs[group]={conserved_terminal(spacers,false,options),conserved_terminal(spacers,true,options)};
    std::vector<std::pair<size_t,size_t>> proposals;
    std::map<std::string,std::set<std::pair<size_t,size_t>>> bounds;
    for(const auto& s:input) {
        const auto& m=motifs[s.group];size_t begin=0,end=s.spacer.size();
        if(!m.first.empty() && s.spacer.compare(0,m.first.size(),m.first)==0)begin=m.first.size();
        if(!m.second.empty() && s.spacer.compare(s.spacer.size()-m.second.size(),m.second.size(),m.second)==0)end-=m.second.size();
        proposals.push_back({begin,end});
        auto key=std::min(s.spacer,boundary_rc(s.spacer));
        auto oriented=key==s.spacer?std::make_pair(begin,end):std::make_pair(s.spacer.size()-end,s.spacer.size()-begin);
        bounds[key].insert(oriented);
    }
    std::vector<BoundaryResult> result;
    for(size_t i=0;i<input.size();++i) {
        const auto& s=input[i];auto [begin,end]=proposals[i];BoundaryResult r{s.repeat,s.spacer};
        auto key=std::min(s.spacer,boundary_rc(s.spacer));
        if(!options.enabled)r.status="disabled";
        else if(bounds[key].size()>1)r.status="conflicting_source_boundaries";
        else if(end<begin || end-begin<k)r.status="corrected_sequence_below_graph_k";
        else if(begin || end<s.spacer.size()) {
            r.left=begin;r.right=s.spacer.size()-end;
            r.spacer=s.spacer.substr(begin,end-begin);
            // R + A X B becomes (B R A) + X, a rotation of the same cycle.
            r.repeat=s.spacer.substr(end)+s.repeat+s.spacer.substr(0,begin);
            r.status="conserved_terminal_reassigned";
        }
        result.push_back(std::move(r));
    }
    return result;
}
} // namespace mcaat::proto
#endif
