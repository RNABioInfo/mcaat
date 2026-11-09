#ifndef MCAAT_PROTO_EXTENSION_EVIDENCE_H
#define MCAAT_PROTO_EXTENSION_EVIDENCE_H
#include "protospacer_detection/extension_engine.h"
#include <functional>
#include <cctype>
#include <sstream>
#include <unordered_map>

namespace mcaat::proto {
inline bool is_dna(const std::string& s) {return !s.empty() && s.find_first_not_of("ACGT")==std::string::npos;}
inline std::string upper(std::string s) {for(char& c:s)c=std::toupper(static_cast<unsigned char>(c));return s;}
inline bool iupac(char motif,char base) {
    std::string allowed;
    switch(motif) {
        case 'A':case 'C':case 'G':case 'T':allowed=motif;break;
        case 'R':allowed="AG";break;case 'Y':allowed="CT";break;case 'S':allowed="GC";break;
        case 'W':allowed="AT";break;case 'K':allowed="GT";break;case 'M':allowed="AC";break;
        case 'B':allowed="CGT";break;case 'D':allowed="AGT";break;case 'H':allowed="ACT";break;
        case 'V':allowed="ACG";break;case 'N':allowed="ACGT";break;
        default:throw std::invalid_argument("invalid IUPAC PAM motif");
    }
    return allowed.find(base)!=std::string::npos;
}
struct PamRule {std::string array,subtype,strand,side,motif;};
inline std::string check_pam(const std::string& sequence,size_t begin,size_t end,
                             bool query_reverse,const PamRule& rule,bool full_match) {
    if(rule.side=="none")return "not_applicable";
    if(!full_match)return "unknown_partial_match";
    if(rule.strand=="?")return "unknown_orientation";
    bool flip=query_reverse^(rule.strand=="-");
    std::string s=flip?reverse_complement(sequence):sequence;
    if(flip){size_t b=s.size()-end;end=s.size()-begin;begin=b;}
    if((rule.side=="5prime" && begin<rule.motif.size()) ||
       (rule.side=="3prime" && s.size()-end<rule.motif.size()))return "unknown_short_flank";
    size_t position=rule.side=="5prime"?begin-rule.motif.size():end;
    for(size_t i=0;i<rule.motif.size();++i)if(!iupac(rule.motif[i],s[position+i]))return "mismatch";
    return "match";
}
struct ReadRegion {std::string sequence;std::vector<std::string> junctions;size_t pairs=0,unique_pairs=0;};
struct ReadMatch {size_t region,position,length;bool reverse;};
class ReadEvidence {
    size_t k_,max_entries_,max_work_;
    std::unordered_map<std::string,std::vector<std::pair<size_t,size_t>>> anchors_;
    std::unordered_map<std::string,size_t> junction_counts_;
public:
    std::vector<ReadRegion> regions;
    size_t reads=0,pairs=0;
    ReadEvidence(size_t k,size_t max_entries,size_t max_work=1000000):k_(k),max_entries_(max_entries),max_work_(max_work) {}
    void charge(size_t& work) const {
        if(++work>max_work_)throw std::runtime_error("ambiguous read pair exceeds --max-read-work; increase it or reduce retained alternatives");
    }
    void index(std::vector<ReadRegion> input) {
        regions=std::move(input);size_t entries=0;
        for(size_t id=0;id<regions.size();++id) {
            const auto& r=regions[id];
            for(size_t pos=0;pos+k_<=r.sequence.size();++pos) {
                if(++entries>max_entries_)throw std::runtime_error("read-support index exceeds --max-read-index; increase it or reduce candidate/search limits");
                anchors_[r.sequence.substr(pos,k_)].push_back({id,pos});
            }
            for(const auto& word:r.junctions)junction_counts_.emplace(canonical(word),0);
        }
    }
    void read(const std::string& sequence) {
        ++reads;std::unordered_set<std::string> seen;
        for(size_t i=0;i+k_+2<=sequence.size();++i) {
            auto word=sequence.substr(i,k_+2);if(!is_dna(word))continue;
            auto key=canonical(word);if(junction_counts_.count(key))seen.insert(std::move(key));
        }
        for(const auto& key:seen)++junction_counts_.at(key);
    }
    std::vector<ReadMatch> matches(const std::string& read,size_t& work) const {
        std::vector<ReadMatch> result;if(read.size()<k_ || !is_dna(read))return result;
        for(bool reverse:{false,true}) {
            auto s=reverse?reverse_complement(read):read;
            auto it=anchors_.find(s.substr(0,k_));if(it==anchors_.end())continue;
            for(auto [id,pos]:it->second) {
                charge(work);
                if(regions[id].sequence.compare(pos,s.size(),s)==0)result.push_back({id,pos,s.size(),reverse});
            }
        }
        return result;
    }
    void pair(const std::string& r1,const std::string& r2,size_t min_insert,size_t max_insert) {
        ++pairs;size_t work=0;auto a=matches(r1,work),b=matches(r2,work);std::set<size_t> supported;
        // Join placements within each region, avoiding a Cartesian product over
        // unrelated candidates. Bound repetitive-read work explicitly as well.
        std::unordered_map<size_t,std::vector<const ReadMatch*>> mates;
        for(const auto& y:b)mates[y.region].push_back(&y);
        for(const auto& x:a) {
          auto it=mates.find(x.region);if(it==mates.end() || supported.count(x.region))continue;
          for(const auto* mate:it->second) {
            charge(work);const auto& y=*mate;
            if(x.reverse==y.reverse)continue;
            const auto& left=x.reverse?y:x;const auto& right=x.reverse?x:y;
            if(left.position>right.position)continue;
            size_t span=std::max(left.position+left.length,right.position+right.length)-left.position;
            if(span>=min_insert && span<=max_insert){supported.insert(x.region);break;}
          }
        }
        for(auto id:supported)++regions[id].pairs;
        if(supported.size()==1)++regions[*supported.begin()].unique_pairs;
    }
    size_t junction_count(const std::string& s) const {auto it=junction_counts_.find(canonical(s));return it==junction_counts_.end()?0:it->second;}
};
} // namespace mcaat::proto
#endif
