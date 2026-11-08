#include "crispr_array/string_structure_filter.h"
#include <iostream>
#include <random>
#include <set>
using namespace mcaat::arrays;
size_t checks=0;
void check(bool ok,const std::string& message) {
    ++checks;if(!ok)throw std::runtime_error(message);
}
template<class F> void invalid(F f) {
    bool failed=false;try {f();}catch(const std::invalid_argument&) {failed=true;}
    check(failed,"invalid string-filter input must fail");
}
bool near(double a,double b) {return std::abs(a-b)<1e-10;}
// Independent exhaustive-template oracle, deliberately limited to short DNA.
double exhaustive_identity(const std::string& s) {
    size_t best=0;
    for(size_t p=1;p<=s.size()/3;++p)
        for(size_t code=0;code<(size_t(1)<<(2*p));++code) {
            size_t matches=0;
            for(size_t i=0;i<s.size();++i)matches+=s[i]=="ACGT"[(code>>(2*(i%p)))&3];
            best=std::max(best,matches);
        }
    return double(best)/s.size();
}
Quality separate_families(const std::vector<std::string>& spacers) {
    Quality q;for(const auto& s:std::set<std::string>(spacers.begin(),spacers.end()))q.families.push_back({s,{s}});
    return q;
}
StringEvidence score(const std::string& repeat,const std::vector<std::string>& spacers,
                     StringFilterOptions options={},size_t limit=500) {
    std::vector<std::pair<std::string,std::string>> rows;
    for(const auto& s:spacers)rows.emplace_back(repeat,s);
    return assess_strings(rows,separate_families(spacers),options,limit);
}
int main() {
    try {
        StringFilterOptions options;
        check(options.enabled && options.column_support==.85,"native filters enabled by default");
        check(set_string_filter_option("--string-filter","off",options) && !options.enabled,"explicit off");
        check(set_string_filter_option("string-filter","on",options) && options.enabled,"settings on");
        check(set_string_filter_option("--string-column-support",".6",options) && options.column_support==.6,"custom column support");
        check(!set_string_filter_option("unrelated","1",options),"unrelated parser option ignored");
        for(auto value:{"","0","-1","1.01","nan","inf","85%",".85junk"})
            invalid([&]{set_string_filter_option("string-motif-support",value,options);});
        invalid([&]{set_string_filter_option("string-filter","maybe",options);});
        check(set_string_filter_option("string-max-contrast","0",options),"zero contrast allowed");
        options.column_support=std::numeric_limits<double>::quiet_NaN();invalid([&]{options.validate();});
        options={};
        std::mt19937 rng(20261108);
        for(size_t n=3;n<=12;++n)for(int trial=0;trial<10;++trial) {
            std::string s;for(size_t i=0;i<n;++i)s+="ACGT"[rng()%4];
            check(near(short_motif_identity(s),exhaustive_identity(s)),"motif score agrees with exhaustive templates");
        }
        check(short_motif_identity("ACGTACGTACGTACGT")==1,"exact periodic unit");
        check(short_motif_identity("TCGTACGTACGTACGT")==15./16,"one substitution retained in motif score");
        for(auto rows:{std::vector<std::string>{},std::vector<std::string>{""},std::vector<std::string>{"ACTN"}})
            invalid([&]{observed_alignment(rows);});
        std::vector<std::string> originals={"ATGCGTACCGTA","ATGCGTACGTA","ATGCGTACCGTA","ATGCGTTACCGTA"};
        auto aligned=observed_alignment(originals);
        for(size_t i=0;i<originals.size();++i) {
            auto s=aligned[i];s.erase(std::remove(s.begin(),s.end(),'-'),s.end());
            check(s==originals[i] && aligned[i].size()==aligned[0].size(),"alignment preserves every observation");
        }
        check(aligned[0]==aligned[2],"duplicate strings have identical alignment");
        auto reversed=originals;std::reverse(reversed.begin(),reversed.end());
        auto reversed_alignment=observed_alignment(reversed);std::reverse(reversed_alignment.begin(),reversed_alignment.end());
        check(aligned==reversed_alignment,"alignment independent of input order");
        auto profile=base_conservation({"A-","--","--","--"},{.25,.25,.25,.25});
        check(profile==std::vector<double>({.25,0}),"gaps cannot supply conserved-base evidence");
        invalid([]{base_conservation({"AA","A"},{1,1});});
        invalid([]{base_conservation({"A"},{0});});
        invalid([]{base_conservation({"N"},{1});});
        const std::string repeat="AGTCGATGCATCAGGATCCGATGA";
        const std::vector<std::string> spacers={"ATGCAACGTTGC","CGATGGCATACG","GTCATTAGCGTA","TACGCTTACGAT"};
        auto evidence=score(repeat,spacers);
        check(evidence.evaluated && evidence.alignment_scored && evidence.reason().empty(),"variable spacers remain accepted");
        check(near(evidence.repeat_conservation,1),"identical extracted repeats score one");
        auto permuted=spacers;std::reverse(permuted.begin(),permuted.end());
        auto other=score(repeat,permuted);
        check(near(evidence.spacer_conservation,other.spacer_conservation) && evidence.reason()==other.reason(),"decision invariant under row order");
        auto duplicates=spacers;for(int i=0;i<20;++i)duplicates.push_back(spacers.front());
        other=score(repeat,duplicates);
        check(near(evidence.spacer_conservation,other.spacer_conservation) && near(evidence.motif_fraction,other.motif_fraction),"duplicate observations cannot overweight a family");
        evidence=score(repeat,{repeat+"A",repeat+"C",repeat+"G"});
        check(evidence.conserved_spacer && evidence.reason()=="conservation_across_spacer","conserved spacer columns trigger rule");
        evidence=score("ACGTACGTACGT",{"ACGTACGTACGT","ACGTACGTACGA"});
        check(evidence.short_motif,"short motifs trigger rule");
        evidence=score("ACGTACGTACGT",{"ACGTACGTACGT"});
        check(evidence.reason().empty(),"one similarity family alone is not rejected");
        evidence=score("ACGTACGTACGT",{"ACGTACGTACGT","ACGTACGTACGA"},{},2);
        check(evidence.alignment_skipped && !evidence.alignment_scored && evidence.short_motif,"bounded alignment work still permits motif scoring");
        options.enabled=false;evidence=score("ACGTACGTACGT",{"ACGTACGTACGT","ACGTACGTACGA"},options);
        check(!evidence.evaluated && evidence.reason().empty(),"off disables both rules");
        std::cout<<checks<<" native string-filter checks passed\n";
    }catch(const std::exception& e){std::cerr<<"FAILED: "<<e.what()<<'\n';return 1;}
}
