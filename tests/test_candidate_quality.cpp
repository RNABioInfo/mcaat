#include "crispr_array/candidate_quality.h"
#include <iostream>
#include <random>
using namespace mcaat::arrays;
size_t checks=0;
void check(bool result,const char* message) {++checks;if(!result)throw std::runtime_error(message);}
size_t reference(const std::string& a,const std::string& b) {
    std::vector<std::vector<size_t>> d(a.size()+1,std::vector<size_t>(b.size()+1));
    for(size_t i=0;i<=a.size();++i)d[i][0]=i;
    for(size_t j=0;j<=b.size();++j)d[0][j]=j;
    for(size_t i=1;i<=a.size();++i)for(size_t j=1;j<=b.size();++j)
        d[i][j]=std::min({d[i-1][j]+1,d[i][j-1]+1,d[i-1][j-1]+size_t(a[i-1]!=b[j-1])});
    return d.back().back();
}
int main() {
    try {
        std::mt19937 random(41);
        for(size_t i=0;i<2500;++i) {
            std::string a,b;
            for(size_t j=random()%101;j--;)a+="ACGT"[random()%4];
            for(size_t j=random()%101;j--;)b+="ACGT"[random()%4];
            check(edit_distance(a,b)==reference(a,b),"bit-vector distance differs from independent DP");
        }
        QualityOptions options;
        std::string a(30,'A'),b=a,c=a;b[0]=b[1]='C';c[0]=c[1]=c[2]=c[3]='C';
        auto bridge=assess_quality(std::string(30,'G'),{a,b,c},3,options);
        check(bridge.families.size()==2,"similarity bridge incorrectly merged dissimilar endpoints");
        size_t variants=0;for(auto& f:bridge.families)variants+=f.variants.size();
        check(variants==3,"variant grouping lost an observed sequence");
        auto low=assess_quality(std::string(30,'G'),{a,b},2,options);
        check(low.reason=="fewer_than_two_spacer_families","one spacer plus variants accepted as independent spacers");
        auto tandem=assess_quality(a,{std::string(9,'C')+std::string(21,'A'),std::string(21,'A')+std::string(9,'G')},2,options);
        check(tandem.reason=="repeat_like_spacers","repeat-like spacers accepted");
        auto diverse=assess_quality(std::string(30,'G'),{a,std::string(30,'T')},2,options);
        check(diverse.reason.empty() && diverse.families.size()==2,"dissimilar spacer families rejected");
        check(assess_quality(a,{a,b},500,options).reason=="cycle_limit","merged cycle bound ignored");
        options.variant_similarity=std::numeric_limits<double>::quiet_NaN();
        bool rejected=false;try{options.validate();}catch(const std::invalid_argument&){rejected=true;}
        check(rejected,"nonfinite quality threshold accepted");
        std::cout<<checks<<" candidate-quality checks passed\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
