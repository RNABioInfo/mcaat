#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>
using Row=std::map<std::string,std::string>;
std::vector<std::string> split(std::string s){if(!s.empty()&&s.back()=='\r')s.pop_back();std::vector<std::string> r;size_t a=0;while(true){auto b=s.find('\t',a);r.push_back(s.substr(a,b-a));if(b==std::string::npos)return r;a=b+1;}}
std::vector<Row> table(const std::string& path){std::ifstream f(path);if(!f)throw std::runtime_error("cannot open "+path);std::string line;std::getline(f,line);auto keys=split(line);std::vector<Row> rows;while(std::getline(f,line)){auto v=split(line);if(v.size()!=keys.size())throw std::runtime_error("bad table");Row r;for(size_t i=0;i<keys.size();++i)r[keys[i]]=v[i];rows.push_back(r);}return rows;}
int main(int argc,char** argv){try{
 if(argc!=4)throw std::runtime_error("usage: score_flanks CASES SELECTED OUTPUT");auto cases=table(argv[1]);std::map<std::string,Row> selected;for(auto r:table(argv[2]))selected[r["seed_id"]]=r;
 std::ofstream out(argv[3]);out<<"case_id\tseed_id\tside\tcheckpoint_ref_bp\tselected_flank_bp\tbest_prefix_query_bp\tweighted_cost\tqualifies_90pct\n";
 for(const auto& c:cases)for(bool back:{true,false}){
  auto ref=c.at(back?"left_outward":"right_outward");auto n=ref.find_first_not_of("ACGT");if(n!=std::string::npos)ref.resize(n);
  auto s=selected.at(c.at("seed_id"));auto seq=s.at("sequence");if(!seq.empty())seq=back?seq.substr(0,std::stoull(s.at("seed_begin"))):seq.substr(std::stoull(s.at("seed_end")));if(back)std::reverse(seq.begin(),seq.end());
  std::set<size_t> targets;for(size_t t:{24,100,250,500,750})if(t<=ref.size())targets.insert(t);if(!ref.empty())targets.insert(ref.size());
  std::vector<size_t> row(seq.size()+1),next(row.size());for(size_t j=0;j<row.size();++j)row[j]=9*j;
  for(size_t i=1;i<=ref.size();++i){next[0]=10*i;for(size_t j=1;j<row.size();++j)next[j]=std::min({row[j]+10,next[j-1]+9,row[j-1]+(ref[i-1]==seq[j-1]?0:10)});row.swap(next);
   if(targets.count(i)){auto it=std::min_element(row.begin(),row.end());out<<c.at("case_id")<<'\t'<<c.at("seed_id")<<'\t'<<(back?"left":"right")<<'\t'<<i<<'\t'<<seq.size()<<'\t'<<it-row.begin()<<'\t'<<*it<<'\t'<<(*it<=i)<<'\n';}
  }
 }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
