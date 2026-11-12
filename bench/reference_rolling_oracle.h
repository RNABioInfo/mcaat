#ifndef MCAAT_REFERENCE_ROLLING_ORACLE_H
#define MCAAT_REFERENCE_ROLLING_ORACLE_H
#include "reference_graph_oracle.h"
#include <cstdio>
#include <memory>
namespace mcaat::proto::audit {
// The reference position never decreases. Keep just one alignment column in
// memory; close insertion edges with Dijkstra before advancing that column.
// Immutable traceback records go to an automatically removed temporary file.
// This is diagnostic alignment, never a production path-selection method.
template<class Graph> Search align_rolling(Graph& graph,Node start,const std::string& reference,bool back,
        const std::set<size_t>& targets,size_t max_states,
        const std::function<size_t(size_t,size_t)>& suffix_bound={},const Path& extra_starts={},
        bool first_only=false,PrefixContext* context=nullptr,CycleContext* cycles=nullptr) {
    Search out;if(reference.empty() || (start==missing&&extra_starts.empty()))return out;
    std::vector<int64_t> admissible(reference.size()+1,-1);
    for(size_t pos=0;pos<=reference.size();++pos)for(auto target:targets)if(target>=pos)
        admissible[pos]=std::max(admissible[pos],int64_t(target)-int64_t(suffix_bound?suffix_bound(pos,target):0));
    struct Record {Node node;size_t parent;char operation;};
    struct Cell {size_t cost,record;};
    std::unique_ptr<FILE,decltype(&std::fclose)> file(std::tmpfile(),&std::fclose);
    if(!file)throw std::runtime_error("cannot create temporary audit traceback");
    size_t records=0;
    auto save=[&](Node node,size_t parent,char operation) {
        Record record{node,parent,operation};
        if(std::fwrite(&record,sizeof(record),1,file.get())!=1)throw std::runtime_error("cannot write audit traceback");
        return records++;
    };
    using Column=std::unordered_map<Key,Cell,KeyHash>;
    Column current,next;
    auto insert=[&](Column& column,Key key,size_t cost,size_t parent,char op,size_t pos) {
        if(int64_t(cost)>admissible[pos])return false;
        auto it=column.find(key);if(it!=column.end()&&it->second.cost<=cost)return false;
        if(records>=max_states){out.limited=true;out.lowest_unexplored_cost=0;return false;}
        column[key]={cost,save(key.node,parent,op)};return true;
    };
    std::set<Node> roots(extra_starts.begin(),extra_starts.end());if(start!=missing)roots.insert(start);
    for(auto root:roots)insert(current,{root,0,context?context->initial():0,0},0,size_t(-1),'S',0);
    std::unordered_map<Node,std::vector<Arc>> adjacency;
    auto witness=[&](size_t record,size_t pos) {
        std::vector<Record> chain;
        while(record!=size_t(-1)) {
            Record value;
            if(fseeko(file.get(),off_t(record*sizeof(Record)),SEEK_SET)!=0 || std::fread(&value,sizeof(value),1,file.get())!=1)
                throw std::runtime_error("cannot read audit traceback");
            chain.push_back(value);record=value.parent;
        }
        if(fseeko(file.get(),0,SEEK_END)!=0)throw std::runtime_error("cannot resume audit traceback");
        std::reverse(chain.begin(),chain.end());Witness w;w.ref_bp=pos;w.nodes={chain.front().node};
        for(size_t i=1;i<chain.size();++i) {
            auto& step=chain[i];w.operations+=step.operation;++w.columns;
            if(step.operation=='M')++w.matches;else ++w.edits;
            if(step.operation!='D'){w.nodes.push_back(step.node);auto label=graph.label(step.node);w.sequence+=back?label.front():label.back();}
        }
        if(w.matches*100<90*w.columns)throw std::runtime_error("rolling oracle witness below identity threshold");
        auto path=w.nodes;if(back)std::reverse(path.begin(),path.end());auto spelled=graph.spell(path);
        auto emitted=back?spelled.substr(0,spelled.size()-graph.k()):spelled.substr(graph.k());if(back)std::reverse(emitted.begin(),emitted.end());
        if(emitted!=w.sequence)throw std::runtime_error("rolling oracle witness not spelled by graph");
        return w;
    };
    bool complete=false;
    for(size_t pos=0;pos<=reference.size()&&!out.limited&&!complete;++pos) {
        using Item=std::pair<size_t,size_t>; // cost, immutable traceback record
        std::priority_queue<Item,std::vector<Item>,std::greater<Item>> queue;
        std::unordered_map<size_t,Key> queued;
        for(const auto& [key,cell]:current){queue.emplace(cell.cost,cell.record);queued.emplace(cell.record,key);}
        while(!queue.empty()&&!out.limited) {
            auto [cost,record]=queue.top();queue.pop();auto item=queued.find(record);auto key=item->second;queued.erase(item);
            auto cell=current.find(key);if(cell->second.record!=record)continue;
            ++out.popped;
            if(!key.context && targets.count(pos)&&cost<=pos&&!out.found.count(pos)) {
                out.found.emplace(pos,witness(record,pos));
                if(first_only||out.found.size()==targets.size()){complete=true;break;}
            }
            if(pos==reference.size())continue;
            insert(next,key,cost+10,record,'D',pos+1);
            auto cached=adjacency.find(key.node);
            if(cached==adjacency.end()) {
                if(adjacency.size()>100000)adjacency.clear();
                cached=adjacency.emplace(key.node,graph.neighbors(key.node,back)).first;
            }
            for(const auto& arc:cached->second) {
                auto base=graph.extension_base(arc,back);
                auto prefix=context?context->advance(key.context,base):0;if(prefix==PrefixContext::reject)continue;
                auto cycle=cycles?cycles->advance(graph,key.cycle,key.node,arc,back):0;if(cycle==CycleContext::reject)continue;
                Key dest{arc.to,0,prefix,cycle};bool same=base==reference[pos];
                insert(next,dest,cost+(same?0:10),record,same?'M':'X',pos+1);
                if(insert(current,dest,cost+9,record,'I',pos)) {
                    auto updated=current.at(dest);queue.emplace(updated.cost,updated.record);queued.emplace(updated.record,dest);
                }
            }
        }
        current=std::move(next);next.clear();
        if(current.empty())break;
    }
    out.discovered=records;out.exhausted=!out.limited&&!complete;return out;
}
}
#endif
