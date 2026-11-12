#ifndef MCAAT_REFERENCE_PREFIX_CONTEXT_H
#define MCAAT_REFERENCE_PREFIX_CONTEXT_H
#include "protospacer_detection/extension_engine.h"
#include <array>
namespace mcaat::proto::audit {
// Lazy exact automaton for the production initial source-repeat rule.
// Pending strings are separate alignment states; after an external context is
// established they merge into state 0. This preserves path-dependent context
// during graph-node/reference-position dynamic programming.
class PrefixContext {
    struct Prefix {std::string text;std::array<size_t,4> next{{missing,missing,missing,missing}};};
    std::vector<Prefix> prefixes_{{""},{""}};
    std::vector<FlankContext> models_;
    double identity_;
public:
    static constexpr size_t reject=missing-1;
    PrefixContext(std::vector<FlankContext> models,double identity=.85):models_(std::move(models)),identity_(identity){}
    size_t initial()const{return context_state("",models_,identity_)==ContextState::external?0:1;}
    size_t advance(size_t state,char base) {
        if(!state)return 0;
        auto letter=std::string("ACGT").find(base);
        auto cached=prefixes_.at(state).next.at(letter);if(cached!=missing)return cached;
        auto text=prefixes_.at(state).text+base;
        auto status=context_state(text,models_,identity_);size_t next=0;
        if(status==ContextState::array)next=reject;
        else if(status==ContextState::pending){next=prefixes_.size();prefixes_.push_back({std::move(text)});}
        prefixes_.at(state).next.at(letter)=next;return next;
    }
    size_t size()const{return prefixes_.size();}
};
}
#endif
