#ifndef node_h
#define node_h
#include <cstdint>

using namespace std;


struct node
{
    uint64_t id;
    uint64_t multiplicity;

    node(uint64_t id, uint64_t multiplicity) : id(id), multiplicity(multiplicity) {}

    bool operator<(const node &rhs) const
    {
        return multiplicity < rhs.multiplicity;
    }

    bool operator>(const node &rhs) const
    {
        return multiplicity > rhs.multiplicity;
    }

    bool operator<=(const node &rhs) const
    {
        return multiplicity <= rhs.multiplicity;
    }

    bool operator>=(const node &rhs) const
    {
        return multiplicity >= rhs.multiplicity;
    }

    bool operator==(const node &rhs) const
    {
        return multiplicity == rhs.multiplicity;
    }

    bool operator!=(const node &rhs) const
    {
        return multiplicity != rhs.multiplicity;
    }

    
};

#endif