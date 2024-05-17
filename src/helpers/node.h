#ifndef node_h
#define node_h
#include <cstdint>

using namespace std;

/// @brief The helper node struct serves as a container for the node id and its multiplicity
/// @brief and is used to sort the nodes based on their multiplicity.
struct node
{
    /// @brief Node ID
    uint64_t id;

    /// @brief Node Multiplicity
    uint64_t multiplicity;

    /// @brief Constructor
    /// @param id Node ID
    /// @param multiplicity Node Multiplicity
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