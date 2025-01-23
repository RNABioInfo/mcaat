#ifndef NODE_H
#define NODE_H
#include <cstdint>
#include <functional>

struct node {
    // Add members for node
    uint64_t id; // Example member
    uint64_t multiplicity;

    bool operator==(const node& other) const {

        return id == other.id;
    }
};

// Specialization of std::hash for the node structure
namespace std {
    template<>
    struct hash<node> {
        std::size_t operator()(const node& n) const {
            return std::hash<uint64_t>()(n.id); // Use appropriate member(s) for hashing
        }
    };
}
#endif