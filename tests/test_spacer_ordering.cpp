#include <gtest/gtest.h>
#include "spacer_ordering.h"

#include <unordered_set>
#include <vector>

using namespace std;


// Border case: Empty universe
TEST(SolveMinCoverProblemTest, EmptyUniverse) {
    unordered_set<uint32_t> universe;
    vector<vector<uint32_t>> sets = {{0, 1}, {2, 3}};
    
    auto result = solve_min_cover_problem(universe, sets);
    
    EXPECT_TRUE(result.empty());
}

// Border case: Empty sets
TEST(SolveMinCoverProblemTest, EmptySets) {
    unordered_set<uint32_t> universe = {0, 1, 2, 3};
    vector<vector<uint32_t>> sets;
    
    auto result = solve_min_cover_problem(universe, sets);
    
    EXPECT_TRUE(result.empty());
}

// Border case: No solution possible (sets don't cover universe)
TEST(SolveMinCoverProblemTest, NoSolutionPossible) {
    unordered_set<uint32_t> universe = {0, 1, 2};
    vector<vector<uint32_t>> sets = {{0, 1}, {3, 4}};
    
    auto result = solve_min_cover_problem(universe, sets);
    
    EXPECT_TRUE(result.empty());
}

// Special case: Single element universe, single set
TEST(SolveMinCoverProblemTest, SingleElementSingleSet) {
    unordered_set<uint32_t> universe = {0};
    vector<vector<uint32_t>> sets = {{0}};
    
    auto result = solve_min_cover_problem(universe, sets);
    
    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(result[0], 0);
}

// Normal case: Simple optimal solution
TEST(SolveMinCoverProblemTest, SimpleOptimalSolution) {
    unordered_set<uint32_t> universe = {0, 1, 2, 3, 4};
    vector<vector<uint32_t>> sets = {{0, 1, 2}, {3, 4}, {1, 3}, {2, 4}};
    
    auto result = solve_min_cover_problem(universe, sets);
    
    EXPECT_FALSE(result.empty());
    EXPECT_LE(result.size(), 2);
    
    unordered_set<uint32_t> covered;
    for (size_t idx : result) {
        for (uint32_t element : sets[idx]) {
            covered.insert(element);
        }
    }
    EXPECT_EQ(covered, universe);
}

// Normal case: Redundant sets
TEST(SolveMinCoverProblemTest, RedundantSets) {
    unordered_set<uint32_t> universe = {0, 1, 2, 3};
    vector<vector<uint32_t>> sets = {{0, 1, 2, 3}, {0, 1}, {2}, {3}, {1, 2}};
    
    auto result = solve_min_cover_problem(universe, sets);
    
    EXPECT_FALSE(result.empty());
    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(result[0], 0);
}

// Normal case: Complex overlapping sets
TEST(SolveMinCoverProblemTest, ComplexOverlappingSets) {
    unordered_set<uint32_t> universe = {0, 1, 2, 3, 4, 5, 6};
    vector<vector<uint32_t>> sets = {
        {0, 1, 2, 3},
        {0, 2, 4, 5},
        {0, 3, 5, 6},
        {0, 1, 4},
        {0, 6}
    };
    
    auto result = solve_min_cover_problem(universe, sets);
    
    EXPECT_FALSE(result.empty());
    
    unordered_set<uint32_t> covered;
    for (size_t idx : result) {
        EXPECT_LT(idx, sets.size());
        for (uint32_t element : sets[idx]) {
            covered.insert(element);
        }
    }
    EXPECT_EQ(covered, universe);
    
    EXPECT_LE(result.size(), 3);
}
