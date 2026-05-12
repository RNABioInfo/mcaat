#ifndef EVALUATION_H_
#define EVALUATION_H_

#include <unordered_map>
#include <vector>
#include <string>
#include <stdint.h>

using namespace std;

/**
 * @internal
 * @brief Get the levenshtein distance of both strings
 * 
 * Using the Levenshtein distance (matrix implementation in O(|s1| |s2|)) with:
 *  - insert(c)   costs 1
 *  - update(a,b) costs 1
 *  - delete(c)   costs 1
 * 
 * @param s1 First string
 * @param s2 Second string
 * @return Distance according to the algorithm (uint16_t)
 */
uint16_t get_levenshtein_distance(const string& s1, const string& s2);

/**
 * @internal
 * @brief Compares strings using levenshtein distance
 * 
 * With |s1| = n, |s2| = m, levenshtein distance d
 * The result is 1 - (d / max(n, m))
 * 
 * @param s1 First sequence
 * @param s2 Second sequence
 * 
 * @return Value between 0 and 1 for how similar the sequences are (float)
 */
float get_string_similarity(const string& s1, const string& s2);

/**
 * @brief Get the number of duplicate spacers
 * 
 * Simply counts the amount of occurences of every spacer.
 * If it occurs more than once, it adds it to the result
 * 
 * @param spacers 
 * @param expected_sequence 
 * 
 * @return Amount of spacer duplicates (int)
 */
int get_number_of_duplicate_spacers(
    const vector<string>& spacers,
    const string& expected_sequence
);

/**
 * @internal
 * @brief Chooses the most similar sequence and removes it from the vector
 * 
 * @post choices will be modified and one elements will be removed (if possible)
 * 
 * @param sequence Used to compare and find the most similar
 * @param choices Where to pop the sequence from
 * 
 * @return The most similar sequence (string)
 * @return "", if choices is empty
 */
string get_most_similar_sequence(
    const string& sequence,
    const vector<string>& choices
);

#endif
