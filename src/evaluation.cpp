#include "evaluation.h"

uint16_t get_levenshtein_distance(const string& s1, const string& s2) {
    vector<vector<uint16_t>> dist;

    /* Setup matrix dist */
    vector<uint16_t> first_row;
    for (int i = 0; i < s1.size() + 1; ++i) {
        first_row.push_back(i);
    }
    dist.push_back(first_row);

    for (int i = 1; i < s2.size() + 1; ++i) {
        vector<uint16_t> row;
        row.push_back(i);
        dist.push_back(row);
    }

    /* Fill matrix with computation */
    for (int x = 1; x < s1.size() + 1; ++x) {
        for (int y = 1; y < s2.size() + 1; ++y) {
            char s1_char = s1.at(x - 1);
            char s2_char = s2.at(y - 1);

            uint16_t min;
            if (s1_char == s2_char) { // No extra distance
                min = dist.at(y - 1).at(x - 1);
            } else { // update(s1_char, s2_char)
                min = dist.at(y - 1).at(x - 1) + 1;
            }

            uint16_t insert_cost = dist.at(y).at(x - 1) + 1;
            if (min > insert_cost) { // insert(s1_char)
                min = insert_cost;
            }

            uint16_t delete_cost = dist.at(y - 1).at(x) + 1;
            if (min > delete_cost) { // delete()
                min = delete_cost;
            }

            dist.at(y).push_back(min);
        }
    }

    /* Return result of matrix */
    return dist.at(s2.size()).at(s1.size());
}

float get_string_similarity(const string& s1, const string& s2) {
    uint16_t d = get_levenshtein_distance(s1, s2);
    size_t max_size = s1.size() > s2.size() ? s1.size() : s2.size();
    
    return 1.0 - (static_cast<float>(d) / static_cast<float>(max_size));
}

int get_number_of_duplicate_spacers(
    const vector<string>& spacers,
    const string& expected_sequence
) {
    int result = 0;

    for (const auto& spacer : spacers) {
        int spacer_count = 0;

        size_t pos = 0;
        while ((pos = expected_sequence.find(spacer, pos)) != string::npos) {
            ++spacer_count;
            ++pos;
        }

        if (spacer_count > 1) {
            result += spacer_count - 1;
        }
    }

    return result;
}

string get_most_similar_sequence(
    const string& sequence,
    const vector<string>& choices
) {
    if (choices.empty()) {
        return "";
    }

    float best_similarity = -1.0;
    int best_idx = -1;
    string result_sequence = "";
    for (int i = 0; i < choices.size(); ++i) {
        const string choice_sequence = choices.at(i);
        float current_similarity = get_string_similarity(sequence, choice_sequence);

        if (current_similarity > best_similarity) {
            best_similarity = current_similarity;
            best_idx = i;
            result_sequence = choice_sequence;
        }
    }

    // Remove the choice from the compares:
    // choices.erase(choices.begin() + best_idx);

    return result_sequence;
}
