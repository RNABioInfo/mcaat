#include <gtest/gtest.h>
#include "reads.h"

#include <string>
#include <unordered_map>

using namespace std;


// Normal test
TEST(ReadsTest, ReversePairEndsSequenceBasic) {
    const string input = "ACGT";
    const string expected = "ACGT";

    string actual = input;
    reverse_pair_ends_sequence(actual);

    EXPECT_EQ(actual, expected);
}

// Normal test
TEST(ReadsTest, ReversePairEndsSequenceBasic2) {
    string input = "ACGTT";
    const string expected = "AACGT";

    string actual = input;
    reverse_pair_ends_sequence(actual);

    EXPECT_EQ(actual, expected);
}

// Normal test
TEST(ReadsTest, ReversePairEndsSequenceBasic3) {
    string input = "AAGCT";
    const string expected = "AGCTT";

    string actual = input;
    reverse_pair_ends_sequence(actual);

    EXPECT_EQ(actual, expected);
}

// Special test 'X' is in string
TEST(ReadsTest, ReversePairEndsSequenceNonStandardChars) {
    string input = "AXGT";
    const string expected = "ACXT";

    string actual = input;
    reverse_pair_ends_sequence(actual);

    EXPECT_EQ(actual, expected);
}

// Special test: Many non 'ACGT' chracters are in the string
TEST(ReadsTest, ReversePairEndsSequenceNonStandardChars2) {
    const string input = "XXXQUIOCPOPYM";
    const string expected = "MYPOPGOIUQXXX";

    string actual = input;
    reverse_pair_ends_sequence(actual);

    EXPECT_EQ(actual, expected);
}
