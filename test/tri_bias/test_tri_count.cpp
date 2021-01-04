#include "../src/tri_bias/tri_count.cpp"
#define BOOST_TEST_MODULE "Tri count test"
#include <boost/test/unit_test.hpp>

using namespace std;

BOOST_AUTO_TEST_CASE(TriCountTest)
{
    string sequence = "TTGTGCTGCAAAAGCTTGCATGGGGCCGGAGGGCATGCCTTCCTGCACACGCCGTCCACAGACCAAA";
    BOOST_CHECK_EQUAL(CountTriNucleotide(sequence, "AAA"), 3);
    BOOST_CHECK_EQUAL(CountTriNucleotide(sequence, "TTG"), 2);
    BOOST_CHECK_EQUAL(CountTriNucleotide(sequence, "GCA"), 4);

    auto counts = CountTriNucleotides(sequence);

    BOOST_CHECK_EQUAL(counts.size(), 64);
    BOOST_CHECK_EQUAL(counts[0], 3);   // AAA
    BOOST_CHECK_EQUAL(counts[63], 0);  // TTT
}
