#include "../src/count/tri_count.cpp"
#define BOOST_TEST_MODULE "Tri count test"
#include <boost/test/unit_test.hpp>

using namespace std;

BOOST_AUTO_TEST_CASE(TriCountTest)
{
    string sequence = "TTGTGCTGCAAAAGCTTGCATGGGGCCGGAGGGCATGCCTTCCTGCACACGCCGTCCACAGACAAA";

    BOOST_CHECK_EQUAL(CountKmer(sequence, "AAA", true), 3);
    BOOST_CHECK_EQUAL(CountKmer(sequence, "TTG", true), 2);
    BOOST_CHECK_EQUAL(CountKmer(sequence, "GCA", true), 4);

    BOOST_CHECK_EQUAL(CountKmer(sequence, "AAA", false), 2);
    BOOST_CHECK_EQUAL(CountKmer(sequence, "TTG", false), 2);
    BOOST_CHECK_EQUAL(CountKmer(sequence, "GCA", false), 4);

    auto counts = CountTriNucleotides(sequence, true, true);

    BOOST_CHECK_EQUAL(counts.size(), 64);
    BOOST_CHECK_EQUAL(counts[0], 3);   // AAA
    BOOST_CHECK_EQUAL(counts[63], 0);  // TTT

    auto counts_no_overlap = CountTriNucleotides(sequence, false, true);

    BOOST_CHECK_EQUAL(counts_no_overlap.size(), 64);
    BOOST_CHECK_EQUAL(counts_no_overlap[0], 2);   // AAA
    BOOST_CHECK_EQUAL(counts_no_overlap[63], 0);  // TTT

    auto counts_sequential = CountTriNucleotides(sequence, true, false);

    BOOST_CHECK_EQUAL(counts_sequential.size(), 64);
    BOOST_CHECK_EQUAL(counts_sequential[0], 3);   // AAA
    BOOST_CHECK_EQUAL(counts_sequential[63], 0);  // TTT
}
