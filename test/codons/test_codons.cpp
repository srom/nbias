#include "../src/codons/swap_codons.cpp"
#define BOOST_TEST_MODULE "Codons module test"
#include <boost/test/unit_test.hpp>

using namespace std;

BOOST_AUTO_TEST_CASE(SwapSynonymousCodonsTest)
{
    string sequence = "TTGTGCATG";
    
    // Basic test
    string swapped = SwapSynonymousCodons(sequence, 444);
    BOOST_CHECK_EQUAL(swapped, "CTTTGCATG");

    // Different seed = different outcome
    string swapped2 = SwapSynonymousCodons(sequence, 123);
    BOOST_CHECK_EQUAL(swapped2, "CTATGCATG");

    // Try with no seed
    string swapped3 = SwapSynonymousCodons(sequence);
    BOOST_CHECK_EQUAL(swapped3.size(), sequence.size());
    BOOST_CHECK_EQUAL(swapped3.substr(6), "ATG");

    // Test with short string
    string swapped4 = SwapSynonymousCodons("AA");
    BOOST_CHECK_EQUAL(swapped4, "AA");

    // Test when string length is not a multiple of 3
    string swapped5 = SwapSynonymousCodons("TTGTGCATGC", 444);
    BOOST_CHECK_EQUAL(swapped5, "CTTTGCATGC");

    // Test with stop codon
    string swapped6 = SwapSynonymousCodons("CTTTGCTAG", 666);
    BOOST_CHECK_EQUAL(swapped6, "CTGTGCTAA");
}

BOOST_AUTO_TEST_CASE(ShuffleCodonsTest)
{
    string sequence = "TTGTGCATGCCC";

    // Basic test
    string shuffled = ShuffleCodons(sequence, 444);
    BOOST_CHECK_EQUAL(shuffled, "CCCTGCTTGATG");

    // Different seed = different outcome
    string shuffled2 = ShuffleCodons(sequence, 123);
    BOOST_CHECK_EQUAL(shuffled2, "ATGTTGTGCCCC");

    // Try with no seed
    string shuffled3 = ShuffleCodons(sequence);
    BOOST_CHECK_EQUAL(shuffled3.size(), sequence.size());

    // Test with short string
    string shuffled4 = ShuffleCodons("AA");
    BOOST_CHECK_EQUAL(shuffled4, "AA");

    // Test when string length is not a multiple of 3
    BOOST_CHECK_EQUAL(ShuffleCodons("TTGTGCATGCCCA", 444), "CCCTGCTTGATGA");
    BOOST_CHECK_EQUAL(ShuffleCodons("TTGTGCATGCCCAT", 444), "CCCTGCTTGATGAT");
}
