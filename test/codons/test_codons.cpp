#include <random>
#include "../src/codons/swap_codons.cpp"
#define BOOST_TEST_MODULE "Codons module test"
#include <boost/test/unit_test.hpp>

using namespace std;

BOOST_AUTO_TEST_CASE(SwapSynonymousCodonsTest)
{
    string sequence = "TTGTGCATG";
    
    // Basic test
    string swapped = SwapSynonymousCodons(sequence, mt19937(444));
    BOOST_CHECK_EQUAL(swapped, "CTTTGCATG");

    // Different seed = different outcome
    string swapped2 = SwapSynonymousCodons(sequence, mt19937(123));
    BOOST_CHECK_EQUAL(swapped2, "CTATGCATG");

    // Test with short string
    string swapped4 = SwapSynonymousCodons("AA", mt19937(444));
    BOOST_CHECK_EQUAL(swapped4, "AA");

    // Test when string length is not a multiple of 3
    string swapped5 = SwapSynonymousCodons("TTGTGCATGC", mt19937(444));
    BOOST_CHECK_EQUAL(swapped5, "CTTTGCATGC");

    // Test with stop codon
    string swapped6 = SwapSynonymousCodons("CTTTGCTAG", mt19937(666));
    BOOST_CHECK_EQUAL(swapped6, "CTGTGCTAA");

    // Test with unknown codon
    string swapped7 = SwapSynonymousCodons("CTTXTCTAG", mt19937(666));
    BOOST_CHECK_EQUAL(swapped7, "CTGXTCTAA");
}

BOOST_AUTO_TEST_CASE(ShuffleCodonsTest)
{
    string sequence = "TTGTGCATGCCC";

    // Basic test
    string shuffled = ShuffleCodons(sequence, mt19937(444));
    BOOST_CHECK_EQUAL(shuffled, "CCCTGCTTGATG");

    // Different seed = different outcome
    string shuffled2 = ShuffleCodons(sequence, mt19937(123));
    BOOST_CHECK_EQUAL(shuffled2, "ATGTTGTGCCCC");

    // Test with short string
    string shuffled4 = ShuffleCodons("AA", mt19937(444));
    BOOST_CHECK_EQUAL(shuffled4, "AA");

    // Test when string length is not a multiple of 3
    BOOST_CHECK_EQUAL(ShuffleCodons("TTGTGCATGCCCA", mt19937(444)), "CCCTGCTTGATGA");
    BOOST_CHECK_EQUAL(ShuffleCodons("TTGTGCATGCCCAT", mt19937(444)), "CCCTGCTTGATGAT");
}
