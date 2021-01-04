#include "../src/dna/dna.cpp"
#define BOOST_TEST_MODULE "DNA test"
#include <boost/test/unit_test.hpp>

using namespace std;

BOOST_AUTO_TEST_CASE(TestReverseComplement)
{
    string sequence = "TTGTGCTGCAAAAG";
    ReverseComplement(sequence);
    BOOST_CHECK_EQUAL(sequence, "CTTTTGCAGCACAA"s);
}
