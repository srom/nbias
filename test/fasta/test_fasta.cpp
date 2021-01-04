#include <fstream>
#include "../src/fasta/fasta.cpp"
#define BOOST_TEST_MODULE "Fasta test"
#include <boost/test/unit_test.hpp>

using namespace std;


BOOST_AUTO_TEST_CASE(FastaGetAll)
{
	string path = "../fixtures/GCA_000010525.1_cds_from_genomic.fna";
	ifstream f(path, ios_base::in);
	BOOST_REQUIRE(f.good());

	FastaParser parser(f);
    auto records = parser.GetAll();

    BOOST_CHECK_EQUAL(records.size(), 4'717);

    auto record = records[0];

    BOOST_CHECK_EQUAL(record.id, "BAF85999.1"s);
    BOOST_CHECK_EQUAL(record.description, "GCA_000010525.1;ASM1052v1;ncbi;AP009384.1;351;1268;+"s);
    BOOST_CHECK_EQUAL(record.content.substr(0, 10), "TTGTGCTGCG"s);

    auto lastRecord = records[records.size() - 1];

    BOOST_CHECK_EQUAL(lastRecord.id, "BAF90715.1"s);
    BOOST_CHECK_EQUAL(lastRecord.description, "GCA_000010525.1;ASM1052v1;ncbi;AP009384.1;5368807;5369772;-"s);
    BOOST_CHECK_EQUAL(lastRecord.content.substr(0, 10), "ATGCGTCAGG"s);

    auto contentSize = lastRecord.content.size();
    BOOST_CHECK_EQUAL(lastRecord.content.substr(contentSize - 10, contentSize), "AACCCTCTGA"s);

    auto recordsEmpty = parser.GetAll();
    BOOST_CHECK_EQUAL(recordsEmpty.empty(), true);
}

BOOST_AUTO_TEST_CASE(FastaGet)
{
	string path = "../fixtures/GCA_000010525.1_cds_from_genomic.fna";
	ifstream f(path, ios_base::in);
	BOOST_REQUIRE(f.good());

	FastaParser parser(f);
	
	FastaRecord record1{};
	BOOST_CHECK_EQUAL(parser.Get(record1), true);
    BOOST_CHECK_EQUAL(record1.id, "BAF85999.1"s);
    BOOST_CHECK_EQUAL(record1.description, "GCA_000010525.1;ASM1052v1;ncbi;AP009384.1;351;1268;+"s);
    BOOST_CHECK_EQUAL(record1.content.substr(0, 10), "TTGTGCTGCG"s);

    FastaRecord record2{};
    BOOST_CHECK_EQUAL(parser.Get(record2), true);
    BOOST_CHECK_EQUAL(record2.id, "BAF86000.1"s);
    BOOST_CHECK_EQUAL(record2.description, "GCA_000010525.1;ASM1052v1;ncbi;AP009384.1;1265;1876;+"s);
    BOOST_CHECK_EQUAL(record2.content.substr(0, 10), "ATGAGCCTGT"s);
}
