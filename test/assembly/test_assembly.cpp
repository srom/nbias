#include "../src/assembly/assembly.cpp"
#define BOOST_TEST_MODULE "Assembly test"
#include <boost/test/unit_test.hpp>

using namespace std;

BOOST_AUTO_TEST_CASE(TestAssemblyClass)
{
	string path = "../fixtures/assemblies.csv";

	Assemblies assemblies(path);

    BOOST_CHECK_EQUAL(assemblies.Size(), 10);
    BOOST_CHECK_EQUAL(assemblies.GetIds().size(), 10);

    BOOST_CHECK_EQUAL(assemblies.Has("GCA_001735525.1"s), true);

    auto assembly = assemblies.Get("GCA_001735525.1"s);

    BOOST_CHECK_EQUAL(assembly.organism_name, "Shewanella colwelliana"s);
    BOOST_CHECK_EQUAL(assembly.species_taxid, 23);

    BOOST_CHECK_EQUAL(assemblies.Has("does not exist"s), false);

    auto missingAssembly = assemblies.Get("does not exist"s);

    BOOST_CHECK_EQUAL(missingAssembly.assembly_accession.empty(), true);
}
