#include "boost/test/unit_test.hpp"

#include "ref_context.hh"


BOOST_AUTO_TEST_SUITE( ref_context )


static
void
single_snp_hpol_test(const pos_t pos,
                     const reference_contig_segment& ref,
                     const unsigned expect) {

    const unsigned result(get_snp_hpol_size(pos,ref));
    BOOST_CHECK_EQUAL(result,expect);
}


BOOST_AUTO_TEST_CASE( test_snp_hpol_size ) {
    reference_contig_segment ref;
    ref.seq() = "TGTTTGAGATTT";

    single_snp_hpol_test(0,ref,2);
    single_snp_hpol_test(1,ref,5);
}



static
void
single_ihpol_test(const pos_t pos,
                  const reference_contig_segment& ref,
                  const unsigned expect) {

    const unsigned result(get_interupted_hpol_size(pos,ref));
    BOOST_CHECK_EQUAL(result,expect);
}


BOOST_AUTO_TEST_CASE( test_interupted_hpol_size ) {
    reference_contig_segment ref;
    ref.seq() = "TGTTTGAGATTT";

    single_ihpol_test(0,ref,4);
    single_ihpol_test(1,ref,4);
}


BOOST_AUTO_TEST_SUITE_END()

