#include "boost/test/unit_test.hpp"

#include "indel_util.hh"


BOOST_AUTO_TEST_SUITE( indel_util )


BOOST_AUTO_TEST_CASE( test_indel_intersect_deletion ) {

    // 5-base deletion starting at zero-indexed position 9:
    const indel_key ikdel1(9,INDEL::DELETE,5);

    // left-side deletion:
    BOOST_CHECK(! is_indel_conflict(ikdel1,indel_key(7,INDEL::DELETE,1)) );
    BOOST_CHECK(  is_indel_conflict(ikdel1,indel_key(8,INDEL::DELETE,1)) );

    // right-side deletion:
    BOOST_CHECK(  is_indel_conflict(ikdel1,indel_key(14,INDEL::DELETE,1)) );
    BOOST_CHECK(! is_indel_conflict(ikdel1,indel_key(15,INDEL::DELETE,1)) );

    // super-set deletion
    BOOST_CHECK(  is_indel_conflict(ikdel1,indel_key(4,INDEL::DELETE,15)) );

    // sub-set deletion
    BOOST_CHECK(  is_indel_conflict(ikdel1,indel_key(12,INDEL::DELETE,1)) );
}


BOOST_AUTO_TEST_SUITE_END()

