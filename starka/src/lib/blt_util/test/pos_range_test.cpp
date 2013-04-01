#include "boost/test/unit_test.hpp"

#include "pos_range.hh"


BOOST_AUTO_TEST_SUITE( test_pos_range )


BOOST_AUTO_TEST_CASE( test_pos_range_is_pos_intersect ) {

    // this corresponds to zero-index range [9,19] :
    const pos_range pr(9,20);

    BOOST_CHECK(! pr.is_pos_intersect(8));
    BOOST_CHECK(  pr.is_pos_intersect(9));
    BOOST_CHECK(  pr.is_pos_intersect(19));
    BOOST_CHECK(! pr.is_pos_intersect(20));
}


BOOST_AUTO_TEST_CASE( test_pos_range_semibound_is_pos_intersect ) {

    // this corresponds to zero-index range [-inf,19] :
    pos_range pr;
    pr.set_end_pos(20);

    BOOST_CHECK(  pr.is_pos_intersect(8));
    BOOST_CHECK(  pr.is_pos_intersect(9));
    BOOST_CHECK(  pr.is_pos_intersect(19));
    BOOST_CHECK(! pr.is_pos_intersect(20));
}

BOOST_AUTO_TEST_CASE( test_pos_range_is_range_intersect ) {

    // this corresponds to zero-index range [9,19] :
    const pos_range pr(9,20);

    // left-side:
    BOOST_CHECK(! pr.is_range_intersect(pos_range(0,9)));
    BOOST_CHECK(  pr.is_range_intersect(pos_range(0,10)));

    // right side:
    BOOST_CHECK(  pr.is_range_intersect(pos_range(19,30)));
    BOOST_CHECK(! pr.is_range_intersect(pos_range(20,30)));

    // superset:
    BOOST_CHECK(  pr.is_range_intersect(pos_range(0,30)));

    // subset:
    BOOST_CHECK(  pr.is_range_intersect(pos_range(12,15)));
}

BOOST_AUTO_TEST_SUITE_END()

