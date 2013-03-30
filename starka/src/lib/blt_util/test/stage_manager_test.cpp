#include "boost/test/unit_test.hpp"

#include "stage_manager.hh"


BOOST_AUTO_TEST_SUITE( test_stage_manager )


// create a standard stage size arrangement for testing:
//
// This returns a tree with:
//   stage 1 following 10 bases behind the root
//   stage 2 following 20 bases behind stage 1
//   stage 3 following 20 bases behind the root
//
static
stage_data
get_test_stage_data() {

    stage_data sd;
    sd.add_stage(0);
    sd.add_stage(1,0,10);
    sd.add_stage(2,1,20);
    sd.add_stage(3,0,20);

    return sd;
}


BOOST_AUTO_TEST_CASE( test_stage_data_dist ) {

    const stage_data sd(get_test_stage_data());

    BOOST_CHECK_EQUAL(sd.get_stage_id_shift(0),0);
    BOOST_CHECK_EQUAL(sd.get_stage_id_shift(1),10);
    BOOST_CHECK_EQUAL(sd.get_stage_id_shift(2),30);
    BOOST_CHECK_EQUAL(sd.get_stage_id_shift(3),20);
}


BOOST_AUTO_TEST_SUITE_END()

