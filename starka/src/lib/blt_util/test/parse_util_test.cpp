#include "boost/test/unit_test.hpp"

#include "parse_util.hh"

#include <string>

BOOST_AUTO_TEST_SUITE( parse_util )

using namespace casava::blt_util;

BOOST_AUTO_TEST_CASE( test_parse_int ) {
    const char* two = "2";
    const int val(parse_int(two));
    BOOST_CHECK_EQUAL(val, 2);
}

BOOST_AUTO_TEST_CASE( test_parse_int_str ) {
    static const char two[] = "2";
    const int val(parse_int_str(std::string(two)));
    BOOST_CHECK_EQUAL(val, 2);
}

BOOST_AUTO_TEST_CASE( test_parse_int_str_bad_input ) {
    static const std::string junk("ABCD");
    BOOST_CHECK_THROW(parse_int_str(junk), std::exception);
}

BOOST_AUTO_TEST_SUITE_END()

