
// no unit test infrastructure in place to run this yet...
#if 0
// hack compile: g++ test_align_path.cpp align_path.cpp ../blt_util/blt_exception.cpp ../blt_util/parse_util.cpp ../blt_util/log.cpp ../blt_util/seq_util.cpp -I. -I.. -I/home/csaunders/opt/x86_64-linux/boost-1.49.0/include

#include "align_path.hh"

#include <cstdlib>

#include <iostream>


void
test_string_clean(const char* cigar, const char* expect) {
    using namespace ALIGNPATH;
    path_t apath;
    cigar_to_apath(cigar,apath);
    apath_cleaner(apath);

    path_t expect_path;
    cigar_to_apath(expect,expect_path);
    std::cerr << "cigar in: " << cigar << " apath_out: " << apath << " expect_out: " << expect << "\n";

    const bool is_pass(apath==expect_path);
    if(! is_pass) exit(EXIT_FAILURE);
}


int
main() {
    test_string_clean("29M2I5M1D0M3I20M16S","29M2I5M1D3I20M16S");
    test_string_clean("1M1P1M","2M");
    test_string_clean("1M1D0M1I1M","1M1D1I1M");
    test_string_clean("0H1H1S0I1M1D0M1I1I1D1I1M1D1M","29M2I5M1D3I20M16S");
}
#endif
