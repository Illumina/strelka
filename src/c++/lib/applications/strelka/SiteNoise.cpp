// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "SiteNoise.hh"
#include "blt_util/parse_util.hh"


void
set_noise_from_vcf(
    const char* line,
    SiteNoise& sn)
{
    const char* p(line);
    if (p == nullptr) return;

    unsigned column=0;
    while (true)
    {
        if (*p=='\0' || *p=='\n') break;

        if (*p=='\t')
        {
            column++;
            if (column>=9)
            {
                sn.total++;
                const char* p2(p+1);
                if(*p2=='.')
                {
                    ++p;
                    continue;
                }

                sn.noise++;
                while(true)
                {
                    if (*p2=='\0' || *p2=='\n' || *p2=='\t') break;
                    if (*p2==':')
                    {
                        const char* p3(p2+1);
                        while(true)
                        {
                            if (*p3=='\0' || *p3=='\n' || *p3=='\t' || *p3==':') break;
                            if (*p3==',')
                            {
                                const char* p4(p3+1);
                                unsigned alt=illumina::blt_util::parse_unsigned(p4);
                                if (alt>1) sn.noise2++;
                                break;
                            }
                            p3++;
                        }
                        break;
                    }
                    p2++;
                }
            }
        }
        ++p;
    }
}

