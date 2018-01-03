//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
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
                if (*p2=='.')
                {
                    ++p;
                    continue;
                }

                sn.noise++;
                while (true)
                {
                    if (*p2=='\0' || *p2=='\n' || *p2=='\t') break;
                    if (*p2==':')
                    {
                        const char* p3(p2+1);
                        while (true)
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

