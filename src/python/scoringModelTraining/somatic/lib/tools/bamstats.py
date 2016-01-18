#!/illumina/development/haplocompare/hc-virtualenv/bin/python
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2016 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

import pysam
import pandas


def bamStats(bamfile):
    """ Extract average depths + idxstats data from BAM file, return data frame """
    istats = pysam.idxstats(bamfile)
    result = []
    samfile = pysam.Samfile(bamfile, "rb")
    for x in istats:
        xs = x.replace("\n", "").split("\t")
        rec = {
            "CHROM": xs[0],
            "NT": int(xs[1]),
            "MAPPED": int(xs[2]),
            "UNMAPPED": int(xs[3]),
            "READLEN": 0,
            "COVERAGE": 0.0,
        }
        count = 0
        rls = 0.0
        try:
            for read in samfile.fetch(xs[0]):
                rls += float(read.rlen)
                count += 1
                if count > 10000:
                    break

            rls /= count
            rec["READLEN"] = rls
            rec["COVERAGE"] = float(rec["MAPPED"] * rec["READLEN"])/float(rec["NT"])
        except:
            pass
        result.append(rec)

    if result:
        return pandas.DataFrame(result, columns=["CHROM", "NT", "MAPPED", "UNMAPPED", "READLEN", "COVERAGE"])
    else:
        return pandas.DataFrame(columns=["CHROM", "NT", "MAPPED", "UNMAPPED", "READLEN", "COVERAGE"])
