#!/illumina/development/haplocompare/hc-virtualenv/bin/python
#
# Starka
# Copyright (c) 2009-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
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
