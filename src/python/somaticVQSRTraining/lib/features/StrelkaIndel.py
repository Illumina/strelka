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

import pandas
import logging

from vqsr.tools.vcfextract import vcfExtract, extractHeaders
from . import FeatureSet


@FeatureSet.register("strelka.indel")
class StrelkaAdmixSNVFeatures(FeatureSet):
    """ Collect SNV features from Strelka-to-admixture comparison """

    def collect(self, vcfname):
        """ Return a data frame with features collected from the
            given VCF, tagged by given type """
        features = ["CHROM", "POS", "REF", "ALT", "FILTER",
                    "I.NT", "I.SOMATIC", "I.QSI_NT",
                    "I.SGT", "I.RC", "I.RU", "I.IC", "I.IHP",
                    "I.MQ", "I.MQ0",
                    "I.H200", "I.RC_HPOL_200", "I.RC_DINUC_200", "I.RC_TRIPLET_200",
                    "S.1.DP", "S.2.DP",
                    "S.1.TAR", "S.2.TAR",
                    "S.1.TIR", "S.2.TIR",
                    "S.1.TOR", "S.2.TOR",
                    "S.1.DP50", "S.2.DP50",
                    "S.1.FDP50", "S.2.FDP50",
                    "S.1.SUBDP50", "S.2.SUBDP50"]

        records = []

        avg_depth = self.chr_depth
        if not avg_depth:
            avg_depth = {}

            for l in list(extractHeaders(vcfname)):
                x = str(l).lower()
                if '##maxdepth_' in x:
                    xl = str(l).split('=')
                    xchr = xl[0][11:]
                    avg_depth[xchr] = float(xl[1])
                    # logging.info("Maxdepth for %s depth from VCF header is %f" % (xchr, avg_depth[xchr]))

        has_warned = {}
        for vr in vcfExtract(vcfname, features):
            rec = {}
            for i, ff in enumerate(features):
                rec[ff] = vr[i]

            # fix missing features
            for q in ["I.QSI_NT", "I.RC", "I.IC", "I.IHP",
                      "S.1.DP", "S.2.DP", "I.H200", "I.RC_HPOL_200",
                      "I.RC_DINUC_200", "I.RC_TRIPLET_200",
                      "S.1.FDP50", "S.2.FDP50",
                      "S.1.SUBDP50", "S.2.SUBDP50"]:
                if q not in rec or rec[q] is None:
                    rec[q] = 0
                    if not ("feat:" + q) in has_warned:
                        logging.warn("Missing feature %s" % q)
                        has_warned["feat:" + q] = True

            for q in ["S.1.TAR", "S.2.TAR",
                      "S.1.TIR", "S.2.TIR",
                      "S.1.TOR", "S.2.TOR"]:
                if q not in rec or rec[q] is None:
                    rec[q] = [0, 0]
                    if not ("feat:" + q) in has_warned:
                        logging.warn("Missing feature %s" % q)
                        has_warned["feat:" + q] = True

            NT = rec["I.NT"]
            NT_is_ref = int(NT == "ref")
            QSI_NT = int(rec["I.QSI_NT"])

            n_D_total_1 = float(
                rec["S.1.TIR"][0]) + float(rec["S.1.TAR"][0]) + float(rec["S.1.TOR"][0])
            t_D_total_1 = float(
                rec["S.2.TIR"][0]) + float(rec["S.2.TAR"][0]) + float(rec["S.2.TOR"][0])
            n_D_total_2 = float(
                rec["S.1.TIR"][1]) + float(rec["S.1.TAR"][1]) + float(rec["S.1.TOR"][1])
            t_D_total_2 = float(
                rec["S.2.TIR"][1]) + float(rec["S.2.TAR"][1]) + float(rec["S.2.TOR"][1])

            n_TOR_ratio_1 = float(
                rec["S.1.TOR"][0]) / n_D_total_1 if n_D_total_1 != 0 else 0
            t_TOR_ratio_1 = float(
                rec["S.2.TOR"][0]) / t_D_total_1 if t_D_total_1 != 0 else 0
            n_TOR_ratio_2 = float(
                rec["S.1.TOR"][1]) / n_D_total_2 if n_D_total_2 != 0 else 0
            t_TOR_ratio_2 = float(
                rec["S.2.TOR"][1]) / t_D_total_2 if t_D_total_2 != 0 else 0

            n_DP = float(rec["S.1.DP"])
            t_DP = float(rec["S.2.DP"])

            in_del = 0

            max_len = len(rec["REF"])
            min_len = len(rec["REF"])

            for a in rec["ALT"]:
                if len(a) > len(rec["REF"]):
                    in_del |= 1
                else:
                    in_del |= 2
                min_len = min(len(a), min_len)
                max_len = max(len(a), max_len)

            ilen = max_len - min_len

            n_DP_ratio = 0
            t_DP_ratio = 0

            if avg_depth:
                if rec["CHROM"] in avg_depth:
                    n_DP_ratio = n_DP / float(avg_depth[rec["CHROM"]])
                    t_DP_ratio = t_DP / float(avg_depth[rec["CHROM"]])
                elif not rec["CHROM"] in has_warned:
                    logging.warn("Cannot normalize depths on %s" % rec["CHROM"])
                    has_warned[rec["CHROM"]] = True
            elif "DPnorm" not in has_warned:
                logging.warn("Cannot normalize depths.")
                has_warned["DPnorm"] = True

            # Ref and alt allele counts for tier1 and tier2
            t_allele_ref_counts = map(float, rec['S.2.TAR'])
            t_allele_alt_counts = map(float, rec['S.2.TIR'])

            # Compute the tier1 and tier2 alt allele rates.
            if t_allele_alt_counts[0] + t_allele_ref_counts[0] == 0:
                t_tier1_allele_rate = 0
            else:
                t_tier1_allele_rate = t_allele_alt_counts[
                    0] / float(t_allele_alt_counts[0] + t_allele_ref_counts[0])

            if t_allele_alt_counts[1] + t_allele_ref_counts[1] == 0:
                t_tier2_allele_rate = 0
            else:
                t_tier2_allele_rate = t_allele_alt_counts[
                    1] / float(t_allele_alt_counts[1] + t_allele_ref_counts[1])

            # Ref and alt allele counts for tier1 and tier2
            n_allele_ref_counts = map(float, rec['S.1.TAR'])
            n_allele_alt_counts = map(float, rec['S.1.TIR'])

            # Compute the tier1 and tier2 alt allele rates.
            if n_allele_alt_counts[0] + n_allele_ref_counts[0] == 0:
                n_tier1_allele_rate = 0
            else:
                n_tier1_allele_rate = n_allele_alt_counts[
                    0] / float(n_allele_alt_counts[0] + n_allele_ref_counts[0])

            if n_allele_alt_counts[1] + n_allele_ref_counts[1] == 0:
                n_tier2_allele_rate = 0
            else:
                n_tier2_allele_rate = n_allele_alt_counts[
                    1] / float(n_allele_alt_counts[1] + n_allele_ref_counts[1])

            bcn = 0

            try:
                bcn = rec["S.1.FDP50"] / rec["S.1.DP50"]
            except:
                pass

            try:
                bcn = max(bcn, rec["S.2.FDP50"] / rec["S.2.DP50"])
            except:
                pass

            # Gather the computed data into a dict
            qrec = {
                "CHROM": rec["CHROM"],
                "POS": int(rec["POS"]),
                "REF": rec["REF"],
                "ALT": ",".join(rec["ALT"]),
                "LENGTH": ilen,
                "LENGTHGT5": 0 if ilen <= 5 else 1,
                "INDELTYPE": in_del,
                "FILTER": ",".join(rec["FILTER"]),
                "NT": NT,
                "NT_REF": NT_is_ref,
                "QSI_NT": QSI_NT,
                "N_TOR_RATE_TIER1": n_TOR_ratio_1,
                "N_TOR_RATE_TIER2": n_TOR_ratio_2,
                "T_TOR_RATE_TIER1": t_TOR_ratio_1,
                "T_TOR_RATE_TIER2": t_TOR_ratio_2,
                "N_DP": n_DP,
                "T_DP": t_DP,
                "N_DP_RATE": n_DP_ratio,
                "T_DP_RATE": t_DP_ratio,
                "T_TIER1_ALT_RATE": t_tier1_allele_rate,
                "T_TIER2_ALT_RATE": t_tier2_allele_rate,
                "N_TIER1_ALT_RATE": n_tier1_allele_rate,
                "N_TIER2_ALT_RATE": n_tier2_allele_rate,
                "SGT": rec["I.SGT"],
                "entropy": rec["I.H200"],
                "hpol": rec["I.RC_HPOL_200"],
                "dinuc": rec["I.RC_DINUC_200"],
                "triplet": rec["I.RC_TRIPLET_200"],
                "bcn": bcn
            }

            try:
                qrec["RC"] = int(rec["I.RC"])
            except:
                qrec["RC"] = 0

            try:
                qrec["RU"] = rec["I.RU"]
            except:
                qrec["RU"] = ""

            try:
                qrec["RU_LEN"] = len(rec["I.RU"])
            except:
                qrec["RU_LEN"] = 0

            try:
                qrec["IC"] = int(rec["I.IC"])
            except:
                qrec["IC"] = 0

            try:
                qrec["IHP"] = int(rec["I.IHP"])
            except:
                qrec["IHP"] = 0

            try:
                qrec["S.1.FDP50"] = float(rec["S.1.FDP50"])
            except:
                qrec["S.1.FDP50"] = 0

            try:
                qrec["S.2.FDP50"] = float(rec["S.2.FDP50"])
            except:
                qrec["S.2.FDP50"] = 0

            try:
                qrec["S.1.SUBDP50"] = float(rec["S.1.SUBDP50"])
            except:
                qrec["S.1.SUBDP50"] = 0

            try:
                qrec["S.2.SUBDP50"] = float(rec["S.2.SUBDP50"])
            except:
                qrec["S.2.SUBDP50"] = 0

            try:
                qrec["MQ"] = float(rec["I.MQ"])
            except:
                qrec["MQ"] = 0

            try:
                qrec["MQ0"] = float(rec["I.MQ0"])
            except:
                qrec["MQ0"] = 0

            records.append(qrec)

        cols = ["CHROM",
                "POS",
                "REF",
                "ALT",
                "LENGTH",
                "LENGTHGT5",
                "INDELTYPE",
                "FILTER",
                "NT",
                "NT_REF",
                "QSI_NT",
                "N_TOR_RATE_TIER1",
                "T_TOR_RATE_TIER1",
                "N_DP",
                "T_DP",
                "N_DP_RATE",
                "T_DP_RATE",
                "T_TIER1_ALT_RATE",
                "T_TIER2_ALT_RATE",
                "N_TIER1_ALT_RATE",
                "N_TIER2_ALT_RATE",
                "SGT",
                "RC",
                "RU",
                "RU_LEN",
                "IC",
                "IHP",
                "S.1.FDP50",
                "S.1.SUBDP50",
                "MQ",
                "MQ0",
                "bcn"
                ]

        if records:
            df = pandas.DataFrame(records, columns=cols)
        else:
            df = pandas.DataFrame(columns=cols)

        return df

    def trainingfeatures(self):
        """ Return a list of columns that are features to use for VQSR training """
        return ["QSI_NT",
                "TIER1_ALLELE_RATE",
                "RC",
                "IHP",
                "MQ",
                "bcn",
                "T_TOR_RATE_TIER1"]
