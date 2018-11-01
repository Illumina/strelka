#!/usr/bin/env python2
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2018 Illumina, Inc.
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

from __future__ import print_function

import argparse
import csv
import errno
import itertools
import os
import re
import sys
import distutils.spawn
import multiprocessing
import copy

import subprocess as sp

import cyvcf2 as vcf
import numpy as np

# check for python 2.7
if sys.version_info[:2] != (2, 7):
    msg = "Python 2.7 is required for execution."  # pragma: no cover
    raise ImportError(msg)  # pragma: no cover

# constants

# 0) program name and version

_name = 'denovo.py'
_version = "0.7.0"


# 1) output variables

_field_name = 'DQ'
_sample_filter_name = 'FT'
_mendel_conflict_name = 'ME'

# upper limit for DQ score
_max_score = np.float(100)

# number of reported decimals for DQ score
_default_dq_decimals = 1

# header line definition
_dq_format = {
    'ID': _field_name,
    'Number': 1,
    'Type': 'Float',
    'Description': 'Denovo score'
}


# 2) input requirements

_pedigree_names = ('proband', 'mother', 'father')

_required_format_fields = ('PL', )


# 3) model specifications and defaults

_default_model_name = 'dng'

_select_methods = ('call', 'best')
_default_select_method = _select_methods[0]

# filter settings only used on allosomes in a male proband
_min_read_depth = 20  # minimum read depth for diploid chromosome
_min_gqx = 10

_indel_mrate = 1e-5
_indel_mu_scale = 20.0


# 4) general

_indel_gt_types = {
    (True, True): 'RR',
    (False, True): 'RV',
    (True, False): 'RV',
    (False, False): 'VV'
}

_default_params = {
    'model': _default_model_name,
    'min_read_depth': _min_read_depth,
    'min_gqx': _min_gqx,
    'pass_sample_filters': True,
    'select_method': 'call',
    'score_decimals': _default_dq_decimals,
    'print_vcf_header': True,
    'print_vcf_body': True,
    'filter_pedphase': True,
    'par_regions': None
}


# functions
def get_data_path(file_name, dir_name='', space_name='data', check=True):
    """Compute full path to a distributed file"""

    script_dir = os.path.abspath(os.path.dirname(__file__))
    rel_path = os.path.join(space_name, dir_name, file_name)
    path = os.path.join(script_dir, rel_path)

    if check and not os.path.exists(path):
        msg = "File '{0}' does not exist.".format(path)
        raise IOError(msg)

    return path


def is_field_in_vcf_header(vf, field_name, header_section):
    """Check if field is defined in section of the VCF header"""

    allowed_header_sections = {'GENERIC', 'CONTIG', 'INFO', 'FORMAT', 'FILTER', 'STR'}

    if header_section not in allowed_header_sections:
        raise ValueError("Field '{}' not defined in the '{}' section of the VCF header".format(field_name, header_section))

    for header_line in vf.header_iter():
        header_line = header_line.info()
        if header_line.get('HeaderType') == header_section and header_line.get('ID') == field_name:
            return True

    return False


_prior_types = ('auto', 'xx', 'xy')
_var_types = ('snv', 'indel')
_prior_paths = {t: {v: get_data_path('{}_lookup_{}.tsv'.format(v, t), 'prior') for v in _var_types} for t in _prior_types}


def parse_pedigree_from_vcf_header(vcf):
    """Extract the pedigree sample names from SPW VCF header"""

    regex = re.compile('^##PEDIGREE=<(.*)>$')
    for line in vcf.raw_header.split('\n'):
        match = regex.match(line)
        if match:
            value = match.groups()[0]
            pedigree = (s.split('=') for s in value.split(','))
            pedigree = dict((k.lower(), v) for k, v in pedigree)
            return pedigree


def parse_karyotypes_from_vcf_header(vcf, valid_karyotypes=('XX', 'XY')):
    """Extract predicted samples karyotypes from SPW VCF header"""

    regex = re.compile('##PredictedSexChromosomeKaryotype(.+)=(.{2})')
    karyotypes = dict()
    for line in vcf.raw_header.split('\n'):
        match = regex.match(line)
        if match:
            sample, karyo = match.groups()
            if karyo not in valid_karyotypes:
                msg = "Invalid karyotype '{}' for '{}'".format(sample, karyo)  # pragma: no cover
                raise ValueError(msg)  # pragma: no cover
            karyotypes[sample] = karyo

    karyotypes = karyotypes if karyotypes else None

    return karyotypes


def paste(elements, sep=''):
    """Simplified joining of objects with conversion"""

    return sep.join(map(str, elements))


def select_output(out_vcf, mode='w', default=sys.stdout):
    """Select the output stream for the final VCF"""

    output_vcf = open(out_vcf, mode) if out_vcf is not None else default

    return output_vcf


def check_pedigree_input(pedigree, samples):
    """Check if pedigree inputs are consistent and match the VCF"""

    for k in _pedigree_names:
        if k not in pedigree:
            msg = "Pedigree information for {0} is missing.".format(k)
            raise ValueError(msg)

    # check if all pedigree sample names are referenced in the header
    for k, v in pedigree.iteritems():
        if v not in samples:
            msg = "Sample '{1}' for {0} missing in the input VCF.".format(k, v)
            raise ValueError(msg)

    return pedigree


def check_vcf_input(vcf_in):
    """Check for valid VCF input and if all required fields are present"""

    # open VCF stream
    try:
        vf_in = vcf.VCF(vcf_in)
    except BaseException:
        msg = "Invalid VCF file or input stream."
        raise IOError(msg)

    format_fields = {x['ID']: x for x in vf_in.header_iter() if x.type == 'FORMAT'}

    for field in _required_format_fields:
        if field not in format_fields:
            msg = "Required format field {0} missing in VCF.".format(field)
            raise ValueError(msg)

    return vf_in


def check_param_input(param, samples, pedigree):
    """Check if model inputs are consistent"""

    if param['model'] == 'dng' and param['min_read_depth'] < 1:
        msg = "The minimum read depth per sample must be positive."
        raise ValueError(msg)

    if param['select_method'] not in _select_methods:
        msg = "Invalid select method {0}".format(param['select_method'])
        raise ValueError(msg)

    param['is_select_call'] = param['select_method'] == 'call'

    param['sample_index'] = [samples.index(pedigree[s])
                             for s in _pedigree_names]
    param['n_samples'] = len(samples)
    param['sample_index_chrX'] = [param['sample_index'][0], param['sample_index'][1]]  # w/o father
    param['sample_index_chrY'] = [param['sample_index'][0], param['sample_index'][2]]  # w/o mother

    return param


def add_denovo_vcf_header_lines(vf, dq_format=_dq_format, name=_name, version=_version):
    """Add denovo and DQ header lines to the VCF"""

    new_header_line = '##denovo_program={} {}'.format(name, version)

    existing_header_lines = vf.raw_header.splitlines()

    # only add header lines if they do not already exist
    for line in existing_header_lines:
        if line == new_header_line:
            break
        vf.add_to_header(new_header_line)

    if 'DQ' not in vf:
        vf.add_format_to_header(dq_format)

    return vf


def denovo(in_vcf, pedigree, param=_default_params, out_vcf=None):
    """Calculates DQ score for each variant, returns VCF stream"""

    vf_in = check_vcf_input(in_vcf)

    karyotypes = parse_karyotypes_from_vcf_header(vf_in)

    pedigree = check_pedigree_input(pedigree, vf_in.samples)

    param = check_param_input(param, vf_in.samples, pedigree)

    model_fun = _models[param['model']]
    if param['model'] is 'dng':
        prior = {}
        prior['auto'] = {
            'snv': reorder_prior(read_prior(_prior_paths['auto']['snv'], 'snv'), _full_gt_idx_snv),
            'indel': reorder_prior(read_prior(_prior_paths['auto']['indel'], 'indel'), _full_gt_idx_indel)
        }
        if karyotypes is not None and pedigree['proband'] in karyotypes:
            gender = karyotypes[pedigree['proband']].lower()
            param['proband_karyotype'] = gender
            prior['chrX'] = {
                'snv': reorder_prior(read_prior(_prior_paths[gender]['snv'], 'snv'), _full_gt_idx_snv),
                'indel': reorder_prior(read_prior(_prior_paths[gender]['indel'], 'indel'), _full_gt_idx_indel)
            }
        else:
            # cannot decide if male or female -> fallback: use autosomal priors
            param['proband_karyotype'] = None
            prior['chrX'] = {
                'snv': prior['auto']['snv'],
                'indel': prior['auto']['indel']
            }
    else:
        prior = None

    param['female_proband'] = param.get('proband_karyotype') == 'xx'
    param['male_proband'] = param.get('proband_karyotype') == 'xy'
    param['filter_pedphase'] = param.get('filter_pedphase', True) and \
                               is_field_in_vcf_header(vf_in, _mendel_conflict_name, 'FORMAT')
    param['pass_sample_filters'] = param.get('pass_sample_filters', True) and \
                                   is_field_in_vcf_header(vf_in, _sample_filter_name, 'FORMAT')

    decimals = param.get('score_decimals', _default_dq_decimals)
    print_vcf_header = param.get('print_vcf_header', _default_params['print_vcf_header'])
    print_vcf_body = param.get('print_vcf_body', _default_params['print_vcf_body'])

    # template for output VCF: add header lines
    vf_in = add_denovo_vcf_header_lines(vf_in, _dq_format, _name, _version)

    output_file = select_output(out_vcf)

    try:
        if print_vcf_header:
            print(vf_in.raw_header, end='', file=output_file)
        if print_vcf_body:
            for variant in vf_in:
                score = model_fun(variant, prior, param)
                variant_entry = add_DQ_score_to_variant(variant, param, score)
                print(variant_entry, end='', file=output_file)
    # handle premature closing of stdout
    except IOError as e:
        if e.errno == errno.EPIPE:
            pass

    if out_vcf is not None:
        output_file.close()

    return True


# dummy model

def calculate_dummy_DQ(variant, prior, param):
    """Return a score of 0.0 for each variant, without any calculation"""

    score = 0.0

    return score


# DNG model

def bases_index(bases):
    """Map single sample genotypes to DNG indexing schema"""

    gts = list(itertools.product(bases, bases))
    gtsu = [y for y in gts if y[0] <= y[1]]
    gtsu_idx = {gt: i for i, gt in enumerate(gtsu)}
    gts_idx = {gt: gtsu_idx[tuple(sorted(gt))] for gt in gts}

    return gts_idx, gtsu_idx


def full_genotype_mapping(gtsu_idx):
    """Map trio sample genotypes to DNG indexing schema"""

    n = len(gtsu_idx)
    full_idx = itertools.product(xrange(n), xrange(n), xrange(n))
    idx_gtsu = {v: k for k, v in gtsu_idx.iteritems()}  # index->GT mapping
    full_gt = [
        "/".join([paste(idx_gtsu[k]), paste(idx_gtsu[j]), paste(idx_gtsu[i])])
        for i, j, k in full_idx
    ]
    full_gt = {k: i for i, k in enumerate(full_gt)}

    return full_gt


def genotype_maps(bases):
    """Construct mappings of genotypes to internal representation"""

    gts_idx, gtsu_idx = bases_index(bases)
    gts_gtsu = {
        sep.join(g): paste(sorted(g))
        for g in list(itertools.product(bases, bases)) for sep in ["/", "|"]
    }
    full_gt_idx = full_genotype_mapping(gtsu_idx)

    return gts_idx, gtsu_idx, gts_gtsu, full_gt_idx


def alleles_integers(variant):
    """Returns the integer representation of the VCF alleles:
    0: REF, 1,2,...,N: ALTs, -1:missing"""

    gt = variant.format("GT", int)
    if gt is not None:
        gt = gt / 2 - 1

    return gt


def all_genotypes_complete(variant, sample_indices):
    """Returns true if all sample genotypes are complete,
    i.e. contain no missing alleles"""

    gt = variant.format('GT', int)
    if gt is None or np.any(gt[sample_indices] < 0):
        return False
    else:
        return True


# snv
_alleles_snv = ('A', 'C', 'G', 'T')
_gts_idx_snv, _gtsu_idx_snv, _gts_gtsu_snv, _full_gt_idx_snv = genotype_maps(_alleles_snv)
_snv_possible_genotypes = tuple(_gts_gtsu_snv.keys())

for a in _alleles_snv:
    _gts_gtsu_snv[a] = a*2


# indel
_alleles_indel = ('R', 'V')  # R: ref, V: variant
_gts_idx_indel, _gtsu_idx_indel, _gts_gtsu_indel, _full_gt_idx_indel = genotype_maps(_alleles_indel)


# maximum indel length in prior lookup table
_max_indel_len = 200

# parameters for length-dependent indel prior
_indel_prior_params = {
    'ins': {'slope': -0.2994, 'intercept': -22.8689},
    'del': {'slope': -0.2856, 'intercept': -21.9313}
}

_indel_mrate_scale = 1.0


def precompute_indel_prior(param, max_length):
    """Precompute a lookup table for a length-dependent indel prior"""

    dx = np.arange(0, max_length+1)
    mut_rate = np.exp(dx * param['slope'] + param['intercept'])

    return mut_rate


def read_prior(path, vartype):
    """Read prior lookup tables generated by DNG"""

    if vartype.lower() == 'snv':
        col_names = ('n_u_alleles', 'case', 'transmission_prob', 'gt', 'mrate',
                     'denovo_flag', 'normal_flag', 'A', 'C', 'G', 'T')
        col_formats = ('i4', 'i4', 'f4', 'S8', 'f4', 'i4', 'i4', 'f4', 'f4',
                       'f4', 'f4')
    elif vartype.lower() == 'indel':
        col_names = ('n_u_alleles', 'case', 'transmission_prob', 'gt', 'mrate',
                     'denovo_flag', 'normal_flag', 'prior')
        col_formats = ('i4', 'i4', 'f4', 'S8', 'f4', 'i4', 'i4', 'f4')
    else:
        msg = "Unknown variant type for prior: " + vartype
        raise ValueError(msg)
    assert len(col_names) == len(col_formats)

    prior = np.loadtxt(
        path, unpack=True, dtype={
            'names': col_names,
            'formats': col_formats
        })
    prior = dict(zip(col_names, prior))

    prior['denovo_flag'].astype(bool)

    if vartype.lower() in 'snv':
        pp = prior['transmission_prob']
        pp *= prior['mrate']
        for ref in _alleles_snv:
            prior[ref] = prior[ref] * pp
    elif vartype.lower() in 'indel':
        for indel_type in ('ins', 'del'):
            mrate = precompute_indel_prior(_indel_prior_params[indel_type], _max_indel_len)
            pp = prior['transmission_prob']
            pp *= prior['prior']
            ppx = np.outer(pp, mrate * _indel_mu_scale)
            n_alleles = np.tile(prior['n_u_alleles'].reshape(-1, 1), _max_indel_len+1)
            ppx = np.power(ppx, n_alleles)
            assert ppx.shape == (27, _max_indel_len+1)
            prior[indel_type] = ppx
            assert prior[indel_type].shape == (27, _max_indel_len+1)

    return prior


def reorder_prior(prior, full_gt_idx):
    """Reorder imported prior to match the DNG genotype order"""

    p_gt = {g: i for i, g in enumerate(prior['gt'])}
    idx_full_gt = {v: k for k, v in full_gt_idx.iteritems()}  # revert mapping
    new_prior_ord = [p_gt[idx_full_gt[i]] for i in xrange(len(idx_full_gt))]

    prior = {k: np.take(v, new_prior_ord, axis=0) for k, v in prior.iteritems()}

    return prior


def get_PL_field(variant):
    """Extract and rectify PL field from a variant record, return as numpy array"""

    pl = variant.format('PL', int)
    pl = rectify_pl_field(pl)

    return pl


def map_sample_pl(pl, alleles_idx):
    """Map sample PL field to DNG genotypes, return probabilities"""

    g = np.zeros(10)
    g[alleles_idx] = _score2prob_np[pl]

    return g


def map_sample_pl_indel(pl, alleles_idx):
    """Map sample PL field to DNG genotypes, return probabilities"""

    g = np.zeros(3)
    g[alleles_idx] = _score2prob_np[pl]

    return g


def allele_indices(alleles, genotype_index):
    """Return ordered indices for allele combinations for PL field"""

    indices = [genotype_index[alleles[j], alleles[i]]
               for i in xrange(len(alleles)) for j in xrange(i + 1)]

    return indices


def build_PL_indices_lookup(alleles, genotype_index, ref=None):
    """Build lookup table for mapping of alleles to PL field"""

    n_uniq_alleles = len(set(alleles))
    alleles_comb = itertools.chain(
        itertools.product(
            alleles, repeat=2),
        itertools.product(
            alleles, repeat=3),
        itertools.product(
            alleles, repeat=4))
    alleles_comb = itertools.ifilter(
        lambda x: len(set(x)) == min(len(x), n_uniq_alleles), alleles_comb)
    alleles_comb = itertools.ifilter(lambda x: ref not in x[1:], alleles_comb)
    alleles_lookup = {
        alleles: allele_indices(alleles, genotype_index)
        for alleles in alleles_comb
    }

    return alleles_lookup


def get_PL_indices(alleles, lookup, genotype_index):
    """Lookup mapping of alleles to PL field"""

    return lookup[alleles] if alleles in lookup else allele_indices(
        alleles, genotype_index)


def get_PL_indices_snv(alleles):
    """Lookup mapping of SNV alleles to PL field"""

    return get_PL_indices(alleles, _PL_indices_snv, _gts_idx_snv)


def get_PL_indices_indel(alleles):
    """Lookup mapping of indel alleles to PL field"""

    return get_PL_indices(alleles, _PL_indices_indel, _gts_idx_indel)


_PL_indices_snv = build_PL_indices_lookup(_alleles_snv, _gts_idx_snv)

_PL_indices_indel = build_PL_indices_lookup(_alleles_indel, _gts_idx_indel, 'R')


def score2prob(x):
    """Convert score to probability"""

    return np.power(10.0, -x / 10.0)


score2prob_vec = np.vectorize(score2prob)


def prob2score(x, max_score=_max_score):
    """Converts probability to score, optimised for speed"""

    if x > 0.0:
        y = np.abs(-10.0 * np.log10(x))
        if y > max_score:
            return max_score
        else:
            return y
    else:
        return max_score


prob2score_vec = np.vectorize(prob2score)

_score2prob_np = score2prob_vec(np.arange(0, 10001))


def is_denovo_candidate(variant, param, mendel_conflict_name=_mendel_conflict_name):
    """Check if a variant is a potential denovo candidate"""

    # use pedigree phasing and Mendelian conflict information
    if param['filter_pedphase']:
        sample_index_proband = param['sample_index'][0]
        try:
            # a potential denovo variants must have ME != 0
            # missing is encoded as typemin(int32): also consider as a candidate
            return variant.format(mendel_conflict_name)[sample_index_proband][0] != 0
        # safeguard if format field is missing or invalid
        except KeyError:
            pass

    # treat each variant as a candidate unless proven otherwise
    return True


def calculate_dng_DQ_male_allosome(variant, chrom, param):
    """Calculate denovo score for allosomal variants in a male proband"""

    if variant.is_snp and not can_compute_dng_DQ_snv(variant, param):
        return None

    if variant.is_indel and not can_compute_dng_DQ_indel(variant, param):
        return None

    # Only compare to the one parent from which the allosome can originate
    parent_id = 1 if chrom == 'X' else 2  # compare: mother(1)<>X, father(2)<>Y
    if not can_compute_dng_DQ_male_allosome(variant, param, parent_id):
        return None

    # genotypes for proband and relevant parent
    gts = variant.genotypes
    gt_proband = set(gts[param['sample_index'][0]][:-1])
    gt_parent = set(gts[param['sample_index'][parent_id]][:-1])

    if gt_proband.isdisjoint(gt_parent):
        score = _max_score
    else:
        score = 0.0

    return score


def calculate_dng_DQ(variant, prior, param):
    """Calculate variant denovo score"""

    chrom = variant.CHROM.lstrip('chr')

    if param.get('female_proband') and chrom == 'Y':
        score = None
    elif is_in_regions(variant, param.get('par_regions'), chrom):
        # use autosomal model in PAR
        if variant.is_snp:
            score = calculate_dng_DQ_snv(variant, prior['auto']['snv'], param)
        elif variant.is_indel:
            score = calculate_dng_DQ_indel(variant, prior['auto']['indel'], param)
        elif variant.is_sv:
            score = calculate_dng_DQ_sv(variant, prior['auto']['indel'], param)
        else:
            score = None
    elif param.get('male_proband') and chrom in ('X', 'Y'):
        if variant.is_snp or variant.is_indel:
            score = calculate_dng_DQ_male_allosome(variant, chrom, param)
        else:
            score = None
    elif chrom == 'X':
        if variant.is_snp:
            score = calculate_dng_DQ_snv(variant, prior['chrX']['snv'], param)
        elif variant.is_indel:
            score = calculate_dng_DQ_indel(variant, prior['chrX']['indel'], param)
        elif variant.is_sv:
            score = calculate_dng_DQ_sv(variant, prior['chrX']['indel'], param)
        else:
            score = None
    else:
        if variant.is_snp:
            score = calculate_dng_DQ_snv(variant, prior['auto']['snv'], param)
        elif variant.is_indel:
            score = calculate_dng_DQ_indel(variant, prior['auto']['indel'], param)
        elif variant.is_sv:
            score = calculate_dng_DQ_sv(variant, prior['auto']['indel'], param)
        else:
            score = None

    return score


def rectify_pl_field(pl, max_pl=500):
    """Reformat PL matrix in case that samples have different ploidy"""

    if pl.min() < 0:
        for (i, p) in enumerate(pl):
            if p[2] < 0:
                p[2] = p[1]
                p[1] = max_pl
                pl[i] = p

    return pl


def calculate_dng_DQ_snv(variant, prior, param):
    """Calculate SNV denovo score"""

    if not can_compute_dng_DQ_snv(variant, param):
        return None

    if not is_denovo_candidate(variant, param, _mendel_conflict_name):
        return 0.0

    sample_index = param['sample_index']

    if param['is_select_call']:
        idx_trio = [_gts_gtsu_snv[variant.gt_bases[i]] for i in sample_index]
        idx_pick = _full_gt_idx_snv['/'.join(idx_trio)]
        if not prior['denovo_flag'][idx_pick]:
            return 0.0

    try:
        pl = get_PL_field(variant)
        alleles = tuple([variant.REF] + variant.ALT)
        alleles_idx = get_PL_indices_snv(alleles)
        C = map_sample_pl(pl[sample_index[0]], alleles_idx)
        M = map_sample_pl(pl[sample_index[1]], alleles_idx)
        D = map_sample_pl(pl[sample_index[2]], alleles_idx)
        P = np.outer(M, D).ravel()  # P
        PP = np.outer(P, C).ravel()  # F
        PP *= prior[variant.REF]
    except BaseException:
        return None

    if param['is_select_call']:
        ml_denovo = PP[idx_pick]
    else:
        ml_denovo = max(PP * prior['denovo_flag'])

    if ml_denovo > 0.0:
        pp_denovo = ml_denovo / sum(PP)
        score = prob2score(1.0 - pp_denovo)
    else:
        score = 0.0  # pragma: no cover

    return score


def indel_length(variant):
    """Get absolute length of indel"""

    ref_length = len(variant.REF)
    if len(variant.ALT) == 1:
        # fast track for most variants
        indel_len = abs(ref_length - len(variant.ALT[0]))
    else:
        # conservative estimate: smallest indel length
        indel_len = min([abs(ref_length - len(alt)) for alt in variant.ALT])

    return indel_len


def calculate_dng_DQ_indel(variant, prior, param):
    """Calculate indel denovo score"""

    if not can_compute_dng_DQ_indel(variant, param):
        return None

    if not is_denovo_candidate(variant, param, _mendel_conflict_name):
        return 0.0

    sample_index = param['sample_index']
    if param['is_select_call']:
        idx_trio = [_indel_gt_types[tuple(x == 0)]
                    for x in alleles_integers(variant)]
        idx_trio = [idx_trio[i] for i in sample_index]
        idx_pick = _full_gt_idx_indel['/'.join(idx_trio)]
        if not prior['denovo_flag'][idx_pick]:
            return 0.0

    try:
        pl = get_PL_field(variant)
        alleles = tuple(['R'] + ['V' for i in xrange(len(variant.ALT))])
        alleles_idx = get_PL_indices_indel(alleles)
        C = map_sample_pl_indel(pl[sample_index[0]], alleles_idx)
        M = map_sample_pl_indel(pl[sample_index[1]], alleles_idx)
        D = map_sample_pl_indel(pl[sample_index[2]], alleles_idx)
        P = np.outer(M, D).ravel()  # P
        PP = np.outer(P, C).ravel()  # F
    except:
        return None

    # indel prior
    indel_len = min(indel_length(variant), _max_indel_len)
    indel_type = 'del' if variant.is_deletion else 'ins'
    PP *= prior[indel_type][:, indel_len]

    if param['is_select_call']:
        ml_denovo = PP[idx_pick]
    else:
        ml_denovo = max(PP * prior['denovo_flag'])

    if ml_denovo > 0.0:
        pp_denovo = ml_denovo / sum(PP)
        score = prob2score(1.0 - pp_denovo)
    else:
        score = 0.0  # pragma: no cover

    return score


def calculate_dng_DQ_sv(variant, prior, param):
    """Calculate SV denovo score"""

    if not can_compute_dng_DQ_sv(variant, param):
        return None

    sample_index = param['sample_index']
    if param['is_select_call']:
        idx_trio = [_indel_gt_types[tuple(x == 0)]
                    for x in alleles_integers(variant)]
        idx_trio = [idx_trio[i] for i in sample_index]
        idx_pick = _full_gt_idx_indel['/'.join(idx_trio)]
        if not prior['denovo_flag'][idx_pick]:
            return 0.0

    try:
        pl = get_PL_field(variant)
        alleles = tuple(['R'] + ['V' for i in xrange(len(variant.ALT))])
        alleles_idx = get_PL_indices_indel(alleles)
        C = map_sample_pl_indel(pl[sample_index[0]], alleles_idx)
        M = map_sample_pl_indel(pl[sample_index[1]], alleles_idx)
        D = map_sample_pl_indel(pl[sample_index[2]], alleles_idx)
        P = np.outer(M, D).ravel()  # P
        PP = np.outer(P, C).ravel()  # F
    except:
        return None

    if param['is_select_call']:
        ml_denovo = PP[idx_pick]
    else:
        ml_denovo = max(PP * prior['denovo_flag'])

    if ml_denovo > 0.0:
        pp_denovo = ml_denovo / sum(PP)
        score = prob2score(1.0 - pp_denovo)
    else:
        score = 0.0  # pragma: no cover

    return score


def add_DQ_score_to_variant(variant, param, score,
                           decimals=_default_dq_decimals,
                           null_value='.'):
    """Add DQ field as a proband FORMAT field to VCF entry"""

    fields = str(variant).split()
    formats = variant.FORMAT

    sample_indexes = param['sample_index']
    sample_index_child = sample_indexes[0]
    null_value = '.'

    if _field_name in formats:
        # if variant already has a DQ format field
        # update the DQ score of the proband
        idx_dq = formats.index(_field_name)
        # only process the trio samples, leave additional samples untouched
        for i in sample_indexes:
            field_value, field_index = \
                select_sample_score(score, i, sample_index_child, decimals, null_value)
            sample_fields = fields[field_index].split(':')
            sample_fields[idx_dq] = field_value
            fields[field_index] = ':'.join(sample_fields)
    else:
        # if variant has no DQ format field yet
        # append the DQ score to the proband and 'missing' to all other samples
        fields[8] += ':' + _field_name
        for i in xrange(param['n_samples']):
            field_value, field_index = \
                select_sample_score(score, i, sample_index_child, decimals, null_value)
            fields[field_index] = fields[field_index] + ':' + field_value

    variant_string = '\t'.join(fields) + '\n'

    return variant_string


def select_sample_score(score, sample_index, sample_index_child, decimals, null_value='.'):
    """Assign and format the DQ score to each sample"""

    if sample_index == sample_index_child and score is not None:
        field_value = '%g' % round(score, decimals)
    else:
        field_value = null_value

    field_index = sample_index + 9

    return field_value, field_index


_models = {'dng': calculate_dng_DQ, 'dummy': calculate_dummy_DQ}


def can_compute_dng_DQ_common(variant, param):
    """Determines if a DQ score can be computed for any variant type"""

    if variant.FILTER is not None:  # not PASS
        return False
    if param['pass_sample_filters']:
        if _sample_filter_name not in variant.FORMAT:
            return False
        sample_index = param['sample_index']
        if param.get('male_proband'):
            chrom = variant.CHROM.lstrip('chr')
            if chrom == 'X':
                sample_index = param['sample_index_chrX']
            elif chrom == 'Y':
                sample_index = param['sample_index_chrY']
            else:
                pass
        ft = variant.format(_sample_filter_name)
        for i in sample_index:
            if ft[i] != 'PASS':
                return False
    pl = variant.format('PL', int)
    if pl is None:
        return False

    return True


def can_compute_dng_DQ_snv(variant, param):
    """Determines if a DQ score can be computed for the SNV"""

    if not variant.is_snp:
        return False  # pragma: no cover

    return can_compute_dng_DQ_common(variant, param)


def can_compute_dng_DQ_indel(variant, param):
    """Determines if a DQ score can be computed for the indel"""

    if not variant.is_indel:
        return False  # pragma: no cover

    return can_compute_dng_DQ_common(variant, param)


def can_compute_dng_DQ_sv(variant, param):
    """Determines if a DQ score can be computed for the SV"""

    if not variant.is_sv:
        return False  # pragma: no cover

    return can_compute_dng_DQ_common(variant, param)


def can_compute_dng_DQ_male_allosome(variant, param, parent_id):
    """Determines if a DQ score can be computed for an allosomal variant in a male proband"""

    proband_index = param['sample_index'][0]
    parent_index = param['sample_index'][parent_id]

    depths = variant.gt_depths
    min_depth = param['min_read_depth']
    min_proband_depth = min_depth / 2  # adjust threshold for haploid samples
    min_parent_depth = min_proband_depth if parent_id == 2 else min_depth
    if (depths is None) or (depths[proband_index] < min_proband_depth) or (depths[parent_index] < min_parent_depth):
        return None

    min_gqx = param['min_gqx']
    gqx = variant.format('GQX')
    if (gqx is None) or (gqx[proband_index] < min_gqx) or (gqx[parent_index] < min_gqx):
        return None

    return True


def import_bed_regions(path, fieldnames=('contig', 'start', 'end', 'annotation')):
    """Import the regions from a BED file"""

    try:
        regions = dict()
        with open(path) as f:
            reader = csv.DictReader(f, fieldnames, delimiter='\t')
            for entry in reader:
                contig = entry['contig'].lstrip('chr')
                coords = (int(entry['start']), int(entry['end']))
                regions.setdefault(contig, list()).append(coords)
        # convert values from list to tuple
        regions = {contig: tuple(coords) for contig, coords in regions.iteritems()}
    except (IOError, BaseException):
        regions = None

    return regions


def is_in_regions(variant, regions=None, chrom=None):
    """Check if the variant is located within the annotated regions"""

    if chrom is None:
        chrom = variant.CHROM.lstrip('chr')

    if regions and chrom in regions:
        for coords in regions[chrom]:
            if variant.POS >= coords[0] and variant.POS <= coords[1]:
                return True

    return False


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        'in_vcf',
        metavar='INPUT_VCF',
        help="""Multi-sample VCF input file""")

    parser.add_argument(
        'out_vcf',
        metavar='OUTPUT_VCF',
        nargs='?',
        default=None,
        help="""VCF file with DQ score""")

    parser.add_argument(
        '--proband',
        required=True,
        help="""Proband's sample ID""")

    parser.add_argument(
        '--mother',
        required=True,
        help="""Mother's sample ID""")

    parser.add_argument(
        '--father',
        required=True,
        help="""Father's sample ID""")

    parser.add_argument(
        '--model',
        required=False,
        choices=_models.keys(),
        default=_default_model_name,
        help=argparse.SUPPRESS)

    parser.add_argument(
        '--select-method',
        required=False,
        choices=_select_methods,
        default=_default_select_method,
        help=argparse.SUPPRESS)

    parser.add_argument(
        '--ignore-sample-filters',
        action='store_true',
        help="""Analyse also variants that have not set a 'FT' format field or fail this filter.""")

    parser.add_argument(
        '--ignore-pedigree-phasing',
        action='store_false',
        help=argparse.SUPPRESS)

    parser.add_argument(
        '--parallel',
        action='store_true',
        help="""Run computation in parallel, using GNU parallel"""
    )

    parser.add_argument(
        '--ncores',
        default=multiprocessing.cpu_count(),
        help="""Number of compute cores to use""")

    parser.add_argument(
        '--score-decimals',
        type=int,
        default=_default_dq_decimals,
        help="""Number of decimals to report for DQ score""")

    parser.add_argument(
        '--vcf-output',
        choices=('full', 'header', 'body'),
        default='full',
        help=argparse.SUPPRESS)

    parser.add_argument(
        '--par-bed-path',
        default='/illumina/sync/software/released/Isas/Genomes/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/PARv5.bed',
        help="""Path to the annotation BED files that defines the PAR (pseudoautosomal regions)""")

    parser.add_argument(
        '--min-read-depth',
        type=int,
        default=_min_read_depth,
        help="""Minimum required read depth for allosomal variants in a male proband""")

    parser.add_argument(
        '--min-gqx',
        type=int,
        default=_min_gqx,
        help="""Minimum required GQX score for allosomal variants in a male proband""")


    parser.add_argument('--version', action='version', version=_version)

    args = parser.parse_args()

    pedigree = {  # None if not specified
        'proband': args.proband,
        'mother': args.mother,
        'father': args.father
    }

    # import PAR coordinates
    par_regions = import_bed_regions(args.par_bed_path)

    param = {
        'model': args.model,
        'min_read_depth': args.min_read_depth,
        'pass_sample_filters': not args.ignore_sample_filters,
        'select_method': args.select_method,
        'score_decimals': args.score_decimals,
        'print_vcf_header': args.vcf_output in ('full', 'header'),
        'print_vcf_body': args.vcf_output in ('full', 'body'),
        'min_gqx': args.min_gqx,
        'filter_pedphase': args.ignore_pedigree_phasing,
        'par_regions': par_regions
    }

    if args.parallel:
        cmd_header = "{} {} --proband {} --mother {} --father {} --vcf-output {} {}".format(
            'python -E -s', __file__, args.proband, args.mother, args.father, 'header', copy.copy(args.in_vcf))
        denovo_cmd = "'{} {} --proband {} --mother {} --father {} --vcf-output {} -'".format(
            'python -E -s', __file__, args.proband, args.mother, args.father, 'body')
        cmd_unzip = "gunzip -cdf {}".format(
            args.in_vcf)
        cmd_body = "parallel -j {} --header '(\#.*\n)*' --keep-order --pipe --no-notice {}".format(
            args.ncores, denovo_cmd)

        try:
            output_file = select_output(args.out_vcf, 'w')
            p0 = sp.Popen(cmd_header, stdout=output_file, shell=True)
            p0.wait()
            output_file = select_output(args.out_vcf, 'a')
            p1 = sp.Popen(cmd_unzip, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
            p2 = sp.Popen(cmd_body, stdin=p1.stdout, stdout=output_file, shell=True)
            p2.wait()
            # handle premature closing of stdout
        except (KeyboardInterrupt, BaseException):
            p2.terminate()

    else:
        denovo(args.in_vcf, pedigree, param, args.out_vcf)
