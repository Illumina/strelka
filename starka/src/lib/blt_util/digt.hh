// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file

/// \author Chris Saunders
///
#ifndef __DIGT_HH
#define __DIGT_HH


namespace DIGT {
    enum index_t {
        AA,
        CC,
        GG,
        TT,
        AC,
        AG,
        AT,
        CG,
        CT,
        GT,
        SIZE
    };

    inline
    const char*
    label(const unsigned idx){
        switch(idx){
        case AA: return "AA";
        case CC: return "CC";
        case GG: return "GG";
        case TT: return "TT";
        case AC: return "AC";
        case AG: return "AG";
        case AT: return "AT";
        case CG: return "CG";
        case CT: return "CT";
        case GT: return "GT";
        default: return "XX";
        }
    }

    inline
    bool
    is_het(const unsigned idx){
        switch(idx){
        case AA:
        case CC:
        case GG:
        case TT: return false;
        default: return true;
        }
    }

    inline
    char
    hom_base(const unsigned idx){
        switch(idx){
        case AA: return 'A';
        case CC: return 'C';
        case GG: return 'G';
        case TT: return 'T';
        default: return 'X';
        }
    }

    inline
    double
    expect(const int base_id,
           const int gt) {

        static const unsigned N_BASE(4);

        static const double ex[SIZE][N_BASE] = {{ 1.0, 0.0, 0.0, 0.0},
                                                { 0.0, 1.0, 0.0, 0.0},
                                                { 0.0, 0.0, 1.0, 0.0},
                                                { 0.0, 0.0, 0.0, 1.0},
                                                { 0.5, 0.5, 0.0, 0.0},
                                                { 0.5, 0.0, 0.5, 0.0},
                                                { 0.5, 0.0, 0.0, 0.5},
                                                { 0.0, 0.5, 0.5, 0.0},
                                                { 0.0, 0.5, 0.0, 0.5},
                                                { 0.0, 0.0, 0.5, 0.5}};

        return ex[gt][base_id];
    }

    /// coded form of expect function:
    /// 0 is 0.0
    /// 1 is 0.5
    /// 2 is 1.0
    inline
    unsigned
    expect2(const int base_id,
            const int gt) {

        static const unsigned N_BASE(4);

        static const unsigned ex[SIZE][N_BASE] = {{ 2, 0, 0, 0},
                                                  { 0, 2, 0, 0},
                                                  { 0, 0, 2, 0},
                                                  { 0, 0, 0, 2},
                                                  { 1, 1, 0, 0},
                                                  { 1, 0, 1, 0},
                                                  { 1, 0, 0, 1},
                                                  { 0, 1, 1, 0},
                                                  { 0, 1, 0, 1},
                                                  { 0, 0, 1, 1}};

        return ex[gt][base_id];
    }

    /// expect2_bias is a copy of the expect2 function for biased het
    /// allele calculations
    ///
    /// 0 is 0.0
    /// 1 is het-ratio  (defined by client-code)
    /// 2 is (1.-het-ratio)
    /// 3 is 1.0
    ///
    /// for consistency, state 1 is always applied to the lower allele
    /// value, and state 2 to the higher.
    ///
    inline
    unsigned
    expect2_bias(const int base_id,
                 const int gt) {

        static const unsigned N_BASE(4);

        static const unsigned ex[SIZE][N_BASE] = {{ 3, 0, 0, 0},
                                                  { 0, 3, 0, 0},
                                                  { 0, 0, 3, 0},
                                                  { 0, 0, 0, 3},
                                                  { 1, 2, 0, 0},
                                                  { 1, 0, 2, 0},
                                                  { 1, 0, 0, 2},
                                                  { 0, 1, 2, 0},
                                                  { 0, 1, 0, 2},
                                                  { 0, 0, 1, 2}};

        return ex[gt][base_id];
    }

}


#endif
