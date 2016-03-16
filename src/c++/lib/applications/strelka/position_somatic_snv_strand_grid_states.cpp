//// -*- mode: c++; indent-tabs-mode: nil; -*-
////
//// Strelka - Small Variant Caller
//// Copyright (c) 2009-2016 Illumina, Inc.
////
//// This program is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// at your option) any later version.
////
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//// GNU General Public License for more details.
////
//// You should have received a copy of the GNU General Public License
//// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////
////
//
///// \author Chris Saunders
/////
//
//#include "position_somatic_snv_strand_grid_states.hh"
//#include "blt_util/digt.hh"
//
//#include <iostream>
//
//namespace DDIGT_SGRID
//{
//
//void
//write_state(const DDIGT_SGRID::index_t dgt,
//            const unsigned /* ref_gt */,
//            std::ostream& os)
//{
//    unsigned normal_gt;
//    unsigned tumor_gt;
//    DDIGT_SGRID::get_digt_grid_states(dgt,normal_gt,tumor_gt);
//
//    os << DIGT_SIMPLE::label(normal_gt);
//    os << "->";
//    os << DIGT_SIMPLE::label(tumor_gt);
//}
//
//void
//write_alt_alleles(unsigned alt_gt,
//                  std::ostream& os)
//{
//    os << id_to_base(alt_gt);
//}
//
//}
