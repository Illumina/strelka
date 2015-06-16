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
/*
 * assembler.hh
 *
 *  Test for assembler.
 *
 *  Created on: Aug 10, 2014
 *  Author: Morten Kallberg
 */


#pragma once

#include <assembly/assembly_common/AssemblyReadInfo.hh>
#include <assembly/assembly_common/BamAddOns.hh>
#include <assembly/assembly_common/Contig.hh>
#include "assembly/kmer/DfsAssemblerOptions.hh"

#include "assembly/graph/GraphBase.hh"

#include "boost/foreach.hpp"
#include "boost/unordered_map.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/utility.hpp"

#include <cassert>
#include <vector>
#include <fstream>
#include <iostream>



