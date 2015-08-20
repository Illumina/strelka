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
 * LogisticRegressionModel.cpp
 *
 *  Created on: Jun 23, 2015
 *      Author: mkallberg
 */

#include <LogisticRegressionModel.hh>

LogisticRegressionModel::LogisticRegressionModel()
{
    // TODO Auto-generated constructor stub

}

LogisticRegressionModel::~LogisticRegressionModel()
{
    // TODO Auto-generated destructor stub
}

void LogisticRegressionModel::Deserialize( const Json::Value& root)
{
    serialized_model::Deserialize(root);
}
