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
 * LogisticRegressionModel.hh
 *
 *  Created on: Jun 23, 2015
 *      Author: mkallberg
 */

#include <calibration/SerializedModel.hh>

#ifndef C___LIB_CALIBRATION_LOGISTICREGRESSIONMODEL_HH_
#define C___LIB_CALIBRATION_LOGISTICREGRESSIONMODEL_HH_

struct LogisticRegressionModel : public serialized_model
{

    LogisticRegressionModel();
    virtual ~LogisticRegressionModel();

    void Deserialize( const Json::Value& root);

};

#endif /* C___LIB_CALIBRATION_LOGISTICREGRESSIONMODEL_HH_ */
