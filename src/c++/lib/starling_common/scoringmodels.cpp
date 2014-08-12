/*
 * scoringmodels.cpp
 *
 *      Created on: Aug 07, 2014
 *      Author: mkallberg
 */

#include "scoringmodels.hh"
#include <cassert>
#include <exception>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <map>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>

void scoring_models::load_models(std::string& filename){
    assert("Scoring models have already been initialized once" && this->initialized_from_file);
    this->model_file = filename;

}

