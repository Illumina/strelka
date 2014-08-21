/*
 * scoringmodels.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: mkallberg
 */

#include <stdlib.h>     /* atof */
#include <iostream>
#include <sstream>
#include <cassert>
#include <string>
//#include <fstream>
//#include <iterator>
#include <boost/property_tree/ptree.hpp>
//#include <boost/algorithm/string/split.hpp>
//#include <boost/algorithm/string/classification.hpp>
#include "blt_util/scoringmodels.hh"


#define DEBUG_SCORINGMODELS
#ifdef DEBUG_SCORINGMODELS
    #include "blt_util/log.hh"
#endif

// Global static pointer used to ensure a single instance of the class.
scoring_models* scoring_models::m_pInstance = nullptr;

/** This function is called to create an instance of the class.
    Calling the constructor publicly is not allowed. The constructor
    is private and is only called by this Instance function.
*/
scoring_models* scoring_models::Instance()
{
   if (!m_pInstance)   // Only allow one instance of class to be generated.
      m_pInstance = new scoring_models;
   return m_pInstance;
}

void scoring_models::load_models(const std::string& model_file){
//    if (model_file.length()>2)
    log_os << "loading file " << model_file << "\n";
}
