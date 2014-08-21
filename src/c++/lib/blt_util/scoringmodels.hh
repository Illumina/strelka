/*
 * scoringmodels.hh
 *
 *  Created on: Aug 20, 2014
 *      Author: mkallberg
 */

#ifndef SCORINGMODELS_HH_
#define SCORINGMODELS_HH_

class scoring_models{
public:
   static scoring_models* Instance();
//   bool openLogFile(std::string logFile);
   void load_models(const std::string& model_file);
//   bool closeLogFile();

private:
   scoring_models(){};  // Private so that it can  not be called
   scoring_models(scoring_models const&){};             // copy constructor is private
   scoring_models& operator=(scoring_models const&);  // assignment operator is private
   static scoring_models* m_pInstance;
};

#endif /* SCORINGMODELS_HH_ */
