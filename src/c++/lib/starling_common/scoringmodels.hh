/*
 * scoringmodels.hh
 *
 *  Created on: Aug 07, 2014
 *      Author: mkallberg
 */

#ifndef SCORINGMODELS_HH_
#define SCORINGMODELS_HH_

class scoring_models {
   public:
       static scoring_models& getInstance()
       {
           static scoring_models instance;  // Guaranteed to be destroyed. Instantiated on first use.
           return instance;
       }

       void load_models(std::string& filename);
   private:
       scoring_models() {};                   // Constructor? (the {} brackets) are needed here.
       // Dont forget to declare these two. You want to make sure they
       // are unaccessable otherwise you may accidently get copies of
       // your singleton appearing.
       scoring_models(scoring_models const&);              // Don't Implement
       void operator=(scoring_models const&);              // Don't implement
       virtual ~scoring_models();

       bool initialized_from_file=false;                   // Did we load anything from models.json
       const starling_options& opt
};

#endif /* SCORINGMODELS_HH_ */
