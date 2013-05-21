// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file

/// \author Morten Kallberg
///

#include "blt_common/ranksum.hh"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <iostream>
#include <map>


double
get_z_score(int n1, int n2, double w1){
    double mean = n1*(n1+n2+1)/2.0;
    double var  = sqrt(n1*n2*(n1+n2+1)/12.0);
    if (static_cast<int>(var)==0){
        return 0.0;
    }
    double z = (w1-mean)/var;
    return z;
}

// return the U statistic
double
ranksum::get_u_stat()
{
    //	cout << "doing U stat" << endl;
    vector<int> myvector = this->getSpace();
    int current_rank = 1;
    R1 = 0.0;
    R2 = 0.0;
    N1 = 0;
    N2 = 0;
    //loop over all observations
    for (std::vector<int>::iterator it=myvector.begin(); it!=myvector.end(); ++it){
        int key = *it;

        // get the observation counts for reference and alt
        int obs1= this->get_obs_count(true,key);
        int obs2= this->get_obs_count(false,key);
        int obs = obs1 + obs2;
        double rank_weight = (2*current_rank + obs -1)/2.0;
        R1 += rank_weight*obs1;
        R2 += rank_weight*obs2;
        N1 += obs1;
        N2 += obs2;
        current_rank += obs;
    //		cout << key << "\t{"<< obs1 << ", " << obs2 << "}" << " sum: " << obs ;
    //		cout << " weight " << rank_weight << " U1: " << U1 << " U2: " << U2 << endl;
    //		cout << "updated rank " << current_rank << endl;
	}

	//return the z-score for the smaller of U1 and U2 assuming a gaussian approx.
    double z;
    if (R1>R2){
        z = get_z_score(N2,N1,R2);
    }
    else{
        z = get_z_score(N1,N2,R1);
    }

    return z;
}

void
ranksum::add_observation(bool is_ref, int obs){
    space[obs] =1;
    if (is_ref){
        l1[obs]++;
    //		cout << "ref case" << endl;
    }
    else{
        l2[obs]++;
    //		cout << "alt case" << endl;
    //		cout << "ref_base: " << ref_base << " alt_base: " << base <<   endl;
    }
}


// output specification

ostream&
operator<< (ostream &out, map<int, int> &l)
{
    out << "Elements in l: " << endl;
    for (map<int, int>::iterator it = l.begin();it != l.end(); ++it)
        out << "\t" <<(*it).first << " => " << (*it).second <<  endl;
    return out;
}

ostream&
operator<< (ostream &out, ranksum &r)
{
    double z = r.get_u_stat();
    out << endl << "My reference base: " << r.get_refbase() << endl;
    out << "elements in space: " << endl;
    out << "[";
    vector<int> myvector = r.getSpace();
    for (std::vector<int>::iterator it=myvector.begin(); it!=myvector.end(); ++it)
        std::cout << ' ' << *it;
    out << "]" << endl;
    out << r.l1;
    out << r.l2;
    out << "N1: " << r.N1 << "\tR1: " << r.R1 <<  endl;
    out << "N2: " << r.N2 << "\tR2: " << r.R2 <<  endl;
    out << "Z-score: " << z << endl;
    return out;
}



// for testing
//int main()
//{
//  ranksum r('A');
//  r.add_observation('A',44);
//  r.add_observation('A',45); //check tie both
//  r.add_observation('N',45);
//  r.add_observation('N',50);
//  r.add_observation('A',52);
//  r.add_observation('A',53);
//  r.add_observation('A',56);
//  r.add_observation('A',58); //check tie single
//  r.add_observation('A',58);
//  r.add_observation('N',61);
//  r.add_observation('N',63);
//  r.add_observation('A',65);
//  r.add_observation('N',75);
//  r.add_observation('A',79);
//  r.add_observation('N',85);
//  r.add_observation('N',93);
//  cout << r << '\n';
//  cout << r.get_u_stat() << endl;
////  cout << get_z_score(00,12,0.0) << endl;
//  return 0;
//}

