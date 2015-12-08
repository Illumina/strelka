#include "denovo_snv_call_vcf.hh"
#include <algorithm>


void
denovo_snv_call::get_alt(){
		for (unsigned i(0);i<SampleGts.size();++i){
			std::array<unsigned,2> gt = {0,0};
			for (unsigned chrom=0;chrom<2;chrom++){
				if (SampleGts[i][chrom] != this->ref_gt){

					//make sure all recorded alt alles are recorded
					if(std::find(alts.begin(), alts.end(), SampleGts[i][chrom]) == alts.end()){
						alts.push_back(SampleGts[i][chrom]);
						if (alts.size()>1)
							alt_str += ",";
						alt_str += id_to_base(alts[alts.size()-1]);
					}

					//genotype sample according to ALT map
					for (unsigned t=0;t<alts.size();t++)
						if (alts[t]==SampleGts[i][chrom])
							gt[chrom] = (t+1);
				}
			}
			if (gt[0]>gt[1]){
				unsigned temp = gt[1];
				gt[1] = gt[0];
				gt[0] = temp;
			}
			gts_chrom.push_back(gt);
 		}

		if (alts.size()==0)
			alt_str = ".";
}
