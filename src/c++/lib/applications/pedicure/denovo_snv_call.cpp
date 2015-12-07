#include "denovo_snv_call_vcf.hh"
#include <algorithm>
#include "blt_util/log.hh"

void
denovo_snv_call::get_alt(){
		for (unsigned i(0);i<SampleGts.size();++i){
			for (unsigned chrom=0;chrom<2;chrom++){
				if (SampleGts[i][chrom] != this->ref_gt &&
					std::find(alts.begin(), alts.end(), SampleGts[i][chrom]) == alts.end()){
					alts.push_back(SampleGts[i][chrom]);
				}
			}
		}

		if (alts.size()>0)
			alt_str = id_to_base(alts[0]);
		for (unsigned i=1; i<alts.size(); i++){
			alt_str += ",";
			alt_str += id_to_base(alts[i]);
		}

}
