
#include <assembly_common/VcfParser.hh>

bool VcfParser::buildVcfRecordFromString(const std::string& vcfLine, VcfRecord& rec) {
    rec.clear();
    std::vector<std::string> vcfTokens;
    _tokenizeLine(vcfLine,vcfTokens);

    switch (vcfTokens.size()) {
    case 8:
        rec.chr=vcfTokens[0];
        rec.id=vcfTokens[2];
        rec.ref=vcfTokens[3];
        rec.alt=vcfTokens[4];
        rec.filter=vcfTokens[6];
        rec.format="";
        rec.gt="";
        if (!_convertPosAndQual(vcfTokens[1],vcfTokens[5],rec)) return false;
        break;
    case 10:
        rec.chr=vcfTokens[0];
        rec.id=vcfTokens[2];
        rec.ref=vcfTokens[3];
        rec.alt=vcfTokens[4];
        rec.filter=vcfTokens[6];
        rec.format=vcfTokens[8];
        rec.gt=vcfTokens[9];
        if (!_convertPosAndQual(vcfTokens[1],vcfTokens[5],rec)) return false;
        break;
    default:
        std::cout << "Cannot parse line " << vcfLine << " with " << vcfTokens.size() << " tokens." << std::endl;
        return false;
    }
    bool reverse(false);
    size_t idx1 = rec.id.find("bnd_U");
    size_t idx2 = rec.id.find("bnd_X");
    if (idx1 != std::string::npos || idx2 != std::string::npos) {
        reverse=true;
    }

    if (!rec.parseInfoStringAndAssignLength(vcfTokens[7],reverse)) {
        //std::cerr << "WARNING : cannot parse length from info field : " << vcfTokens[7] << std::endl;
    }
    return true;
}

bool VcfParser::_convertPosAndQual(const std::string& posStr, const std::string& qualStr, VcfRecord& rec) {
    try {
        rec.pos  = boost::lexical_cast<int>(posStr);
        if (qualStr != ".") {
            rec.qual = boost::lexical_cast<float>(qualStr);
        } else {
            rec.qual = -1;
        }
    } catch (boost::bad_lexical_cast& e) {
        std::cout << "Cannot cast pos/quality." << std::endl;
        std::cout << "pos  = " << posStr << std::endl;
        std::cout << "qual = " << qualStr << std::endl;
        std::cout << e.what() << std::endl;
        return false;
    }
    return true;
}


void VcfParser::_tokenizeLine(const std::string& vcfLine, std::vector<std::string>& vcfTokens) {
    boost::tokenizer< boost::char_separator<char> > vcfLineTokens(vcfLine,boost::char_separator<char>("\t"));
    BOOST_FOREACH(const std::string& t, vcfLineTokens) {
        vcfTokens.push_back(t);
    }
}

