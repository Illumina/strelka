
#include <assembly_common/VcfParser.hh>

VcfRecord::VcfRecord(std::string c, int p, int l, std::string i, std::string r, std::string a,
                     int q, std::string fi, std::string in, std::string fo, std::string g)
    : chr(c), pos(p), len(l), id(i), ref(r), alt(a), qual(q), filter(fi),
  /*info(in),*/ format(fo), gt(g)
{
    parseInfoStringAndAssignLength(in,false);
}

void VcfRecord::clear() {
    infoFields.clear();
    chr="";
    pos=0;
    len=0;
    id="";
    ref="";
    alt="";
    qual=0;
    filter="";
    format="";
    gt="";
}

bool VcfRecord::isDeletion() const  {
    vcfInfoMapT::const_iterator svtypeCiter = infoFields.find("SVTYPE");
    if (svtypeCiter == infoFields.end()) {
        return false;
    }
    if (svtypeCiter->second != "DEL") {
        return false;
    }
    return true;
}

bool VcfRecord::isOpenEnded() const  {
    vcfInfoMapT::const_iterator upStreamCiter = infoFields.find("UPSTREAM");
    if (upStreamCiter == infoFields.end()) {
        return false;
    }
    vcfInfoMapT::const_iterator downStreamCiter = infoFields.find("DOWNSTREAM");
    if (downStreamCiter == infoFields.end()) {
        return false;
    }
    if (upStreamCiter->second != "." && downStreamCiter->second != ".") {
        return false;
    }
    return true;
}

bool VcfRecord::isLeftOpenEnded() const  {
    vcfInfoMapT::const_iterator upStreamCiter = infoFields.find("UPSTREAM");
    if (upStreamCiter == infoFields.end()) {
        // no upstream sequence is a no
        return false;
    }
    if (upStreamCiter->second != ".") {
        return false;
    }
    return true;
}

bool VcfRecord::isRightOpenEnded() const  {
    vcfInfoMapT::const_iterator downStreamCiter = infoFields.find("DOWNSTREAM");
    if (downStreamCiter == infoFields.end()) {
        return false;
    }
    if (downStreamCiter->second != ".") {
        return false;
    }
    return true;
}

bool VcfRecord::isTranslocation() const {
    std::string bndOpenEndedTag         = "BND:OPENEND";
    std::string bndOpenEndedTranslocTag = "BND:OPENEND:TRN";
    vcfInfoMapT::const_iterator citer = infoFields.find("SVTYPE");
    if (citer == infoFields.end()) {
        // if no svtype is defined, we cannot say anything
        return false;
    }
    if (citer->second == "BND") {
        size_t foundOpenEnded = alt.find(bndOpenEndedTag);
        if (foundOpenEnded == std::string::npos) {
            // BND event and not open-ended, must be translocation
            return true;
        } else {
            size_t foundOpenEndedTransloc = alt.find(bndOpenEndedTranslocTag);
            if (foundOpenEndedTransloc == std::string::npos) {
                // open-ended BND but not translocation
                return false;
            } else {
                // open-ended BND and translocation
                return true;
            }
        }
    }
    return false;
}

void VcfRecord::updateVcfFilter(const std::string& newFilterTag) {
    if ((filter == "PASS") ||
        (filter == ".") ||
        filter.empty()) {
        filter = newFilterTag;
    } else  {
        filter += ";" + newFilterTag;
    }
}

void VcfRecord::updateVcfGT (const std::string& newGT, const std::string& newGQ) {
    gt = newGT + ":" + newGQ;
}

bool VcfRecord::parseInfoStringAndAssignLength(const std::string& infoString, const bool reverse) {
    bool svLenParsed(false);
    boost::tokenizer< boost::char_separator<char> > infoTokens(infoString, boost::char_separator<char>(";"));
    BOOST_FOREACH (const std::string& tok, infoTokens) {
        //std::cout << "token=" << t << std::endl;
        std::vector<std::string> fields;
        boost::split(fields,tok,boost::is_any_of("="));
        if (fields.size()==1) {
            fields.push_back("N/A");
        }
        assert(fields.size()==2);
        infoFields[fields[0]]=fields[1];
        if (fields[0] == "SVLEN") {
            //std::cout << "SVLEN found in: " << t << std::endl;
            try {
                len = abs(boost::lexical_cast<int>(fields[1]));
                svLenParsed = true;
                if (reverse) {
                    len *= -1;
                }
            } catch (boost::bad_lexical_cast&) {
                std::cout << "Error while trying to convert SVLEN info entry to numerical value: " << fields[1] << std::endl;
            }
            //std::cout << "Extracted event length=" << rec.len << std::endl;
            //assert(rec.len>0);
        }
    }
    return svLenParsed;
}

std::ostream& operator<<(std::ostream& os, const VcfRecord& r) {
    os << r.chr << "\t" << r.pos << "\t" << r.id << "\t" << r.ref << "\t"
       << r.alt << "\t";
    // write quality value only if positive, otherwise set placeholder
    if (r.qual < 0) {
        os << '.';
    } else {
        os << r.qual;
    }
    os << "\t" << r.filter << "\t";
    // print info fields
    for (vcfInfoMapT::const_iterator ct=r.infoFields.begin(); ct!=r.infoFields.end(); ++ct) {
        if (ct!=r.infoFields.begin()) {
            os << ";";
        }
        if (ct->second=="N/A") {
            os << ct->first;
        } else {
            os << ct->first << "=" << ct->second;
        }
    }
    //<< r.info
    return os << "\t" << r.format << "\t" << r.gt;
}

