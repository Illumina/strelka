#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2017 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

from collections import *
#import Counter
# from dynamicModelPlot import *
from dynamicModelUtil import *
import pickle
from pysam import Samfile, Fastafile
from scipy.stats import binom
import os

####
# Intended to emulate the behavior of the Python 2.7.3 module Counter from the packages collections 
# for backwards compatibility with Python 2.4.x installs on Centos 5 
####

#class Counter(dict):
#    def __getitem__(self, item):
#        if item not in self.iterkeys():
#            self[item] = 0
#        return self[item]
#    
#    def __setitem__(self, item):
#        if item not in self.iterkeys():
#            self[item] = 0

class job:
    def __init__(self,opt):
        self.chr            = opt.genomicRegion.split(":")[0]           # chromosomal coordinates to sample
        self.start          = int(opt.genomicRegion.split(":")[1].split('-')[0])
        self.end            = int(opt.genomicRegion.split(":")[1].split('-')[1])
        self.bam            = opt.inputBam          # path of input bam
        self.outputDir      = opt.output            # directory to output temp files
        self.reference      = opt.reference         # Path of reference fasta
        self.countEveryBase = opt.countEveryBase    # count on base-by-base case, rather than each homopolymer
        self.ignoreStrand   = opt.ignoreStrand      # do not consider if a read is forward or reverse strand
        self.naiveProject   = opt.naiveProject      # do fast but simple counting, OFF as default
        self.countUnspanned = opt.countUnspanned    # Consider reads that do not span the full length of the hpol as evidence, OFF by default
        self.maxMotifLen    = opt.maxMotifLen
        self.refOffset      = 500                   # additional up and downstream BPs to buffer reference
        self.verbose        = 1
        self.bufferRef()
        self.setDatastructures(opt)

    def bufferRef(self):
        # load reference sequence
        refStart = self.start - self.refOffset
        if refStart < 0:
            self.refStart = 0
        refEnd = self.end + self.refOffset
        fa = Fastafile(self.reference)
        print "Buffering from fasta: ", self.chr, ":", refStart, "-", refEnd
        self.refInterval = fa.fetch(reference=self.chr, start=refStart, end=refEnd)
        self.buffer = ReferenceBuffer(self.refInterval, self.numericChr(), 
                                      start = refStart, end = refEnd, 
                                      maxMotifLen = self.maxMotifLen)
        self.buffer.buffer()

    def clearEvents(self):
        self.events      = {'ins':{},'del':{}}

    def setDatastructures(self,opt):
        self.calls       = []
        self.eventCount  = {}
        for x in xrange(1, self.maxMotifLen + 1):
            self.eventCount[x] = {'tot':Counter(), 'ins':Counter(), 'del':Counter()}
        self.caseCounter = caseCounter(self.eventCount, opt)
        # for base in ['A','C','T','G','N']: self.eventCount[base] = Counter()
        self.clearEvents()

    def run(self):
        readRecords(self)

    def pickleMe(self):
         f = open(os.path.join(self.outputDir,str(self.chr) +"_"+str(self.start) + ".pickle"),'w')
         pickle.dump(self.caseCounter, f)
         f.close()

    def numericChr(self):
        return int(self.chr.replace('chr',''))

    def __repr__(self):
        return str(self)
    def __str__(self):
        return self.chr + ":" + str(self.start) + "-" + str(self.end) + " - " + str(self.eventCount)  
    
    def printDic(self):
        for type in ['del']:
            for case in self.events[type].keys():
                print case
                print self.events[type][case]
    
  
# cigar str numbering def. by pysam 
#BAM_CMATCH=0, BAM_CINS=1,
#DEF BAM_CDEL       = 2
#DEF BAM_CREF_SKIP  = 3
#DEF BAM_CSOFT_CLIP = 4
#DEF BAM_CHARD_CLIP = 5
#DEF BAM_CPAD       = 6
def readRecords(job,seqContext='All'):
    '''
    Generates the list of indel instances from the which the error rate is estimated
    :param job:
    :param seqContext:
    '''
#    mychr, start, end, countEveryBase, ignoreStrand, naiveProject, countUnspanned = job
    print "Collecting data "
    f  = Samfile(job.bam)
    fa = Fastafile(job.reference)
    i = 0
    for r in f.fetch(job.chr, job.start, job.end):
        if r.mapq > 20 and r.is_paired:# only concerned with reads of high quality so we don't count alignment errors
            myEvent  = 'm'       # treat clipped reads as matches
            myPos    = r.pos
            readPos  = 0
            # for a deletion, the deleted sequence must be a multiple of the motif length and be
            # composed only of the motif sequence. For an insertion, the inserted sequence must be a
            # multiple of the motif length and have the same sequence
            isRepSeq = False
            if len(r.cigar)>1:
                for event in r.cigar:
                    if event[0] < 1 or event[0] > 2:
                        myPos += event[1]
                        readPos += event[1]
                        continue
                    tractLength, motif, motifLength = job.buffer.repMap[myPos]
                    eventSeq   = ''
                    numRepeats = 0
                    isCompleteMotifs = True
                    myEvent    = ''
                    if event[0] == 1: 
                        myEvent = 'ins'
                        eventSeq = r.seq[readPos:readPos + event[1]]
                        isCompleteMotifs = event[1] % motifLength == 0
                        numRepeats = event[1] / motifLength
                    elif event[0] == 2:
                        myEvent  = 'del'
                        eventSeq = fa.fetch(reference = job.chr, start = myPos,
                                            end=(myPos + event[1])).upper()
                        isCompleteMotifs = event[1] % motifLength == 0
                        numRepeats = event[1] / motifLength

                    if isCompleteMotifs:
                        matchSeq = "".join([motif for _ in xrange(numRepeats)])
                        if eventSeq == matchSeq:
                            isRepSeq = True
                    # if myPos > 10008915 and myPos < 10008925:
                    #     print myEvent, tractLength, motif, eventSeq, myPos, isRepSeq
                    #     print myEvent, "event", eventSeq, "match", matchSeq
                        # print r.seq
                    break
                    # else:
                    #     myPos += event[1]
                    #     readPos += event[1]
                # add in information
                if job.buffer.repMap.has_key(myPos):
                    # if a reference STR has occurred at this position
                    strCase, repeat, motifLength = job.buffer.repMap[myPos]
                    if r.is_reverse and not job.ignoreStrand:
                        repeat = RC(repeat)
                    if not myEvent == 'm':
                        # newSeq = ""
                        # if myEvent == "ins":
                        #     newSeq = r.seq[readPos:readPos + event[1]]
                        # print isRepSeq, myEvent, event[1], myPos, strCase, repeat, motifLength, newSeq
                        if isRepSeq:
                            if not job.events[myEvent].has_key(myPos):
                                # if this event at this position has never been seen before, set it
                                # count, event length, STR length, motif length, seq repeated
                                #print [1,event[1],hpoleCase,repeat]
                                if seqContext=='All' or seqContext==repeat.upper():
                                    myLength = event[1]
                                    job.eventCount[motifLength][myEvent][strCase] +=1
                                    job.events[myEvent][myPos] = [1, myLength, strCase, repeat, motifLength]
                            else:
                                if seqContext=='All' or seqContext==repeat.upper():
                                    job.eventCount[motifLength][myEvent][strCase] += 1
                                    job.events[myEvent][myPos][0] += 1
                        else:
                            if not job.events[myEvent].has_key(myPos):
                                if seqContext=='All' or seqContext==repeat.upper():
                                    myLength = event[1]
                                    job.eventCount[1][myEvent][1] +=1
                                    job.events[myEvent][myPos] = [1, myLength, 1, repeat[-1], 1]
                            else:
                                if seqContext=='All' or seqContext==repeat.upper():
                                    job.eventCount[1][myEvent][1] += 1
                                    job.events[myEvent][myPos][0] += 1

            if job.countUnspanned: # i.e. skip logic to check whether alignment extends past the homopolymer
                # we can either count only reference bases actually covered or we can project based on nominal length of the read
                # first, set according to the actual alignment
                countCoverStart = r.pos + 1  # r.pos+1 so we know we have at least one base before the hpol; r.alen handles trimming and indel adjustments
                countCoverEnd   = r.pos + r.alen
                if job.naiveProject:
                    # if user so requests, reset to projection based on read length
                    countCoverStart = r.pos
                    countCoverEnd   = r.pos + 100 # assume a 100 bp readlength

                for ii in xrange(countCoverStart, countCoverEnd):
                    if job.buffer.hpolStartsMap.has_key(ii) or (job.countEveryBase and job.buffer.hpolMap.has_key(ii)):
                        hpoleCase, repeat = job.buffer.hpolMap[ii] # note hpol has same data as hpolStarts where hpolStarts is defined
                        if r.is_reverse and not job.ignoreStrand:
                            repeat = RC(repeat)
                        if not job.eventCount.has_key(repeat): job.eventCount[motifLength][repeat] = Counter()
                        job.eventCount[motifLength]['tot'][strCase]  +=1
                        job.eventCount[motifLength][repeat][strCase] +=1

            else:
                refOffset    = job.start - job.refOffset
                alignedPairs = r.aligned_pairs
                numpairs     = len(alignedPairs)

                prevMatch = False
                for ii in xrange(0, numpairs):
                    (x, y) = alignedPairs[ii]
                    if x != None:
                        readbase = r.query_sequence[x]
                    else:
                        readbase = '-'
                    if y != None:
                        refbase = job.refInterval[y - refOffset].upper()
                    else:
                        refbase = '-'


                    rpos = None
                    qpos = None

                    ###### count only events we can have confidence in:
                    ###### require that the alignment extend beyond a homopolymer by at least one base in each direction
                    ###### require that the base after the homopolymer (in the direction of sequencing) match

                    # if current position is the start of an STR
                    if job.buffer.repStartsMap.has_key(y) or (job.countEveryBase and job.buffer.repMap.has_key(y)):
                        strCase, repeat, motifLength = job.buffer.repMap[y] # note hpol[y] == hpolStarts[y] if the latter is defined

                        # make sure we align past the end of the hpol regions
                        if y + strCase >= r.pos + r.alen:
                            continue

                        ok = False # do not use this alignment position unless we satisfy conditions below

                        # if reverse strand alignment:
                        if r.is_reverse:
                            # require previous position to match (and hence also that the alignment goes beyond the homopolymer
                            if prevMatch:
                                ### we can use it!
                                if not job.ignoreStrand: 
                                    # print repeat, RC(repeat)
                                    repeat = RC(repeat)
                                ok = True
                        else: # forward strand alignment
                            # require this to not be the first aligned base of the read
                            if ii > 0:
                                rpos = None
                                qpos = None
                                # find alignment of first base past homopolymer, or an interior indel
                                for iii in xrange(ii + 1, numpairs):
                                    rpos = alignedPairs[iii][1]
                                    qpos = alignedPairs[iii][0]
                                    if (not rpos) or (not qpos) or rpos == y + strCase:
                                        break
                                # require that we be at the right position and it be an alignment of bases and they match
                                if rpos == y + strCase and qpos and job.refInterval[rpos - refOffset].upper() == r.query_sequence[qpos]:
                                    ### we can use it!
                                    ok = True

                        if ok:
                            # print myEvent, isRepSeq, motifLength, repeat, strCase
                            if myEvent != "m":
                                if not job.eventCount[motifLength].has_key(repeat):
                                    if isRepSeq:
                                        job.eventCount[motifLength][repeat] = Counter()
                                    else:
                                        job.eventCount[1][repeat[-1]] = Counter()
                                if isRepSeq:
                                    job.eventCount[motifLength]['tot'][strCase]  +=1
                                    job.eventCount[motifLength][repeat][strCase] +=1
                                else:
                                    job.eventCount[1]['tot'][1]  +=1
                                    job.eventCount[1][repeat[-1]][1] +=1
                            else:
                                if not job.eventCount[motifLength].has_key(repeat):
                                    job.eventCount[motifLength][repeat] = Counter()
                                job.eventCount[motifLength]['tot'][strCase]  +=1
                                job.eventCount[motifLength][repeat][strCase] +=1

                    if readbase == refbase:
                        prevMatch = True
                    else:
                        prevMatch = False

            i += 1
            if i % 10000 == 0:
                print "validating " + str(i)
                validateRecords(job) #populate the TP/FP counts in the events dict
    print "validating " + str(i)
    validateRecords(job)


def validateRecords(job, mode='BinomTest', freq=0.07, countCut=3, minBinomPval=1e-4):
# def validateRecords(job, mode='freq', freq=0.07, countCut=3, minBinomPval=1e-4):
    '''
    Determine the FP and TP status of a collection of events
    :param calls: List of call events
    :param job: The
    :param mode: eval mode, either count, freq, binomial test, or PG
    :param freq:   
    '''
    if mode=="PG":
        # import PGtools stuff, only for development testing. Not used in any production code
        import sys
        sys.path.append('/illumina/scripts/PGtools-master-dev/src')
        # from scriptsCommon import *
        #
        
        vars = Variant.objects.filter(Type='INDEL',Chr=str(job.chr),Start__gte=(job.start-500),End__lte=(job.end+500))#,Callers__in=callers) # add in callers
        mySample = '77'
        vars = vars.only('isPlat','Start','Callers__callers','Calls__N77__gt','Calls__N79__gt','Calls__N80__gt','Calls__N81__gt','Calls__N82__gt','Calls__N83__gt','Calls__N84__gt','Calls__N85__gt','Calls__N86__gt','Calls__N87__gt','Calls__N88__gt','Calls__N93__gt').as_pymongo(coerce_types=False)
        vars = vars.as_pymongo(coerce_types=False)
        print "Plat loaded "+ str(vars.count())
        varMap = {} # simple plat annotation map
        for v in vars:
            varMap[v['Start']] = int(v['isPlat'] or len(v['Calls'])>6)
    # loop over all observed indel candidates
    if mode=="freq" or mode == "BinomTest":
        f = Samfile(job.bam)

    # if mode == "BinomTest":
    #     print "Assessing TP/FP status of indel events w/ one-sided binomial exact test"

    for type in job.events.keys():
        for k in job.events[type].keys():
            myCall = job.events[type][k]
            isTrue = 0
            origType  = type
            totalObs = 0
            if mode=="PG":
                isTrue = int((varMap.has_key(k) and varMap[k])) #or calls[type][k][0]>10)
            elif mode=="freq":
                # TODO, not doing freq eval right now
                totalObs = 40
                # print f.pileup("chr20", k, k + 1)
                # print myCall
                isTrue = (1.0*myCall[0]/totalObs)>freq #or calls[type][k][0]>10)
            elif mode == "BinomTest":
                # implements a one-sided binomial exact test with the null hypothesis
                # that the binomial probability is 0.5, and the alternate hypothesis
                # that the true binomial probability is less than 0.5 (i.e. noise)
                # N.B. This should assume that if the binomial probability is higher
                # than 0.5, the null hypothesis is not rejected
                
                totalRef = 0
                for pc in f.pileup(job.chr, k, k + 1):
                    if pc.pos == k:
                        for pr in pc.pileups:
                            if pr.alignment.mapq > 20 and pr.alignment.is_paired:
                                totalObs += 1
                                if len(pr.alignment.cigar) == 1:
                                    totalRef += 1
                        break

                eventObs  = myCall[0]
                if eventObs > 0.6 * totalObs:
                    # if the event in question makes up the vast majority
                    # of the reads at the site, assume that this event is in fact
                    # the germline state, and that the reference (and any other )
                    # print "Flipping...", eventObs, totalObs - eventObs, totalObs, totalRef
                    # print type, myCall
                    originalLength = myCall[2]
                    # right now, we're not differing between insertion
                    # or deletion events of different magnitudes, so if
                    # there's already an insertion event, the additional
                    # sequence gets tacked on there.  Reset the length
                    # and the count.  Additionally, in theory, we could
                    # run in to trouble with site that have an insertion
                    # and deletion, if one type has already been processed
                    # and we decide to add more reads.  This hasn't been
                    # been observed yet, but would require some different
                    # implementation to handle
                    if type == 'del':
                        # count, event length, hpol length, seq repeated
                        obsLength = originalLength - myCall[1]
                        type      = "ins"
                        myCall[2] = obsLength
                        if job.events["ins"].has_key(k):
                            # print "has ins key"
                            # print "ins", job.events["ins"][k]
                            myCall[0] = totalRef + job.events["ins"][k][0]
                            # remove key so that it's not double-counterd
                            job.events["ins"].pop(k, None)
                        else:
                            myCall[0] = totalRef
                        # print "new", type, myCall
                    else:
                        # only other type is ins
                        obsLength = originalLength + myCall[1]
                        type      = "del"
                        myCall[2] = obsLength
                        if job.events["del"].has_key(k):
                            # print "has del key"
                            # print "del", job.events["del"][k]
                            myCall[0] = totalRef + job.events["del"][k][0]
                            job.events["del"].pop(k, None)
                        else:
                            myCall[0] = totalRef
                        # print "new", type, myCall

                binomPval = binom.cdf(myCall[0], totalObs, 0.5)
                isTrue = binomPval > minBinomPval
                # print "binom(", myCall[0], ",", totalObs,") = ", binomPval, isTrue
            else:
                isTrue = (myCall[0]>countCut)
            if isTrue: myCall.append('tp')
            else: myCall.append('fp')
            job.caseCounter.addRecord(myCall,type)
            # print job.chr, k, totalObs, myCall, type
            type = origType
    job.clearEvents() # clear the event buffer
    
        

def rec_dd(): return defaultdict(rec_dd)
class caseCounter:
    def __init__(self,totalCounts,opt):
        self.table       = rec_dd()
        self.totals      = totalCounts
        self.opt         = opt
        
    def accumulateCounts(self,cc):
#        for type in
        for ml in cc.totals.keys():
            for k in cc.totals[ml].keys():
                if not self.totals[ml].has_key(k):
                    self.totals[ml][k] = cc.totals[ml][k]
                else:
                    self.totals[ml][k] += cc.totals[ml][k]

        for ml in cc.table.keys():
            for type in cc.table[ml].keys():
                for strLength in cc.table[ml][type].keys():
                    for callLength in cc.table[ml][type][strLength].keys():
                        for repeat in cc.table[ml][type][strLength][callLength].keys():
                            if not self.table[ml][type][strLength][callLength].has_key(repeat): 
                                self.table[ml][type][strLength][callLength][repeat] = cc.table[ml][type][strLength][callLength][repeat]
                            else:
                                self.table[ml][type][strLength][callLength][repeat] += cc.table[ml][type][strLength][callLength][repeat]
        
    def addRecord(self,case,myType):
        counts, callLength, strLength, repeat, motifLength, platCase = case
        if not self.table[motifLength][myType][strLength][callLength].has_key(repeat):
            self.table[motifLength][myType][strLength][callLength][repeat] = Counter()
        self.table[motifLength][myType][strLength][callLength][repeat][platCase] += counts

    def getByHpolandLength(self, strLength, length, myType='del', motifLen = 1):
        res = Counter()
        for repeat in self.table[motifLen][myType][strLength][length].keys():
            res += self.table[motifLen][myType][strLength][length][repeat]
        return res

    def getByHpolandRepeat(self, hpol, repeat, myType='del'):
        res = {}
        for length in self.table[motifLen][myType][hpol].keys():
            for repeat in self.table[motifLen][myType][hpol][length].keys():
                if not res.has_key(repeat): res[repeat] = Counter()
                res[repeat] += self.table[motifLen][myType][hpol][length][repeat]
        return res[repeat]

    def getByStrLen(self, strLen=5, myType='del', motifLen = 1):
        res = Counter()
        for length in self.table[motifLen][myType][strLen].keys():
            res += self.getByHpolandLength(strLen, length, myType, motifLen)
        return res

    def writeTable(self, myType = 'del', motifLen = 1):
        print "TABLE FOR " + myType + " in " + str(motifLen) + " bp STRs"
        print "------------------------------"
        for hpol in sorted(self.table[motifLen][myType].keys()):
            print "Hpol: " + str(hpol)
            for l in sorted(self.table[motifLen][myType][hpol].keys()):
                print "    Indel length: " + str(l)
                print self.getByHpolandLength(hpol, l, myType, motifLen)
                for repeat in sorted(self.table[motifLen][myType][hpol][l]):
                    count  = self.table[motifLen][myType][hpol][l][repeat]
                    print "        " + repeat + ": " +  str(count) + " "  + str(count['fp']*1.0/(count['fp']+count['tp']))
    def writePars(self):
        print "Pars"

    def addModel(self, name, myType, motifLen):
        self.modelList  = ["case","props","conf", "obs","error", "fp", "tp"]
        for m in self.modelList:
            self.parsSTR[motifLen][myType][name][m] = []
        return self.parsSTR[motifLen][myType][name]

    def printModels(self,myType='del'):
        for mName in self.parsHpol[myType].keys():
            mName

    def simplePlot(self, motifLen):
        print "making simple plot"
        base = {'ins':{},'del':{}}
        for k in base.keys():
            print "Case " + k
            print self.parsSTR[motifLen][k]
            base[k]['Nova'] = self.parsSTR[motifLen][k]['simple']
        name = "simple{}".format(motifLen)
        errorPlot(base,name)

    def lengthPlot(self,cuts=[[1,2],[2,3],[3,5],[5,8],[8,12]]):
        base = {'ins':{},'del':{}}
        for k in base.keys():
            for i,c in enumerate(cuts):
                name = self.getModelName(*c)
                legendName = "INDEL("
                if c[0]==(c[1]-1): legendName += str(c[0])
                elif (i+1)==len(cuts): legendName += str(c[0]) + "+"
                else: legendName  += str(c[0]) +"-"+str(c[1]-1)
                legendName +=")"
                base[k][legendName] = self.parsHpol[k][name]
        errorPlot(base,"length",plotFit=0,plotHiseq=1)
#        return self

    def basePlot(self,myBases=["A","C","T","G"]):
        base = {'ins':{},'del':{}}
        for k in base.keys():
            for b in myBases:
                base[k][b] = self.parsHpol[k][b]
        errorPlot(base,"base",plotFit=0,plotHiseq=0,indvObs=1,plotPrev=1)

    def basePlotWithFit(self,myBases=["T","C"]):
        base = {'ins':{},'del':{}}
        for k in base.keys():
            for b in myBases:
                base[k][b] = self.parsHpol[k][b]
        errorPlot(base,"baseFit",plotFit=1,plotHiseq=0,indvObs=1,plotPrev=0)

#    def fitSigmoid(self,x,y,end=12):
#                    # fit curve
#    #        print res[k]['errorProb'].values()[1:20]
#        xdata = np.array(x[0:end])
#        ydata = np.array(y[0:end])
#    #
#        p_guess=(2e-03,1e+01,2e+0,8e-06)
#        (p, cov, infodict, mesg, ier)=scipy.optimize.leastsq(residuals,p_guess,args=(xdata,ydata),full_output=1)
#        xfit = np.linspace(-1,35, 1000)
#        yfit = sigmoid(xfit,*p)
#        return [xfit,yfit]

    def getModelName(self,start,end):
        return "l_"+  str(start) + "_" + str(end)

    def calcLengthProb(self,start,end,myType='del', byBase=''):
        '''
        Determine the error prop 
        :param start: min indel length
        :param end: max indel length
        :param byBase: Only do the model for this base
        '''
        modelName = self.getModelName(start,end)
        if byBase: modelName + "_" + byBase
        model = self.addModel(modelName, myType)
        for hpol in xrange(1,35):
            I_all = self.getByHpol(hpol, myType)
            I = I_all["tp"] #+I_all["fp"]
            temp = Counter()
            for i_temp in xrange(start,end):
                if len(byBase):
                    e=1
                else:
                    temp    += self.getByHpolandLength(hpol,i_temp, myType)
            I_l     = temp["tp"] #+temp["fp"]
            I_l_e   = temp["fp"]
            I_e     = I_all["fp"]
            if len(byBase):
                N       = self.totals[byBase][hpol]*1.0
            else:
                N       = self.totals['tot'][hpol]*1.0
            p_e     = 1.0*I_all["fp"]/N
            p_l = 1.0*(I_l+1)/(I+1) # add in baseline count
            p_l_given_e = 1.0*(I_l_e)/(I_e)
            post = p_l_given_e*p_e/p_l
            num = temp["fp"]
            den = (temp["tp"]+1)
            scale = 1.0*(I_all["tp"]+1)/N
#            tempConf =  binom_interval(num,den)
            post2 = 1.0*scale*num/den

#                print modelName + ": " + str(post)
#            conf = (tempConf[0]*scale,tempConf[1]*scale)
#                print conf
            model['case'].append(hpol)
            model['props'].append(post2)
#            model['conf'].append(conf)
            model['obs'].append(N)
#            model['error'].append(abs(model['conf'][-1][0]-model['props'][-1]))

    def calcPars(self, myType='del', lengthCases=[[1,5],[1,4],[2,4],[2,5],[1,2],[1,3]]):
        try: self.parsSTR
        except: self.parsSTR   = rec_dd()
        for ml, info in self.totals.iteritems():
            # print ml
            # print ml, info
            d = self.addModel('simple', myType, ml)
            for strLen in xrange(1, 35):
                # print strLen, self.totals[ml]['tot'][strLen]
                tots  = self.totals[ml]['tot'][strLen] * 1.0 # normalize observations by hpol length
                count = self.getByStrLen(strLen, myType, ml)
                d['case'].append(strLen)
                if tots > 0: 
                    d['props'].append(1.0 * count['fp'] / (tots - count['tp']))
                    d['conf'].append(binom_interval(count['fp'],tots))
                    d['error'].append(abs(d['conf'][-1][0]-d['props'][-1]))
                    d['obs'].append(tots)
                    d['fp'].append(count['fp'])
                    d['tp'].append(count['tp'])
            
                else: 
                    d['props'].append(0.1)
                    d['conf'].append((0.1,0.1))
                    d['error'].append(0.1)
                    d['obs'].append(10)
                    d['fp'].append(1)
                    d['tp'].append(1)
        
        # custom fitting
#        end =12
#        if myType=='ins': end =20
#        d['fit'] = self.fitSigmoid(d['case'], d['props'],end)
            
        # Calculate the posterior by indel length
#        for i in lengthCases:
#            self.calcLengthProb(i[0],i[1],myType=myType)
        
        # calc base specific rates
#        for base in ["A","C","T","G"]:
#            d = self.addModel(base, myType)                
#            for hpol in xrange(1,35):
#                count = self.getByHpolandRepeat(hpol, base, myType)
#                if self.totals[base][hpol]<10000: break
#                tots = self.totals['tot'][hpol]*1.0/hpol # normalize observations by hpol length
#                d['case'].append(hpol)
#                d['props'].append(1.0*count['fp']/(tots+count['tp']))
#                    
#                d['conf'].append(binom_interval(count['fp'],tots))
#                d['error'].append(abs(d['conf'][-1][0]-d['props'][-1]))
#                d['obs'].append(tots)
#            d['fit'] = self.fitSigmoid(d['case'], d['props'],15)
