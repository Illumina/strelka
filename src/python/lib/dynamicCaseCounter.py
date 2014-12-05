#
# Starka
# Copyright (c) 2009-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

from collections import *
from pysam import *
from dynamicModelUtil import *

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
        self.refOffset      = 500                   # additional up and downstream BPs to buffer reference
        self.verbose        = 1
        self.bufferRef()
        self.setDatastructures()

    def bufferRef(self):
        refStart = self.start-self.refOffset
        if refStart < 0: self.refStart = 0
        refEnd = self.end+self.refOffset
        fa = Fastafile(self.reference)
        self.refInterval = fa.fetch(reference=self.chr,start=refStart,end=refEnd)
        self.buffer = ReferenceBuffer(self.refInterval,self.numericChr(),start=refStart,end=refEnd)
        self.buffer.buffer()

    def setDatastructures(self):
        self.calls       = []
        self.events      = {'ins':{},'del':{}}
        self.eventCount  = {'tot':Counter(),'ins':Counter(),'del':Counter(),}
        for base in ['A','C','T','G','N']: self.eventCount[base] = Counter()

    def run(self):
        readRecords(self)

    def numericChr(self):
        return int(self.chr.replace('chr',''))

    def __repr__(self):
        return str(self)
    def __str__(self):
        return self.chr + ":" + str(self.start) + "-" + str(self.end) + " - " + str(self.eventCount)

# cigar str numbering def. by pysam
#BAM_CMATCH=0, BAM_CINS=1,
#DEF BAM_CDEL       = 2
#DEF BAM_CREF_SKIP  = 3
#DEF BAM_CSOFT_CLIP = 4
#DEF BAM_CHARD_CLIP = 5
#DEF BAM_CPAD       = 6
def readRecords(job,seqContext='All'):

#    mychr, start, end, countEveryBase, ignoreStrand, naiveProject, countUnspanned = job
    print "Collecting data "
    f = Samfile(job.bam)
    i = 0
    for r in f.fetch(job.chr,job.start,job.end):
        if r.mapq>20 and r.is_paired:# only concerned with reads of high quality so we don't count alignment errors
            myEvent = 'm'       # treat clipped reads as matches
            myPos = r.pos
            if len(r.cigar)>1:
                for event in r.cigar:
                    if event[0]==1 or event[0]==2:
                        myEvent = 'ins'
                        if event[0]==2: myEvent = 'del'
#                        print myEvent + "  at " + str(myPos)
                        break
                    myPos += event[1]
#                # add in information
                if job.buffer.hpolMap.has_key(myPos):
                    hpoleCase,repeat = job.buffer.hpolMap[myPos]
                    if r.is_reverse and not job.ignoreStrand:
                        repeat = RC(repeat)
                    if not myEvent=='m':
                        if not job.events[myEvent].has_key(myPos):
#                                                    # count, event length, hpol length, seq repeatted
##                            print [1,event[1],hpoleCase,repeat]
                            if seqContext=='All' or seqContext==repeat.upper():
                                myLength = event[1]
                                job.eventCount[myEvent][hpoleCase] +=1
                                job.events[myEvent][myPos] = [1,myLength,hpoleCase,repeat]
                        else:
                            if seqContext=='All' or seqContext==repeat.upper():
                                job.eventCount[myEvent][hpoleCase] +=1
                                job.events[myEvent][myPos][0] +=1


            if job.countUnspanned: # i.e. skip logic to check whether alignment extends past the homopolymer
                # we can either count only reference bases actually covered or we can project based on nominal length of the read
                # first, set according to the actual alignment
                countCoverStart = r.pos+1  # r.pos+1 so we know we have at least one base before the hpol; r.alen handles trimming and indel adjustments
                countCoverEnd = r.pos+r.alen
                if job.naiveProject:
                    # if user so requests, reset to projection based on read length
                    countCoverStart = r.pos
                    countCoverEnd = r.pos+100 # assume a 100 bp readlength

                for ii in xrange(countCoverStart,countCoverEnd):
                    if job.buffer.hpolStartsMap.has_key(ii) or (job.countEveryBase and job.buffer.hpolMap.has_key(ii)):
                        hpoleCase,repeat = job.buffer.hpolMap[ii] # note hpol has same data as hpolStarts where hpolStarts is defined
                        if r.is_reverse and not job.ignoreStrand:
                            repeat = RC(repeat)
                        if not job.eventCount.has_key(repeat): job.eventCount[repeat] = Counter()
                        job.eventCount['tot'][hpoleCase]  +=1
                        job.eventCount[repeat][hpoleCase] +=1

            else:
                refOffset = job.start-job.refOffset
                alignedPairs = r.aligned_pairs
                numpairs = len(alignedPairs)

                prevMatch = False
                for ii in xrange(0,numpairs):
                    (x,y) = alignedPairs[ii]
                    if x != None:
                        readbase = r.query[x]
                    else:
                        readbase = '-'
                    if y != None:
                        refbase = job.refInterval[y-refOffset].upper()
                    else:
                        refbase = '-'


                    rpos = None
                    qpos = None

                    ###### count only events we can have confidence in:
                    ###### require that the alignment extend beyond a homopolymer by at least one base in each direction
                    ###### require that the base after the homopolymer (in the direction of sequencing) match

                    # if current position is the start of a homopolymer run
                    if job.buffer.hpolStartsMap.has_key(y) or (job.countEveryBase and job.buffer.hpolMap.has_key(y)):
                        hpoleCase,repeat = job.buffer.hpolMap[y] # note hpol[y] == hpolStarts[y] if the latter is defined

                        # make sure we align past the end of the hpol regions
                        if y+hpoleCase >= r.pos+r.alen:
                            continue

                        ok = False # do not use this alignment position unless we satisfy conditions below

                        # if reverse strand alignment:
                        if r.is_reverse:
                            # require previous position to match (and hence also that the alignment goes beyond the homopolymer
                            if prevMatch:
                                ### we can use it!
                                if not job.ignoreStrand: repeat = RC(repeat)
                                ok = True
                        else: # forward strand alignment
                            # require this to not be the first aligned base of the read
                            if ii > 0:
                                rpos = None
                                qpos = None
                                # find alignment of first base past homopolymer, or an interior indel
                                for iii in xrange(ii+1,numpairs):
                                    rpos = alignedPairs[iii][1]
                                    qpos = alignedPairs[iii][0]
                                    if (not rpos) or (not qpos) or rpos == y+hpoleCase:
                                        break
                                # require that we be at the right position and it be an alignment of bases and they match
                                if rpos == y+hpoleCase and qpos and job.refInterval[rpos-refOffset].upper() == r.query[qpos]:
                                    ### we can use it!
                                    ok = True

                        if ok:
                            if not job.eventCount.has_key(repeat): job.eventCount[repeat] = Counter()
                            job.eventCount['tot'][hpoleCase] +=1
                            job.eventCount[repeat][hpoleCase] +=1

                    if readbase == refbase:
                        prevMatch = True
                    else:
                        prevMatch = False

            i+=1
#    events = validateRecords(events,job)
#    return events,eventCount



def rec_dd(): return defaultdict(rec_dd)
class caseCounter:
    def __init__(self,totalCounts,opt):
        self.table      = rec_dd()
        self.totals     = totalCounts
        self.opt        = opt
    def addRecord(self,case,myType):
        counts,callLength,hpolLength,repeat,platCase = case
        if not self.table[myType][hpolLength][callLength].has_key(repeat): self.table[myType][hpolLength][callLength][repeat] = Counter()
        self.table[myType][hpolLength][callLength][repeat][platCase] += counts

    def getByHpolandLength(self,hpol,length,myType='del'):
        res = Counter()
        for repeat in self.table[myType][hpol][length].keys():
            res += self.table[myType][hpol][length][repeat]
        return res

    def getByHpolandRepeat(self,hpol,repeat,myType='del'):
        res = {}
        for length in self.table[myType][hpol].keys():
            for repeat in self.table[myType][hpol][length].keys():
                if not res.has_key(repeat): res[repeat] = Counter()
                res[repeat] += self.table[myType][hpol][length][repeat]
        return res[repeat]

    def getByHpol(self,hpol=5,myType='del'):
        res = Counter()
        for length in self.table[myType][hpol].keys():
            res += self.getByHpolandLength(hpol, length, myType)
        print res
        return res

    def writeTable(self,myType='del'):
        print "TABLE FOR " + myType
        print "------------------------------"
        for hpol in sorted(self.table[myType].keys()):
            print "Hpol: " + str(hpol)
            for l in sorted(self.table[myType][hpol].keys()):
                print "    Indel length: " + str(l)
                print self.getByHpolandLength(hpol, l, myType)
                for repeat in sorted(self.table[myType][hpol][l]):
                    count  = self.table[myType][hpol][l][repeat]
                    print "        " + repeat + ": " +  str(count) + " "  + str(count['fp']*1.0/(count['fp']+count['tp']))
    def writePars(self):
        print "Pars"

    def addModel(self,name,myType):
        self.modelList  = ["case","props","conf", "obs","error"]
        for m in self.modelList:
            self.parsHpol[myType][name][m] = []
        return self.parsHpol[myType][name]

    def printModels(self,myType='del'):
        for mName in self.parsHpol[myType].keys():
            mName

    def simplePlot(self):
        print "making simple plot"
        base = {'ins':{},'del':{}}
        for k in base.keys():
            base[k]['Nova'] = self.parsHpol[k]['simple']
        errorPlot(base,"simple")

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
        errorPlot(base,"base",plotFit=0,plotHiseq=1,indvObs=1,plotPrev=1)

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
                N       = self.totals[byBase][hpol]*1.0/hpol
            else:
                N       = self.totals['tot'][hpol]*1.0/hpol
            p_e     = 1.0*I_all["fp"]/N
            p_l = 1.0*(I_l+1)/(I+1) # add in baseline count
            p_l_given_e = 1.0*(I_l_e)/(I_e)
            post = p_l_given_e*p_e/p_l
            num = temp["fp"]
            den = (temp["tp"]+1)
            scale = 1.0*(I_all["tp"]+1)/N
            tempConf =  binom_interval(num,den)
            post2 = 1.0*scale*num/den

#                print modelName + ": " + str(post)
            conf = (tempConf[0]*scale,tempConf[1]*scale)
#                print conf
            model['case'].append(hpol)
            model['props'].append(post2)
            model['conf'].append(conf)
            model['obs'].append(N)
            model['error'].append(abs(model['conf'][-1][0]-model['props'][-1]))

    def calcPars(self,myType='del'):
        try: self.parsHpol
        except: self.parsHpol   = rec_dd()
        d = self.addModel('simple', myType)
        for hpol in xrange(1,35):
            tots = self.totals['tot'][hpol]*1.0/hpol # normalize observations by hpol length
            count = self.getByHpol(hpol, myType)
            d['case'].append(hpol)
            d['props'].append(1.0*count['fp']/tots)
            d['conf'].append(binom_interval(count['fp'],tots))
            d['error'].append(abs(d['conf'][-1][0]-d['props'][-1]))
            d['obs'].append(tots)

        # custom fitting
        end =12
        if myType=='ins': end =20
        d['fit'] = self.fitSigmoid(d['case'], d['props'],end)

        # Calculate the posterior by indel length
        for i in [[1,2],[2,3],[3,5],[5,8],[8,12]]:
            self.calcLengthProb(i[0],i[1],myType=myType)

        # calc base specific rates
        for base in ["A","C","T","G"]:
            d = self.addModel(base, myType)
            for hpol in xrange(1,35):
                count = self.getByHpolandRepeat(hpol, base, myType)
                if self.totals[base][hpol]<10000: break
                tots = self.totals['tot'][hpol]*1.0/hpol # normalize observations by hpol length
                d['case'].append(hpol)
                d['props'].append(1.0*count['fp']/(tots+count['tp']))

                d['conf'].append(binom_interval(count['fp'],tots))
                d['error'].append(abs(d['conf'][-1][0]-d['props'][-1]))
                d['obs'].append(tots)
            d['fit'] = self.fitSigmoid(d['case'], d['props'],15)
