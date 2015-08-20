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
import numpy as np
from scipy.stats import beta

def verifyInput(opt):
    '''
    Verify that input parameters are as expected
    '''
    # check sample name
    indv = ('77','78','79','80','81','82','83','84','85','86','87','88','93')
    if opt.naiveProject and not opt.countUnspanned:
        return overAndOut("naiveProject = True requires countUnspanned = True; please fix options")

    if not opt.inputSample: return overAndOut("Sample name not specified")
    if not opt.inputSample.replace("NA128",'') in indv: return overAndOut("Unknown sample name: " + opt.sampleName)
    logging.info("Got input from sample: " + opt.inputSample)

    # check reference is there
    if opt.reference and not fileSimple(opt.reference): return 0
    logging.info("Using reference: " + opt.reference)
    # check that we have either bam or restart file
    if not opt.inputBam and not opt.startFromFile: return overAndOut("Must specicify either and input bam or input restart file")
    # check input bam
    if not opt.startFromFile and not opt.inputBam: return overAndOut("Input bam does not specified")
    if not fileCheckBam(opt.inputBam): return 0
    # check input json restart file
    if opt.startFromFile and not fileSimple(opt.startFromFile): return 0
    if opt.inputBam:
        logging.info("Using input bam: " + opt.inputBam)
    else:
        logging.info("Using restart file: " + opt.startFromFile)

    # check if there is not  an output path set spit it out in the bam dir/indelFitOut
    if not opt.output and opt.inputBam:
        d = os.path.dirname(opt.inputBam)
        opt.output = os.path.join(d,'indelErrorOutput')
    if not opt.output and opt.startFromFile:
        d = os.path.dirname(opt.startFromFile)
        opt.output = os.path.join(d,'indelErrorOutput')
    # make a missing dir
    if not os.path.exists(opt.output):
            os.mkdir(opt.output)
    logging.info("Using output dir: " + opt.output)

    try:
        opt.cores = int(opt.cores)
    except:
        return overAndOut("Invalid core number " + opt.cores)
    return 1

class ReferenceBuffer:
    def __init__(self, fetch, c, start, end, maxMotifLen = 2, force=1):
        self.chr             = str(c)
        self.start           = start
        self.end             = end
        self.maxMotifLen     = maxMotifLen
        self.repMap          = {} # store length and motif for each base in an STR
        self.repStartsMap    = {}        # store length and motif by STR start
        self.entropyMap      = {}        # store entropy information TODO
        self.count           = [Counter() for x in xrange(maxMotifLen)] # summary stats, # of STRs observed of each length
        self.fetch           = fetch     # a fetch instance from FastaFile

    def buffer(self):
        # look for repeats of specified types in reference sequence
        i           = self.start
        current     = np.empty(self.maxMotifLen, dtype = "S1")
        match       = np.empty(self.maxMotifLen, dtype = "b")
        inSTR       = np.empty(self.maxMotifLen, dtype = "b")
        motif       = np.empty(self.maxMotifLen, 
                               dtype = "S{}".format(self.maxMotifLen))
        l           = self.resetLength()

        current.fill("N")
        match.fill(False)
        inSTR.fill(False)

        for c in self.fetch:
            # no STR detected yet, so look at everything
            for j, m in enumerate(current):
                if c.upper() == m:
                    # in a repetitive region
                    l[j] += 1
                    match[j] = True
                    if not inSTR[j]:
                        inSTR[j] = True
                        motif[j] = "".join(current[j:])
                else:
                    match[j] = False
            # print i, c, current, l, match, inSTR, motif

            # we have exited an STR if we used to be in an STR and the current
            # base does not match the one at the STR's index
            if np.any(np.logical_and(np.logical_not(match), inSTR)):
                # no longer in some STR
                if not np.any(current == 'N') and not c.upper() == "N":
                    for j in reversed(xrange(self.maxMotifLen)):
                        if inSTR[j] and not match[j]:
                            motifLength = self.maxMotifLen - j
                            length = l[j]
                            start  = np.uint(i - length)
                            if length / motifLength >= 2:
                                # e.g., a dinucleotide of length 3 is meaningless
                                # print "adding", length, motif[j], j, motifLength
                                self.repStartsMap[start]  = [length, motif[j], motifLength]
                                for t in xrange(start, i):
                                    if self.repMap.has_key(t):
                                        # print t, length, motif[j], motifLength, self.repMap[t]
                                        # if the base is already part of an STR, keep the STR
                                        # with the longer period (instead of having both STRS
                                        # marked there)
                                        if self.repMap[t][2] < motifLength:
                                            self.repMap[t] = [length, motif[j], motifLength]
                                        elif self.repMap[t][0] < length:
                                            self.repMap[t] = [length, motif[j], motifLength]
                                    else:
                                        self.repMap[t] = [length, motif[j], motifLength]
                                        self.count[motifLength - 1][length] +=1
                                self.repMap[i] = [1, c.upper(), 1]
                            else:
                                for t in xrange(start, i):
                                    if not self.repMap.has_key(t):
                                        kmerPos = np.uint((t - start) % motifLength)
                                        self.repMap[t] = [1, motif[j][kmerPos], 1]
                                self.repMap[i] = [1, c.upper(), 1]
                            break
                    inSTR[0:(j + 1)] = False
                    motif[0:(j + 1)] = None
                    if j > 0:
                        l[0:(j + 1)] = [self.maxMotifLen - x for x in xrange(0, j + 1)]
                    else:
                        l[0] = self.maxMotifLen
            else:
                self.repMap[i] = [1, c.upper(), 1] #stores non-repeat bases as hpol length 1
            current    = np.roll(current, -1)
            current[1] = c.upper()
            i += 1
        # check to see if we terminate in an STR
        if np.any(inSTR):
            # no longer in an STR
            if not np.any(current == 'N') and not c.upper() == "N":
                for j in reversed(xrange(self.maxMotifLen)):
                    if inSTR[j]:
                        motifLength = self.maxMotifLen - j
                        length = l[j]
                        start  = np.uint(i - length)
                        if length / motifLength >= 2:
                            # e.g., a dinucleotide of length 3 is meaningless
                            self.repStartsMap[start]  = [length, motif[j], motifLength]
                            for t in xrange(start, i):
                                if self.repMap.has_key(t):
                                    # print t, length, motif[j], motifLength, self.repMap[t]
                                    # if the base is already part of an STR, keep the STR
                                    # with the longer period (instead of having both STRS
                                    # marked there)
                                    if self.repMap[t][2] < motifLength:
                                        self.repMap[t] = [length, motif[j], motifLength]
                                    elif self.repMap[t][0] < length:
                                        self.repMap[t] = [length, motif[j], motifLength]

                                else:
                                    self.repMap[t] = [length, motif[j], motifLength]                                    
                        else:
                            for t in xrange(start, i):
                                if not self.repMap.has_key(t):
                                    kmerPos = np.uint((t - start) % motifLength)
                                    self.repMap[t] = [1, motif[j][kmerPos], 1]
                        break
            else:
                self.repMap[i] = [1, c.upper(), 1] #stores non-repeat bases as hpol length 1
        # for k in sorted(self.repMap.keys()):
        #     print k, self.repMap[k]
        # print self.count
    #        print c.upper()

    def resetLength(self):
        l = np.empty(self.maxMotifLen, dtype = np.uint)

        for x in xrange(self.maxMotifLen):
            l[x] = self.maxMotifLen - x

        return l

def binom_interval(success, total, confint=0.95):
    if not total: return (0,0)
    quantile = (1 - confint) / 2.
    lower = beta.ppf(quantile, success, total - success + 1)
    upper = beta.ppf(1 - quantile, success + 1, total - success)
    return (lower, upper)

def sigmoid(x,a,b,c,d):
#    a,b,c,d=p
    y = a / (1 + np.exp((b-x)/c))+d
#    est.p<-coef.sig[1]/(1+exp((coef.sig[2]-x)/coef.sig[3]))
    return y

def residuals(p,x,y):
    return y - sigmoid(x,*p)

def entropy(s):
    '''
    Given a sequence calculate the entropy
    :param s:
    '''
    p, lns = Counter(s), float(len(s))
    return -sum( count/lns * math.log(count/lns, 2) for count in p.values())

def RC(seq):
    if not seq:
        logging.info("Empty RC string")
        sys.exit(-1)
    # if len(seq) > maxMotifLen:
    #     logging.info("RC called on seq longer than ", maxMotifLen, ":", seq)
    #     sys.exit(-1)

    rc = []

    for c in seq:
        if c.upper() == "A":
            rc.append("T")
        if c.upper() == "C":
            rc.append("G")
        if c.upper() == "G":
            rc.append("C")
        if c.upper() == "T":
            rc.append("A")
        if c.upper() == "N":
            rc.append("N")
    rc = reversed(rc)
    return "".join(rc)

def validateRecordByFreq(freqCut=0.07):
    myFreq = 0.01

    return 1
