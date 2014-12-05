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

def verifyInput(opt):
    '''
    Verify that input parameters are legit
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
    def __init__(self,fetch,c,start,end,force=1):
        self.chr             = str(c)
        self.start           = start
        self.end             = end
        self.hpolMap         = {}        # store hpol information
        self.hpolStartsMap   = {}        # store hpol information
        self.dinucMap        = {}        # store di-nuc information TODO
        self.entropyMap        = {}        # store entropy information TODO
        self.count           = Counter() # summary stats
        self.fetch           = fetch

    def buffer(self):
        i           = self.start
        current     = "N"
        l           = 1
        for c in self.fetch:
            if c.upper()==current:
                l +=1
            else:
                if not current=='N':
                    self.hpolStartsMap[i-l]  = [l,current]
                    for t in xrange(i-l,i):
                        self.hpolMap[t]  = [l,current]
                    self.count[l]    +=1
#                    print 'pos ' + str(i-l) + " " +  current + ': ' + str(l)
                l = 1
                current = c.upper()
            i+=1
    #        print c.upper()

def binom_interval(success, total, confint=0.95):
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

def RC(base):
    if len(base) != 1:
        logging.info("RC called on bad string: " + base)
        sys.exit(-1)
    if base.upper() == "A":
        return "T"
    if base.upper() == "C":
        return "G"
    if base.upper() == "G":
        return "C"
    if base.upper() == "T":
        return "A"
    if base.upper() == "N":
        return "N"
    logging.info("RC called on bad character: " + base)
    sys.exit(-1)

def validateRecordByFreq(freqCut=0.07):
    myFreq = 0.01

    return 1
