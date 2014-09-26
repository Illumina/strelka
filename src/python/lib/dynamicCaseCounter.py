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