
# produce values for legacy error model
# indel error model parameters for P(error) = Ax+Bx^C, where x=hpol_len
# note that fit does not cover length 1 deletions,
# for which the estimated value is instead provided directly
def prevModel(stop=35):
    def f_del(x):
#        if x==1:return 3.00057e-6 
        return (1.09814e-5)*x+(5.19742e-10)*x**(6.99256) 
    def f_ins(x):
        return (5.03824e-7)*x+(3.30572e-10)*x**(6.99777)
    res = {'ins':{},'del':{}}
    x = np.linspace(1,stop, 100)
    dels = f_del(x)
    dels[0] = 3.00057e-6 # special case for hpol=1
    res['del'] = [x,dels]  
    res['ins'] = [x,f_ins(x)]
    return res

def errorPlot(data,name,plotPrev=1,plotFit=1,plotHiseq=0,plotNova=0,indvObs=0):#,eventCounts,figDict={},color='lightblue',label='',plotFit=0):
#    if not len(figDict.keys()):
    fig2 = plt.figure(figsize=(15.0, 8.0)) 
#    fig,(figDel, figIns) = plt.subplots(1, 2, 0, 1)#, squeeze, subplot_kw)
    figDel = fig2.add_subplot(121)
    figIns = fig2.add_subplot(122)
    figDict = {"del":figDel,"ins":figIns,"del2":figDel.twinx(),"ins2":figIns.twinx(),}
    colors = ['Red','Blue','Green','Gold',"Cyan"]#brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
    prevM   = prevModel(25)
    hiseqM  = hiseqBase(25)
    novaM   = novaBase(25)
    for case in ['del','ins']:
        case2 = case + "2"
        for index,dset in enumerate(data[case].keys()):
            myData = data[case][dset]
            if index==0 and not indvObs:figDict[case].plot(myData['case'],myData['obs'],'--',color='lightgrey', lw = 1,label="Count")
            if indvObs: figDict[case].plot(myData['case'],myData['obs'],'--',color=colors[index], lw = 1)
            figDict[case2].errorbar(myData['case'],myData['props'],yerr=myData['error'], fmt='o',color=colors[index],label=dset)
            if plotFit:
                figDict[case2].plot(myData['fit'][0],myData['fit'][1],'--',color=colors[index],label='Sigmoid fit ' + dset)
        figDict[case].set_xlabel('Homopolymer length')
        figDict[case].set_ylabel('Observations')
        figDict[case2].set_ylabel('P(error)')
        figDict[case].set_title(case.upper())     
        
        # set plot previous model 
        if plotPrev: 
            figDict[case2].plot(prevM[case][0],prevM[case][1],'--',color='blue',label="Current model")
        if plotHiseq:
            figDict[case2].plot(hiseqM[case][0],hiseqM[case][1],'--',color='red',lw=2,label="HiSeq fit")
        if plotNova:
            figDict[case2].plot(novaM[case][0],novaM[case][1],'--',color='red',lw=2,label="Nova fit")
        # set plot descriptors
            
        figDict[case].set_yscale('log')
        figDict[case2].set_yscale('log')
        figDict[case].set_xlim([0,35])
        if case=='del':figDict[case].legend()
        else: figDict[case2].legend(loc=2)
    plt.savefig('/home/mkallberg/workspace/PGtools/data/indelModel/plots_nova/' + name + '.png')
    plt.show()
    return figDict