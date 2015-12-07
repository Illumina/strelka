#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2016 Illumina, Inc.
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

import numpy as np
import matplotlib.pyplot as plt
import pylab
import scipy
from scipy.interpolate import interp1d
from scipy import interpolate

def hiseqBase(end=35):
    p_del = [  1.16392237e-03,   9.97313200e+00,  7.69505224e-01, 1.74938546e-05]
    p_ins = [  9.74877227e-04,   9.98186051e+00,  8.80113078e-01,  1.49391554e-05]
    xfit = np.linspace(-1,35, 100)
    yfit_del = sigmoid(xfit,*p_del)
    yfit_ins = sigmoid(xfit,*p_ins)
    return {"ins":[xfit,yfit_ins],"del":[xfit,yfit_del],}
    
def novaBase(p=(),end=35):
    p_del = [  2.63234901e-03, 9.97963446e+00, 8.54654966e-01, 7.71920706e-05]
    p_ins = [  2.33934311e-03, 9.70614686e+00, 9.01825761e-01, 4.42294198e-05]
    xfit = np.linspace(-1,35, 100)
    yfit_del = sigmoid(xfit,*p_del)
    yfit_ins = sigmoid(xfit,*p_ins)
    return {"ins":[xfit,yfit_ins],"del":[xfit,yfit_del],}

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

def errorPlot(data,name,plotPrev=0,plotFit=0,plotHiseq=0,plotNova=0,indvObs=0):#,eventCounts,figDict={},color='lightblue',label='',plotFit=0):
#    if not len(figDict.keys()):
    fig2 = plt.figure(figsize=(15.0, 8.0))
#    fig,(figDel, figIns) = plt.subplots(1, 2, 0, 1)#, squeeze, subplot_kw)
    figDel = fig2.add_subplot(121)
    figIns = fig2.add_subplot(122)
    figDict = {"del":figDel,"ins":figIns,"del2":figDel.twinx(),"ins2":figIns.twinx(),}
    colors = ['Red','Blue','Green','Gold',"Cyan"]#brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
#    prevM   = prevModel(25)
#    hiseqM  = hiseqBase(25)
#    novaM   = novaBase(25)
    for case in ['del','ins']:
        case2 = case + "2"
        for index,dset in enumerate(data[case].keys()):
            myData = data[case][dset]
#            if index==0 and not indvObs:figDict[case].plot(myData['case'],myData['obs'],'--',color='lightgrey', lw = 1,label="Count")
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
    plt.savefig('/home/mbekritsky/workspace/starka/error_model/' + name + '.png')
    plt.show()
    return figDict