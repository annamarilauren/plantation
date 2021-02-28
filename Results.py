# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 10:34:00 2016

@author: vagrant
"""
import pickle
import os
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
from dateutil.relativedelta import relativedelta


def C_balance(gy, de, ro, pe):
    #temp
    s=(ro.SH + ro.H + ro.LH + ro.F + ro.L)*10000.0
    #------
    soil = (de.cwdSto + (ro.SH + ro.H + ro.LH + ro.F + ro.L)*10000.0 + pe.peatStorage)*ro.massToC 
    standTot = (gy.BiBark + gy.BiBranch  + gy.BiFoliage + gy.RootMass + gy.weedAbove + gy.weedBelow)*ro.massToC 
    standTot2 = (gy.BiStem + gy.BiBark + gy.BiBranch  + gy.BiFoliage + gy.RootMass + gy.weedAbove + gy.weedBelow)*ro.massToC 
    
    bal = (soil[-1] -soil[0]) + (standTot[-1] - standTot[0])
    soil_balance = soil[-1] -soil[0]
    total_balance =    (soil[-1] -soil[0]) + (standTot2[-1] - standTot2[0])
    #print bal, (s[-1]-s[0])*ro.massToC ,    (pe.peatStorage[-1]-pe.peatStorage[0])*ro.massToC, sum(ro.miner)*ro.massToC 
    return soil_balance, bal, total_balance

def save_results(gy, de,ro, pe, N=None, P=None, K=None, fold= None, fname =None):
    if fold==None:
        fold="C:\Apps\WinPython-64bit-2.7.10.3\IPEWG\\"
    if fname==None:
        fname = "PDSS_results.p"
    data={}
    data['gy'] = gy.__dict__
    data['de'] = de.__dict__
    data['ro'] = ro.__dict__
    data['pe'] = pe.__dict__
    if N != None: data['N'] = N.__dict__
    if P != None: data['P'] = P.__dict__   
    if K != None: data['K'] = K.__dict__    
    with open(fold+fname, 'wb') as wfp:
        pickle.dump(data, wfp)
    print ('PDSS simulation saved to ', fname)
    for i in plt.get_fignums():
        plt.figure(i)
        plt.savefig(fold+'figure%d.png' % i)    
    
def print_test(fold=None, fname=None):    
    if fname==None:
        fname = "PDSS_results.p"
    with open(fold+fname, 'rb') as wfp:
        aaa=pickle.load(wfp)
    return aaa

    
def save_nbal(gy, pe, N=None, P=None, K=None, fold=None):
    if fold==None:
        fold="C:\Apps\WinPython-64bit-2.7.10.3\IPEWG\\"    
    nb=pd.DataFrame(np.zeros((len(N.totalDemand),8)), columns = ['Nsupply', 'Ndemand', 'Psupply', 'Pdemand', 'Ksupply', 'Kdemand','Subs', 'Vol'])
    if N != None: 
        nb['Nsupply']=N.totalSupply; nb['Ndemand']=N.totalDemand; 
    if P != None:
        nb['Psupply']=P.totalSupply; nb['Pdemand']=P.totalDemand
    if K != None:
        nb['Ksupply']=K.totalSupply; nb['Kdemand']=K.totalDemand
    ul = int(gy.g['DiamUpperLimit'][0]); nsteps = len(gy.age)  
    Vdist = np.reshape(gy.Vdist, (nsteps, ul))            
    nb['Vol'] = np.sum(Vdist, axis=1)
    nb['Subs']= pe.subsidence
    nb.to_csv(fold+'NutBal.csv')


def writeToPickle(pe, gy, pg, bal, Mrate, N, P, K, fname=None, measS=None):
    if fname == None:    
        fname = "/home/vagrant/.spyder2/Plantation/MonteCarlo3.p"
    data = {'simNo':[], 'annualSubsidence':[], 'meanSubs':[], 'subsidence':[],'rhoIni':[], 'rhoFin':[], 'siteIndex':[],
            'Mrate': [], 'oxidShare':[], 'dwt':[], 'bal':[], 'measS':[], 
            'Nsupply':[], 'Psupply':[], 'Ksupply':[],'Ndemand':[], 'Pdemand':[], 'Kdemand':[],
            'Ntreesupply':[], 'Ptreesupply':[], 'Ktreesupply':[],'Ntreedemand':[], 'Ptreedemand':[], 'Ktreedemand':[],
            }  
    if os.path.exists(fname):
        with open(fname,'rb') as rfp: 
            data = pickle.load(rfp)
    else:
        fi = open(fname, 'w') 
        pickle.dump(data, fi)
        fi.close()
    data['simNo'].append(len(data['siteIndex']))
    data['annualSubsidence'].append(pe.annualSubsidence) 
    data['meanSubs'].append(np.mean(pe.annualSubsidence)) 
    data['subsidence'].append(pe.subsidence)
    data['rhoIni'].append(pe.peatRhoInit) 
    data['rhoFin'].append(pe.peatRhoFinal) 
    data['siteIndex'].append(gy.appSI)       
    data['Mrate'].append(Mrate)    
    data['oxidShare'].append(pe.oxidShare)       
    data['dwt'].append(pe.dwt)  
    data['bal'].append(bal)       
    data['measS'].append(measS)
    data['Ndemand'].append(N.totalDemand)
    data['Pdemand'].append(P.totalDemand)
    data['Kdemand'].append(K.totalDemand)
    data['Nsupply'].append(N.totalSupply)
    data['Psupply'].append(P.totalSupply)
    data['Ksupply'].append(K.totalSupply)
    data['Ntreedemand'].append(N.treeDemand)
    data['Ptreedemand'].append(P.treeDemand)
    data['Ktreedemand'].append(K.treeDemand)
    data['Ntreesupply'].append(N.treeSupply)
    data['Ptreesupply'].append(P.treeSupply)
    data['Ktreesupply'].append(K.treeSupply)
    
    evans_subs =    np.mean(pe.dwt)*-0.0334-1.88 
    with open(fname,'wb') as wfp:
        pickle.dump(data, wfp)
    print (np.round(np.mean(pe.annualSubsidence),2), np.mean(measS), 
           np.round(np.mean(pe.dwt),2), np.round(evans_subs,2))

def writeToPickle2(pe, gy, pg, bal, Mrate, N, P, K, fname=None, measS=None):
    if fname == None:    
        fname = "/home/vagrant/.spyder2/Plantation/MonteCarlo3.p"
    data = {'simNo':[], 'annualSubsidence':[], 'meanSubs':[], 'subsidence':[],'rhoIni':[], 'rhoFin':[], 'siteIndex':[],
            'Mrate': [], 'oxidShare':[], 'dwt':[], 'bal':[], 'measS':[], 
            'Nsupply':[], 'Psupply':[], 'Ksupply':[],'Ndemand':[], 'Pdemand':[], 'Kdemand':[],
            'Ntreesupply':[], 'Ptreesupply':[], 'Ktreesupply':[],'Ntreedemand':[], 'Ptreedemand':[], 'Ktreedemand':[],
            'standLitter':[], 'weedLitter':[], 'CWD':[]}  
    if os.path.exists(fname):
        with open(fname,'rb') as rfp: 
            data = pickle.load(rfp)
    else:
        fi = open(fname, 'wb') 
        pickle.dump(data, fi)
        fi.close()
    data['simNo'].append(len(data['siteIndex']))
    data['annualSubsidence'].append(pe.annualSubsidence) 
    data['meanSubs'].append(np.mean(pe.annualSubsidence)) 
    data['subsidence'].append(pe.subsidence)
    data['rhoIni'].append(pe.peatRhoInit) 
    data['rhoFin'].append(pe.peatRhoFinal) 
    data['siteIndex'].append(gy.appSI)       
    data['Mrate'].append(Mrate)    
    data['oxidShare'].append(pe.oxidShare)       
    data['dwt'].append(pe.dwt)  
    data['bal'].append(bal)       
    data['measS'].append(measS)
    data['Ndemand'].append(N.totalDemand)
    data['Pdemand'].append(P.totalDemand)
    data['Kdemand'].append(K.totalDemand)
    data['Nsupply'].append(N.totalSupply)
    data['Psupply'].append(P.totalSupply)
    data['Ksupply'].append(K.totalSupply)
    data['Ntreedemand'].append(N.treeDemand)
    data['Ptreedemand'].append(P.treeDemand)
    data['Ktreedemand'].append(K.treeDemand)
    data['Ntreesupply'].append(N.treeSupply)
    data['Ptreesupply'].append(P.treeSupply)
    data['Ktreesupply'].append(K.treeSupply)

    standlitter = gy.FineRLitL + gy.RootLitD + gy.BrLitD + gy.BrLitL + gy.BaLitD + gy.BaLitL + gy.FoLitD + gy.FoLitL                                                # foliage litter from living trees, kg ha-1 timestep-1
    data['standLitter'].append(standlitter)
    weedlitter = gy.weedALitL + gy.weedBLitL + gy.weedALitD + gy.weedBLitD                                            # below ground weed litter from dead weeds, kg ha-1 timestep-1         
    data['weedLitter'].append(weedlitter)
    data['CWD'].append(gy.CWD)
    
    with open(fname,'wb') as wfp:
        pickle.dump(data, wfp)
    evans_subs =    np.mean(pe.dwt)*-0.0334-1.88 
    print (np.round(np.mean(pe.annualSubsidence),2), np.mean(measS), 
        np.round(np.mean(pe.dwt),2), np.round(evans_subs,2))
    
def readFromPickle(fname=None):
    if fname == None:    
        fname = "/home/vagrant/.spyder2/Plantation/MonteCarlo4.p"
   # Re-load our database
    with open(fname,'rb') as rfp:
        data = pickle.load(rfp)
    #print data
    #for p in range(len(data['siteIndex'])):
    #    print data['siteIndex'][p]
    print (np.mean(data['meanSubs']), np.std(data['meanSubs']))
    #print(data[0][1]['annualSubsidence'])
    
def spatial_out(gy,de,ro,pe,P):
    """
    as time series:
        -subsidence
        -potential vol growth
        -c balance (cumulative)
        -p balance (instantanous)
    """
    s=(ro.SH + ro.H + ro.LH + ro.F + ro.L)*10000.0
    soil = (de.cwdSto + (ro.SH + ro.H + ro.LH + ro.F + ro.L)*10000.0 + pe.peatStorage)*ro.massToC 
    standTot = (gy.BiBark + gy.BiBranch  + gy.BiFoliage + gy.RootMass + gy.weedAbove + gy.weedBelow)*ro.massToC 
    Cbal = (soil -soil[0]) + (standTot - standTot[0])
    Pbal = P.totalSupply-P.totalDemand
    Ps=P.totalSupply
    s= pe.subsidence; v=gy.V; mv=gy.MV; sur=gy.Survival
    return s,v, mv, sur, Cbal, Pbal, Ps

def saveSubsidence(pe, datesInp=None, start= None, nsteps= None, midDwtArr= None, sideDwtArr=None, fold = 'C:\Apps\WinPython-64bit-2.7.10.3\IPEWG\\', fname= 'subs.p'):
    data={}
    if datesInp is not None: data['dates']=datesInp
    if datesInp is None:
        end = start + relativedelta(months=+nsteps-1)
        dr=pd.date_range(start,end)
        dfD=pd.DataFrame(range(len(dr)),index=dr); dfD=dfD.resample('M', how='mean')
        datesInp=dfD.index.values; data['dates']=datesInp
    data['subs']=pe.subsidence
    if midDwtArr is not None: data['midDwtArr'] = midDwtArr
    if sideDwtArr is not None: data['sideDwtArr'] = sideDwtArr
    with open(fold+fname,'wb') as wfp:
        pickle.dump(data, wfp)
