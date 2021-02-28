# -*- coding: utf-8 -*-
"""
Created on Mon Aug 07 19:13:54 2017

@author: lauren
"""

import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
from scipy.interpolate import interp1d
import datetime

from hydro_utils import getParams, getRainfall, CWTr, gmeanTr, Hadjacent, Amatrix
from hydro_utils import boundConst, rightSide, boundNoFlow

def run_striphy(hdom, LAI, sim_yrs, optFig=False):
    print ('**************************')    
    print ('STRIPHY hydrological model')
    print ('**************************')     
    #********** Stand parameters and weather forcing*******************
    #Changeable parameters, later to input....
    wpara ={
        'riau': {
        'infolder': 'C:\\Users\\alauren\\Documents\\WinPython-64bit-2.7.10.3\\IPEWG\\Plantation2.0\\',
        'infile':'rainfall.csv',
        'use_yr': 2012, 
        'description': 'Pekanbaru, Riau'},
        }
    spara = {
        'nLyrs':100, 'dzLyr': 0.05, 'L': 150., 'n':50, 'ditch depth': -0.4,
        'initial h': -0.2, 'slope': 0.02, 'peat type':'TPmean', 
        'infolder':'C:\\Users\\alauren\\Documents\\WinPython-64bit-2.7.10.3\\IPEWG\\Plantation2.0\\',
        'infile':'DomeData.xlsx'
        }

    #****** Constants, defined here **************    
    wlocation = 'riau'
    use_yr = wpara[wlocation]['use_yr']
    soilpara = getParams(ParamFile=spara['infolder'] + spara['infile'], shname=spara['peat type'], lrs= spara['nLyrs'])                                                             #soil layer hydraulic parameters 
    rain = getRainfall(rainFile=wpara[wlocation]['infolder']+wpara[wlocation]['infile'])                              #rainfall data
    hToSto, stoToh, hToTra  = CWTr(soilpara['profPara'], soilpara['pF'], direction='negative')  #storage, head and transmissivity tables          
    gwl = list(np.linspace(0.0, -5.0, 150))                                                     
    tmp=list(np.diff(hToSto(gwl))/np.diff(gwl))
    tmp.reverse(), gwl.reverse()
    C = interp1d(np.array(gwl), np.array(np.gradient(hToSto(gwl))/np.gradient(gwl)), fill_value='extrapolate')  #storage coefficient function      

    use_yr = wpara[wlocation]['use_yr']    
    start_date = datetime.datetime(use_yr,1,1); end_date=datetime.datetime(use_yr+1,1,1)
    length = (end_date - start_date).days*sim_yrs 
    deltas = np.zeros(length)#+sim_yrs)                                           # Infliltration-evapotranspiration    
    ets = np.zeros(length)#+sim_yrs)

    if sim_yrs > 1:    
        laif = interp1d(np.linspace(0, length, len(LAI)),LAI) 
        ets = 3. + 0.1*laif(range(length))                                                                # constant evapotranspoiration mm/day    
    else:
        laif = LAI[0]*np.ones(length)
        ets = 3. + 0.1*laif                                                                # constant evapotranspoiration mm/day    

    rain = rain[start_date:datetime.datetime(use_yr,12,31)]
    rain = pd.concat([rain]*sim_yrs)    
    rain=np.ravel(rain.values)
    deltas = rain - ets
    print ('Weather input: ', wpara[wlocation]['description'], ', year ', 
           use_yr, 'repeated for ', sim_yrs, ' years' )
 
    print ('Total infiltration, mm ', np.round(np.cumsum(deltas)[-1],2))
    print ('Total precipitation, mm ', np.round(np.cumsum(rain)[-1]) )    

    #******** Soil and strip parameterization *************************
    nLyrs = spara['nLyrs']                                                     # number of soil layers
    dz = np.ones(nLyrs)*spara['dzLyr']                                          # thickness of layers, m
    z = np.cumsum(dz)-dz/2.                                                     # depth of the layer center point, m 
    L= spara['L']                                                               # compartemnt width, m
    n= spara['n']                                                               # number of computation nodes
    dy = float(L/n)                                                             # node width m
    h0 =  spara['ditch depth']                                                  # Initial conditions, gw at canal
    hini = spara['initial h']                                                   # h in the compartment
    sl= spara['slope']                                                          # slope %
    lev=1.                                                                  # basic level of soil surface
    DrIrr=False
    mound_h = 0.0                                                           # mound height m
    mound_dist = 2                                                          # mound distance from each other m
    mounds = np.zeros(n)
    mounds = np.where(np.cumsum(np.ones(n)*dy)%mound_dist < 0.2, mound_h, 0.0)
    mounds[-1]=0.; mounds[0]=0.
    ele = np.linspace(0,L*sl/100., n) + lev                                 # surface rise in y direction, m
    ele = ele + mounds
    implic =  1.                                                            # 0-forward Euler, 1-backward Euler, 0.5-Crank-Nicolson
    dt= 1                                                                   # time step, days

    #**********Soil computation matrix************
    print ('Computing soil water fluxes ', len(deltas), ' days:')

    A=np.zeros((n,n))                                                       # computation matrix 
    h = np.ones(n)*hini                                                     # right hand side vector
    h[0]=h[n-1]= h0                                                         # symmetrical boundaries, set water level
    H=ele+h                                                                 # head with respect to absolute reference level, m
    hts = np.empty((int(length/dt+sim_yrs),n), dtype=float)                         # result table ndarray(days, number of nodes)
    sruno=0.
    for d, p in enumerate(deltas):
        Htmp = H.copy(); Htmp1=H.copy()
        S = p/1000.*np.ones(n)                                              # source/sink, in m
        airv = hToSto(ele)-hToSto(Htmp-ele)                                 # air volume, in m
        S=np.where(S>airv, airv, S)                        
        sruno += np.where(np.ones(len(airv))*(p)/1000. > airv, np.ones(len(airv))*(p)/1000.-airv, 0.0)  #cut the surface water above to runoff
        Tr0 = hToTra(H-ele)                                                 # Transmissivity from the previous time step 
        Trminus0, Trplus0 = gmeanTr(Tr0)                                    # geometric mean of adjacent node transmissivities
        Hminus, Hplus = Hadjacent(H)                                        # vector of adjacent node H

        for it in range(100):                                               # iteration loop for implicit solution
                                                                            # Htmp from old iteration, Htmp1 new iteration 
            Tr1 = hToTra(Htmp1-ele)                                         # transmissivity in new iteration
            CC = C(Htmp1-ele)                                               # storage coefficient in new iteration
            Trminus1, Trplus1 = gmeanTr(Tr1)                                # geometric mean of adjacent node transmissivity                
            alfa = CC*dy**2/dt
            A = Amatrix(A, n, implic, Trminus1, Trplus1, alfa)              # construct tridiaginal matrix
            A = boundConst(A, n)                                             # constant head boundaries to A matrix
            hs=rightSide(S,dt,dy,implic,alfa, H,Trminus0,Hminus,Trplus0,Hplus,\
                DrIrr, Htmp1, ele, h0)                                          #right hand side of the equation
            Htmp1 = np.linalg.multi_dot([np.linalg.inv(A),hs])              # solve equation                    
            Htmp1=np.where(Htmp1>ele, ele,Htmp1)                            # cut the surface water 
            conv = max(np.abs(Htmp1-Htmp))                                  # define convergence
            Htmp=Htmp1.copy()                                               # new wt to old for new iteration
            if conv < 1.e-7: 
                if d%365==0: print ('  - day #',d, 'iterations', it)                    
                break
        H=Htmp1.copy()                     
        hts[d,:]=H
        
    gwlev =[np.mean(np.abs(ele - hts[d,:])[1:-1]) for d in range(len(deltas))]
    print ('mean gwl', np.round(np.mean(gwlev), 3), ' m')
    print ('Canal water level ', h0, ' m')
    print ('Field drain spacing  ', L, ' m')    
    dfOut = pd.DataFrame(data={'gwl':gwlev}, index=pd.date_range(start_date,periods=length))  #len(deltas)
    dfOut = dfOut[str(use_yr):str(use_yr+sim_yrs-1)]
    dfOut=dfOut.resample('M', how='mean')
    diffL = len(LAI) - len(dfOut)    
    if diffL > 0:
        for d in range(diffL):
            arr= dfOut['gwl'][-1]
            df2 = pd.DataFrame([arr], columns=['gwl'])
            dfOut=dfOut.append(df2, ignore_index=True)
    if optFig: fig_hydro(ele, dfOut, hts, spara, wpara, wlocation, ets, rain,  LAI, sim_yrs)
    import pickle
    data={}
    data['gwl']= [np.abs(ele - hts[d,:]) for d in range(len(deltas))]   
    fname = "gwls_" + str(np.abs(h0))+"_results.p"
    with open(wpara[wlocation]['infolder']+fname, 'wb') as wfp:
        pickle.dump(data, wfp)

    return dfOut


def fig_hydro(ele, dfOut, hts, spara, wpara, wlocation, ets, Prec, lai, sim_yrs):
    from matplotlib.lines import Line2D
    h0=spara['ditch depth']; n=spara['n']; L=spara['L']

    aa, bb = np.shape(hts); x = np.linspace(0,L,n); dy = float(L/n) 
    fig= plt.figure(num = 'Striphy', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(20.0,11.0))   #Figsize(w,h), tuple inches 
    ax = fig.add_axes([0.05, 0.5, 0.55, 0.46]) #left, bottom, width, height
    low = min([ele[0]+h0, ele[n-1]+h0])*0.4; high = max(ele)*1.2
    ax.set_ylim([low,high])
    line2, = ax.plot(x[1:n-1], ele[1:n-1], 'k-', linewidth = 2, label = 'Surface elevation')
    line1 = Line2D([], [], color='blue', marker='o', markeredgecolor='b', markersize = 1, label='Water table')
    limit = ele-0.4
    limitGr = ele-0.2
    plt.plot(x, ele , 'k-', x, limit, 'k--')
    mid=int(n/2); west=int(n/10); east =n-int(n/10)
    plt.fill_between(x, limit, ele, color='green', alpha= 0.3)       
    plt.fill_between(x, limitGr, ele, color='red', alpha=0.3) 
    plt.plot(x[mid], ele[mid], 'ro', markersize=20)
    plt.plot(x[west], ele[west], 'go', markersize=20)
    plt.plot(x[east], ele[east], 'mo', markersize=20)
    plt.xlabel('Distance, m', fontsize=16); plt.ylabel('Elevation, m', fontsize=16)
    
    y = np.mean(hts, axis=0)
    yu=y+np.std(hts, axis=0); yl=y-np.std(hts, axis=0)
    yuu=np.max(hts, axis=0); yll=np.min(hts, axis=0)
    line1.set_data(x, y)
    plt.plot(x, yu, 'b--')
    plt.plot(x, yl, 'b--')
    plt.fill_between(x, yl, yu, color='blue', alpha=0.1)
    plt.plot(x, yuu, linestyle='dotted')
    plt.plot(x, yll, linestyle='dotted')
    plt.fill_between(x, yll, yuu, color='blue', alpha=0.1)
    
    t1='Ditch depth ' + str(h0) + ' m'
    ax.text(0.5, high*0.95 , str(t1),fontsize=16, color='0.25')
    t2 = 'Slope ' + str(spara['slope']) + ' %'
    ax.text(0.5, high*0.9 , str(t2),fontsize=16, color='0.25')
    t3 = 'Peat type: ' + spara['peat type']
    ax.text(0.5, high*0.85 , t3,fontsize=16, color='0.25')
    if type(lai) is float: lai=[lai]    
    t4 = 'LAI ' + str(lai[0]) + '...' + str(lai[-1])    
    ax.text(0.5, high*0.8 , t4,fontsize=16, color='0.25')
       
    ax.add_line(line1)
    ax.legend(loc=1)
    
    ax2=fig.add_axes([0.05, 0.1, 0.55, 0.3]) #left, bottom, width, height)
    ax2.set_ylim([h0*2., 0.2])
    surf = np.zeros(aa)
    limit2=-0.4*np.ones(aa); limit3=-0.2*np.ones(aa)
    plt.plot(range(aa), limit2, 'k--', range(aa), surf, 'k-')
    plt.fill_between(range(aa), limit2, surf, color='green', alpha=0.3)
    plt.fill_between(range(aa), limit3, surf, color='red', alpha=0.3)
    xx = range(np.shape(hts)[0])    
    yy = hts[:,mid] - ele[mid]
    yywest =hts[:,west] - ele[west]    
    yyeast =hts[:,east] - ele[east]    
    line3 = Line2D([], [], color='red', marker='o', markeredgecolor='b', markersize = 1, linewidth=1, label='Mid field wt')
    ax2.add_line(line3)
    line4 = Line2D([], [], color='green', marker='o', markeredgecolor='b', markersize = 1, linewidth=0.5, label='West wt')
    ax2.add_line(line4)
    line5 = Line2D([], [], color='m', marker='o', markeredgecolor='b', markersize = 1, linewidth=0.5, label='East wt')
    ax2.add_line(line5)
    ax2.set_xlim([0,aa])
    plt.xlabel('Time, days', fontsize = 16); plt.ylabel('wt depth, m', fontsize = 16)
    ax2.legend(loc=1)
    line3.set_data(range(len(yy)),yy)
    line4.set_data(range(len(yywest)),yywest)
    line5.set_data(range(len(yyeast)),yyeast)

    ax3 = fig.add_axes([0.68, 0.75, 0.3, 0.21]) #left, bottom, width, height
    ax3.set_ylim([0, sum(Prec)*1.1]); ax3.set_xlim([0,aa])
    line6 = Line2D([], [], color='c', marker='o', markeredgecolor='b', markersize = 1, label='Cumul P')
    ax3.add_line(line6)
    line7 = Line2D([], [], color='k', marker='o', markeredgecolor='b', markersize = 1, label='Cumul ET')
    ax3.add_line(line7)
    plt.xlabel('Time, days', fontsize=16); plt.ylabel('mm', fontsize=16)
    ax3.legend(loc=1)
    line6.set_data(range(len(Prec)),np.cumsum(Prec))
    line7.set_data(range(len(ets)), np.cumsum(ets))

    t1='Weather from '+ wpara[wlocation]['description']
    ax3.text(10, sum(Prec)*1.0 , str(t1),fontsize=14, color='0.25')
    t2='Year ' + str(wpara[wlocation]['use_yr']) + ' repeated ' +str(sim_yrs) + ' times'
    ax3.text(10, sum(Prec)*0.85 , t2,fontsize=14, color='0.25')
    t3 = 'Rain '+ str(np.round(sum(Prec)/sim_yrs)) + ', ET ' +str(np.round(sum(ets)/sim_yrs)) + ' mm yr-1'
    ax3.text(10, sum(Prec)*0.7 , t3,fontsize=14, color='0.25')
    """    
    ax4 = fig.add_axes([0.68, 0.45, 0.3, 0.21]) #left, bottom, width, height
    ax4.set_ylim([-10, 35]); ax4.set_xlim([0,aa])
    line8 = Line2D([], [], color='b', marker='o', markeredgecolor='b', markersize = 1, label='Air temperature, deg C')
    ax4.add_line(line8)
    nolla=np.zeros(aa); lo = np.ones(aa)*-10
    plt.fill_between(range(aa), lo, nolla, color='yellow', alpha=0.2)
    line8.set_data(range(len(T)), T)
    plt.xlabel('Time, days', fontsize = 16); plt.ylabel('deg C', fontsize = 16)
    ax4.legend(loc=1)
    t1 = 'Mean temperature ' + str(np.round(np.mean(T)))    
    ax4.text(10, 30 , t1, fontsize=14, color='0.25')
    
    ax5 = fig.add_axes([0.68, 0.1, 0.3, 0.21]) #left, bottom, width, height
    line9 = Line2D([], [], color='b', marker='o', markeredgecolor='b', markersize = 1, label='CO2 efflux')
    ax5.set_ylim([-10, 100]); ax5.set_xlim([0,aa])
    ax5.add_line(line9)
    ax5.legend(loc=1)
    plt.fill_between(range(aa), lo, nolla, color='yellow', alpha=0.2)
    line9.set_data(range(len(het)), het)
    plt.xlabel('Time, days', fontsize = 16); plt.ylabel('kg/ha', fontsize = 16)
    t1 = 'Annual CO2 efflux ' + str(np.round(np.sum(het)/sim_yrs))
    ax5.text(10, 85 , t1, fontsize=14, color='0.25')
    """
#hdom=[2.,10.,15.]; lai=[1.,10.,12.]; sim_yrs=3
#hdom=[2.]; lai=[1.]; sim_yrs=1
#dfOut = run_striphy(hdom, lai, sim_yrs, optFig=True) 
#print dfOut