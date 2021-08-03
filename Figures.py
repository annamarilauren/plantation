# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 17:10:20 2016

@author: vagrant
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline as spline
from scipy.interpolate import interp1d
fcol =  (232/255.0, 243/255.0, 245.0/255); figsi = (20.0,9.0)       
import seaborn as sns
sns.set()
def fDrawGrowth(gy):
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)    
    dt = gy.g['dt'][0]; dty = dt/12.0
    time = np.arange(0, len(gy.Hdom))*dty
    fig= plt.figure(num = 'Plantation simulator, stand', facecolor=fcol, edgecolor='k',figsize=figsi)   #Figsize(w,h), tuple inches 
    ax1 = fig.add_axes([0.33, 0.47, 0.3, 0.46]) #left, bottom, width, height
    ax1.set_title('Stand characteristics')    
    ax1.set_xlabel('Age, yrs', fontsize = 14)
    ax1.plot(time, gy.Hdom, 'r-', label='$H_{dom}$, m')
    ax1.plot(time, gy.BA, 'b-', label='BA, $m^{2}$')
    ax1.plot(time, gy.WMeanDiam, 'c-', label='$\overline{dbh}$, cm')
    ax1.legend(loc='upper left')
    ax11=ax1.twinx()
    ax11.plot(time, gy.Survival, 'g-', label='Stocking')
    ax11.legend(loc='lower right')
    ax11.set_ylim(0,max(gy.Survival)*1.2)

    ax2 = fig.add_axes([0.33, 0.05, 0.3, 0.33]) #left, bottom, width, height
    ax2.set_xlabel('Diameter, cm', fontsize=14)
    ax2.set_title('Dbh distribution, stems $ha^{-1}$')        
    ul = int(gy.g['DiamUpperLimit'][0]); nsteps = len(gy.age)  
    DbhDist = np.reshape(gy.DbhDist, (nsteps, ul))
    diams = np.arange(0,(gy.g['DiamUpperLimit'][0]), dtype='int')    

    xnew = np.linspace(0,40,300)
    length = DbhDist.shape[0]; dt = gy.g['dt'][0]    
    ages = np.linspace(0, length * dt, 6); lables=[str(x/12.0)[0:3] for x in ages]
    age1 = int(ages[1]); age2 =int(ages[1]); age3 =int(ages[3]); 
    age4=int(ages[4]); age5 =int(ages[5])-1

    cs = spline(diams, DbhDist[age1][0:])
    ax2.plot(xnew, cs(xnew), 'r-', label=lables[1] + ' yr')
    cs = spline(diams, DbhDist[age2][0:])
    ax2.plot(xnew, cs(xnew), 'b-', label=lables[2] + ' yr')
    cs = spline(diams, DbhDist[age3][0:])
    ax2.plot(xnew, cs(xnew), 'y-', label=lables[3] + ' yr')
    cs = spline(diams, DbhDist[age4][0:])
    ax2.plot(xnew, cs(xnew), 'g-', label=lables[4] + ' yr')
    cs = spline(diams, DbhDist[age5][0:])
    ax2.plot(xnew, cs(xnew), 'k-', label=lables[5] + ' yr')
    ax2.legend(loc='upper right')    
    ax2.set_xlim(0,40)

    Vdist = np.reshape(gy.Vdist, (nsteps, ul)); MercDist = np.reshape(gy.MercDist, (nsteps, ul))
    ax3 = fig.add_axes([0.68, 0.40, 0.3, 0.25]) #left, bottom, width, height
    ax3.set_xlabel('Age, yrs', fontsize=14)
    ax3.set_title('Volume development, $m^{3}$ $ha^{-1}$ ')
    ax3.plot(time, np.sum(Vdist, axis=1), 'r-', label='Over bark')
    ax3.plot(time, np.sum(MercDist, axis=1), 'b-', label='Merchantable')
    ax3.fill_between(time,0,np.sum(MercDist, axis=1),facecolor='green', alpha=0.3 )
    ax3.set_xlim(0,time[-1])
    ax3.legend(loc=2)

    ax4 = fig.add_axes([0.68, 0.05, 0.3, 0.25]) #left, bottom, width, height
    ax4.set_xlabel('Age, yrs', fontsize=14)
    ax4.set_title('Biomass, kg $ha^{-1}$ ')
    y = gy.MerchMass
    ax4.plot(time, y , 'g-', label='Stem')
    ax4.fill_between(time,0, y,facecolor='green', alpha=0.3)
    yy = y + gy.BiBark
    ax4.plot(time, yy , 'b-', label='Bark')
    ax4.fill_between(time,y,yy,facecolor='blue', alpha=0.3)
    y = yy + gy.BiBranch
    ax4.plot(time, y , 'r-', label='Branch')
    ax4.fill_between(time,yy,y,facecolor='red', alpha=0.3)
    yy = y + gy.BiFoliage
    ax4.plot(time, yy , 'c-', label='Foliage')
    ax4.fill_between(time,y,yy,facecolor='cyan', alpha=0.3)
    y = yy + gy.weedAbove
    ax4.plot(time, y, 'y-', label = 'Weed')
    ax4.fill_between(time, yy, y, facecolor = 'yellow', alpha=0.3)
    ax4.plot(time, gy.RootMass*-1.0 , 'k-', label='Root')
    ax4.fill_between(time, gy.RootMass*-1.0 ,0,facecolor='gray', alpha=0.3)
    y = -gy.RootMass - gy.weedBelow
    ax4.plot(time, y, 'm-', label='Weed roots')
    ax4.fill_between(time, gy.RootMass*-1.0, y, facecolor='magenta', alpha=0.3)
    ax4.legend(loc='upper left')
    ax4.set_xlim(0,time[-1])

    ax61= fig.add_axes([0.03, 0.74, 0.12, 0.18]) #left, bottom, width, height
    a=gy.BaLitD; b = gy.BaLitL 
    ax61.plot(time, a, '0.75')
    ax61.fill_between(time, 0, a, facecolor='gray', alpha=0.3)
    ax61.plot(time, a+b, 'g')   
    ax61.fill_between(time, a, a+b, facecolor='green', alpha=0.3)
    ax61.set_title('Bark litter')
    ax61.set_xlim(-1, time[-1])    
    ax61.set_ylim(bottom=0)    
    
    ax62= fig.add_axes([0.17, 0.74, 0.12, 0.18]) #left, bottom, width, height
    a=gy.BrLitD; b = gy.BrLitL 
    ax62.plot(time, a, '0.75')
    ax62.fill_between(time, 0, a, facecolor='gray', alpha=0.3)
    ax62.plot(time, a+b, 'g')   
    ax62.fill_between(time, a, a+b, facecolor='green', alpha=0.3)
    ax62.set_title('Branch litter')
    ax62.set_xlim(-1, time[-1])    
    ax62.set_ylim(bottom=0)    
    
    ax63= fig.add_axes([0.03, 0.49, 0.12, 0.18]) #left, bottom, width, height
    a=gy.FoLitD; b = gy.FoLitL; c = gy.weedALitD; d = gy.weedALitL  
    ax63.plot(time, a, '0.75')
    ax63.fill_between(time, 0, a, facecolor='gray', alpha=0.3)
    ax63.plot(time, a+c, '0.25')
    ax63.fill_between(time, a, a+c, facecolor = 'gray', alpha = 0.75)
    ax63.plot(time, a+c+b, 'g')   
    ax63.fill_between(time, a+c, a+c+b, facecolor='green', alpha=0.3)
    ax63.plot(time, a+b+c+d, 'g')
    ax63.fill_between(time, a+b+c, a+b+c+d, facecolor='green', alpha=0.75)
    ax63.set_title('Foliage litter')
    ax63.set_xlim(-1, time[-1])    
    ax63.set_ylim(bottom=0)    

    ax64= fig.add_axes([0.17, 0.49, 0.12, 0.18]) #left, bottom, width, height
    a=gy.RootLitD; b = gy.FineRLitL; c = gy.weedBLitD; d = gy.weedBLitL   
    ax64.plot(time, a, '0.75')
    ax64.fill_between(time, 0, a, facecolor='gray', alpha=0.3)
    ax64.fill_between(time, 0, a, facecolor='gray', alpha=0.3)
    ax64.plot(time, a+c, '0.25')
    ax64.plot(time, a+b+c, 'g')   
    ax64.fill_between(time, a+c, a+c+b, facecolor='green', alpha=0.3)
    ax64.plot(time, a+b+c+d, 'g')
    ax64.fill_between(time, a+b+c, a+b+c+d, facecolor='green', alpha=0.75)
    ax64.set_title('Root litter')
    ax64.set_xlim(-1, time[-1])    
    ax64.set_ylim(bottom=0)    

    ax65= fig.add_axes([0.03, 0.23, 0.12, 0.18]) #left, bottom, width, height
    a=gy.CWD 
    ax65.plot(time, a, '0.75')
    ax65.fill_between(time, 0, a, facecolor='gray', alpha=0.3)
    ax65.set_title('CWD')
    ax65.set_xlim(-1, time[-1])    
    ax65.set_ylim(bottom=0)    

    ax66= fig.add_axes([0.17, 0.23, 0.12, 0.18]) #left, bottom, width, height    
    a = gy.BaLitD + gy.BrLitD + gy.FoLitD + gy.RootLitD + gy.CWD + gy.weedALitD + gy.weedBLitD
    b = gy.BaLitL + gy.BrLitL + gy.FoLitL + gy.FineRLitL + gy.weedALitL + gy.weedBLitL
    ax66.plot(time, a, '0.75')
    ax66.fill_between(time, 0, a, facecolor='gray', alpha=0.3)
    ax66.plot(time, a+b, 'g')   
    ax66.fill_between(time, a, a+b, facecolor='green', alpha=0.3)
    ax66.set_title('Total')
    ax66.set_xlim(-1, time[-1])    
    ax66.set_ylim(bottom=0)    

    ax7= fig.add_axes([0.03, 0.03, 0.26, 0.15]) #left, bottom, width, height    
    a = np.cumsum(a); b= np.cumsum(b)
    ax7.plot(time, a, '0.75')
    ax7.fill_between(time, 0, a, facecolor='gray', alpha=0.3)
    ax7.plot(time, a+b, 'g')   
    ax7.fill_between(time, a, a+b, facecolor='green', alpha=0.3)
    ax7.set_title('Cumulative litter')
    ax7.set_xlim(0, time[-1])    
    ax7.set_ylim(bottom=0)    
    
    #----meatadata text: Table  ---------------  
    
    le = 0.675; w = 0.045; fs =14    
    row0 = 0.93; rh =0.03
    plt.figtext(le+0*w, row0, 'Rotat', fontsize = fs, fontweight='bold')
    plt.figtext(le+1.0*w, row0, 'Species', fontsize = fs, fontweight='bold')        
    iage = str(int(gy.g['Iage'][0]))    
    plt.figtext(le+2.0*w, row0, 'SI at '+ iage+ ' y', fontsize = fs, fontweight='bold')
    plt.figtext(le+3.0*w, row0, 'MercVol', fontsize = fs, fontweight='bold')        
    plt.figtext(le+4.0*w, row0, 'Tons', fontsize = fs, fontweight='bold')
    plt.figtext(le+5.0*w, row0, 'Biomass', fontsize = fs,fontweight='bold')        
    plt.figtext(le+6.0*w, row0, 'Litter', fontsize = fs,fontweight='bold')   
    Nrotat = int(gy.g['Nrotat'][0]); rotLength = int(len(gy.age)/Nrotat); lit = 0 
    for n in range(Nrotat):
        rg = gy.rowGe[n]  
        s = int(gy.g['BiomSpe'][rg])
        spe = str(gy.b['species'][s])[0:9]  
        plt.figtext(le+0*w, row0-(n+1)*rh, str(n+1), fontsize = fs)
        plt.figtext(le+1.0*w, row0-(n+1)*rh, spe, fontsize = fs) 
        plt.figtext(le+2.0*w, row0-(n+1)*rh, str(gy.g['SI'][rg]), fontsize = fs)        
        plt.figtext(le+3.0*w, row0-(n+1)*rh, str(np.sum(MercDist, axis=1)[(n+1)*rotLength-1])[0:5], fontsize = fs)        
        plt.figtext(le+4.0*w, row0-(n+1)*rh, str(gy.MerchMass[(n+1)*rotLength-1] / 1000.0)[0:4], fontsize = fs)
        t = gy.MerchMass[(n+1)*rotLength-1] + gy.BiBark[(n+1)*rotLength-1] + gy.BiBranch[(n+1)*rotLength-1] + \
            gy.BiFoliage[(n+1)*rotLength-1] + gy.RootMass[(n+1)*rotLength-1] + gy.weedAbove[(n+1)*rotLength-1] + \
            gy.weedBelow[(n+1)*rotLength-1]
        plt.figtext(le+5.0*w, row0-(n+1)*rh, str(t/1000.0)[0:5], fontsize = fs)        
        lit = (a[(n+1)*rotLength-2]+b[(n+1)*rotLength-2])/1000.0 - lit 
        plt.figtext(le+6.0*w, row0-(n+1)*rh, str(lit)[0:5], fontsize = fs)   
    plt.show()

def fDrawMAICAI(gy):
    params = {'mathtext.default': 'regular' }          
    fs = 14
    plt.rcParams.update(params)    
    dt = gy.g['dt'][0]; dty = dt/12.0
    time = np.arange(0, len(gy.Hdom))*dty
    fig= plt.figure(num = 'Plantation simulator, MAI and CAI', facecolor=fcol, edgecolor='k',figsize=(15.,6.))   #Figsize(w,h), tuple inches 
    ax1 = fig.add_axes([0.1, 0.15, 0.35, 0.7]) #left, bottom, width, height
    ax1.set_title('Stand characteristics')    
    ax1.set_xlabel('Age, years', fontsize = fs)
    ax1.plot(time, gy.Hdom, 'r-', label='$H_{dom}$, m')
    ax1.plot(time, gy.BA, 'b-', label='BA, $m^{2}$')
    ax1.plot(time, gy.WMeanDiam, 'c-', label='$\overline{dbh}$, cm')
    ax1.legend(loc='upper left', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)

    ax11=ax1.twinx()
    ax11.plot(time, gy.Survival, 'g-', label='Stocking')
    ax11.legend(loc='lower right', fontsize=fs)
    ax11.set_ylim(0,max(gy.Survival)*1.2)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
 
    ul = int(gy.g['DiamUpperLimit'][0]); nsteps = len(gy.age)  
    Vdist = np.reshape(gy.Vdist, (nsteps, ul)) #; MercDist = np.reshape(gy.MercDist, (nsteps, ul))
    volume = np.sum(Vdist, axis=1)
    ax3 = fig.add_axes([0.55, 0.15, 0.35, 0.7]) #left, bottom, width, height
    ax3.set_xlabel('Age, years', fontsize=fs)
    ax3.set_ylabel('Volume increment, $m^{3}$ $ha^{-1}$ $yr^{-1}$', fontsize=fs)
    ax3.plot(time,np.gradient(volume)*12. , 'r-', label='CAI')
    ax3.plot(time[1:], volume[1:]/time[1:] , 'b-', label='MAI')
    ax3.legend(loc='upper left', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    
   

def fAnimateGrowth(gy):
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)    
    
    xage = np.linspace(0, gy.age[-1], 300)    
    fig= plt.figure(num = 'Plantation simulator, stand', facecolor=fcol, edgecolor='k',figsize=figsi)   #Figsize(w,h), tuple inches 
    ax1 = fig.add_axes([0.33, 0.47, 0.3, 0.46]) #left, bottom, width, height
    ax1.set_title('Stand characteristics')    
    ax1.set_xlabel('Age, yrs', fontsize = 14)
    line01, = ax1.plot(gy.age, gy.Hdom, 'r-', label='$H_{dom}$, m')
    ax1.legend(loc='upper left')
    ax11=ax1.twinx()
    line11, = ax11.plot(gy.age, gy.Survival, 'g-', label='Stocking')
    ax11.legend(loc='lower right')
    ax11.set_ylim(0,max(gy.Survival)*1.2)

    ax2 = fig.add_axes([0.33, 0.05, 0.3, 0.33]) #left, bottom, width, height
    ax2.set_xlabel('Diameter, cm', fontsize=14)
    ax2.set_title('Dbh distribution, stems $ha^{-1}$')        
    diams = range(gy.g['DiamUpperLimit'][0])    
    xnew = np.linspace(0,40,300)
    xnew = np.linspace(0,int(gy.g['DiamUpperLimit'][0]),300)
    power_smooth1 = spline(diams, gy.DbhDist[6][0:],xnew)
    line21, = ax2.plot(xnew, power_smooth1, 'k-', linewidth=2)
    ax2.set_xlim(0,30)

    ax3 = fig.add_axes([0.68, 0.40, 0.3, 0.25]) #left, bottom, width, height
    ax3.set_xlabel('Age, yrs', fontsize=14)
    ax3.set_title('Volume development, $m^{3}$ $ha^{-1}$ ')
    line31, = ax3.plot(xage, spline(gy.age, np.sum(gy.Vdist, axis=1), xage), 'r-', label='Over bark')
    line32, = ax3.plot(xage, spline(gy.age, np.sum(gy.MercDist, axis=1), xage), 'b-', label='Merchantable')
    ax3.fill_between(gy.age,0,np.sum(gy.MercDist, axis=1),facecolor='green', alpha=0.3 )
    ax3.legend(loc=2)


    ax4 = fig.add_axes([0.68, 0.05, 0.3, 0.25]) #left, bottom, width, height
    ax4.set_xlabel('Age, yrs', fontsize=14)
    ax4.set_title('Biomass, kg $ha^{-1}$ ')
    y = gy.MerchMass
    line41, = ax4.plot(xage, spline(gy.age, y, xage) , 'g-', label='Stem')
    ax4.fill_between(gy.age,0, y,facecolor='green', alpha=0.3)
    yy = y + gy.BiBark
    line42, = ax4.plot(xage, spline(gy.age, yy, xage) , 'b-', label='Bark')
    ax4.fill_between(gy.age,y,yy,facecolor='blue', alpha=0.3)
    y = yy + gy.BiBranch
    line43, = ax4.plot(xage, spline(gy.age, y, xage) , 'r-', label='Branch')
    ax4.fill_between(gy.age,yy,y,facecolor='red', alpha=0.3)
    yy = y + gy.BiFoliage
    line44, = ax4.plot(xage, spline(gy.age, yy, xage) , 'c-', label='Foliage')
    ax4.fill_between(gy.age,y,yy,facecolor='cyan', alpha=0.3)
    line45, = ax4.plot(xage, spline(gy.age, gy.RootMass*-1.0, xage) , 'k-', label='Root')
    ax4.fill_between(gy.age,gy.RootMass*-1.0 ,0,facecolor='gray', alpha=0.3)
    ax4.legend(loc='upper left')

    
    ax61= fig.add_axes([0.03, 0.74, 0.12, 0.18]) #left, bottom, width, height
    a=gy.BaLitD; b = gy.BaLitL 
    line611, = ax61.plot(xage, spline(gy.age, a, xage), '0.75')
    ax61.fill_between(gy.age, 0, a, facecolor='gray', alpha=0.3)
    line612, = ax61.plot(xage, spline(gy.age, a+b, xage), 'g')   
    ax61.fill_between(gy.age, a, a+b, facecolor='green', alpha=0.3)
    ax61.set_title('Bark litter')
    
    ax62= fig.add_axes([0.17, 0.74, 0.12, 0.18]) #left, bottom, width, height
    a=gy.BrLitD; b = gy.BrLitL 
    line621, = ax62.plot(xage, spline(gy.age, a,xage), '0.75')
    ax62.fill_between(gy.age, 0, a, facecolor='gray', alpha=0.3)
    line622, = ax62.plot(xage, spline(gy.age, a+b, xage), 'g')   
    ax62.fill_between(gy.age, a, a+b, facecolor='green', alpha=0.3)
    ax62.set_title('Branch litter')

    ax63= fig.add_axes([0.03, 0.49, 0.12, 0.18]) #left, bottom, width, height
    a=gy.FoLitD; b = gy.FoLitL 
    line631, = ax63.plot(xage, spline(gy.age, a, xage), '0.75')
    ax63.fill_between(gy.age, 0, a, facecolor='gray', alpha=0.3)
    line632, = ax63.plot(xage, spline(gy.age, a+b, xage),'g')   
    ax63.fill_between(gy.age, a, a+b, facecolor='green', alpha=0.3)
    ax63.set_title('Foliage litter')

    ax64= fig.add_axes([0.17, 0.49, 0.12, 0.18]) #left, bottom, width, height
    a=gy.RootLitD; b = gy.FineRLitL 
    line641, = ax64.plot(xage, spline(gy.age, a, xage), '0.75')
    ax64.fill_between(gy.age, 0, a, facecolor='gray', alpha=0.3)
    line642, = ax64.plot(xage, spline(gy.age, a+b, xage), 'g')   
    ax64.fill_between(gy.age, a, a+b, facecolor='green', alpha=0.3)
    ax64.set_title('Root litter')

    ax65= fig.add_axes([0.03, 0.23, 0.12, 0.18]) #left, bottom, width, height
    a=gy.CWD 
    line651, = ax65.plot(xage, spline(gy.age, a, xage), '0.75')
    ax65.fill_between(gy.age, 0, a, facecolor='gray', alpha=0.3)
    ax65.set_title('CWD')

    ax66= fig.add_axes([0.17, 0.23, 0.12, 0.18]) #left, bottom, width, height    
    a = gy.BaLitD + gy.BrLitD + gy.FoLitD + gy.RootLitD + gy.CWD
    b = gy.BaLitL + gy.BrLitL + gy.FoLitL + gy.FineRLitL
    line661, = ax66.plot(xage, spline(gy.age, a, xage), '0.75')
    ax66.fill_between(gy.age, 0, a, facecolor='gray', alpha=0.3)
    line662, = ax66.plot(xage, spline(gy.age, a+b, xage),'g')   
    ax66.fill_between(gy.age, a, a+b, facecolor='green', alpha=0.3)
    ax66.set_title('Total')

    ax7= fig.add_axes([0.03, 0.03, 0.26, 0.15]) #left, bottom, width, height    
    a = np.cumsum(a); b= np.cumsum(b)
    line71, = ax7.plot(xage, spline(gy.age, a, xage),'0.75')
    ax7.fill_between(gy.age, 0, a, facecolor='gray', alpha=0.3)
    line72, = ax7.plot(xage, spline(gy.age, a+b, xage),'g')   
    ax7.fill_between(gy.age, a, a+b, facecolor='green', alpha=0.3)
    ax7.set_title('Cumulative litter')
    
    #----meatadata text  ---------------  
    s = 'Species: ' + str(gy.species)
    plt.figtext(0.68, 0.9, s, fontsize = 14)
    s = 'Site index: ' + str(gy.g['SI'][0]) + ' m at age of ' + str(gy.g['Iage'][0]) + ' yrs'
    plt.figtext(0.68, 0.87, s, fontsize = 14)

    
    s = 'Litter production: living trees (green), dead (gray) (kg $ha^{-1} $ $%d months^{-1}$ )' %gy.g['dt'][0]
    plt.figtext(0.02, 0.96, s, fontsize = 14)

    #New smoothed diameter distribution
    DD=np.resize(np.zeros(300*40),(300,40))
    for d in range(40):
        dd = interp1d(gy.age, gy.DbhDist[0:, d])
        for a in  range(300):
            DD[a][d] = dd(a/300.0*gy.age[-1])
      
    OptAni=True
    dur = 7; fpers =10
    if OptAni == True:
        # ANIMATE WITH MOVIEPY (UPDATE THE CURVE FOR EACH t). MAKE A GIF.
        import moviepy.editor as mpy  
        from moviepy.video.io.bindings import mplfig_to_npimage
        from moviepy.video.fx.resize import resize
        def make_frame_mpl(t):
            i = int(t/7.0*float(len(xage)))
            line01.set_ydata(spline(gy.age, gy.Hdom, xage)[0:i])  
            line11.set_ydata(spline(gy.age, gy.Survival, xage)[0:i])  
            line01.set_xdata(xage[0:i])        
            line11.set_xdata(xage[0:i])         
            if xage[i]>0.6:            
                power_smooth1 = spline(diams, DD[i,0:],xnew)
                line21.set_ydata(power_smooth1)
            else:
                line21.set_ydata(np.zeros(len(xnew)))
            line31.set_ydata(spline(gy.age, np.sum(gy.Vdist, axis=1), xage)[0:i])
            line31.set_xdata(xage[0:i])             
            line32.set_ydata(spline(gy.age, np.sum(gy.MercDist, axis=1), xage)[0:i])
            line32.set_xdata(xage[0:i])            
            for coll in (ax3.collections):
                ax3.collections.remove(coll)            
            ax3.fill_between(xage[0:i], 0,spline(gy.age, np.sum(gy.MercDist, axis=1), xage)[0:i] ,facecolor='green', alpha=0.3 )

            for coll in (ax4.collections):
                ax4.collections.remove(coll)            
            y = gy.MerchMass
            line41.set_ydata(spline(gy.age, y, xage)[0:i])
            line41.set_xdata(xage[0:i])             
            ax4.fill_between(xage[0:i],0, spline(gy.age, y, xage)[0:i],facecolor='green', alpha=0.3)
            yy = y + gy.BiBark
            line42.set_ydata(spline(gy.age, yy, xage)[0:i])
            line42.set_xdata(xage[0:i])             
            ax4.fill_between(xage[0:i],spline(gy.age, y, xage)[0:i],spline(gy.age, yy, xage)[0:i],facecolor='blue', alpha=0.3)
            y = yy + gy.BiBranch
            line43.set_ydata(spline(gy.age, y, xage)[0:i])
            line43.set_xdata(xage[0:i])             
            ax4.fill_between(xage[0:i],spline(gy.age, yy, xage)[0:i],spline(gy.age, y, xage)[0:i],facecolor='red', alpha=0.3)
            yy = y + gy.BiFoliage
            line44.set_ydata(spline(gy.age, yy, xage)[0:i])
            line44.set_xdata(xage[0:i])             
            ax4.fill_between(xage[0:i],spline(gy.age, y, xage)[0:i],spline(gy.age, yy, xage)[0:i],facecolor='cyan', alpha=0.3)
            line45.set_ydata(spline(gy.age, gy.RootMass*-1.0, xage)[0:i])
            line45.set_xdata(xage[0:i])             
            ax4.fill_between(xage[0:i],spline(gy.age, gy.RootMass*-1.0, xage)[0:i] ,0,facecolor='gray', alpha=0.3)

            for coll in (ax61.collections):
                ax61.collections.remove(coll)            
            a=gy.BaLitD; b = gy.BaLitL 
            line611.set_ydata(spline(gy.age, a, xage)[0:i]) 
            line611.set_xdata(xage[0:i])
            ax61.fill_between(xage[0:i], 0, spline(gy.age, a, xage)[0:i], facecolor='gray', alpha=0.3)
            line612. set_ydata(spline(gy.age, a+b, xage)[0:i])    
            line612.set_xdata(xage[0:i])
            ax61.fill_between(xage[0:i], spline(gy.age, a, xage)[0:i], spline(gy.age, a+b, xage)[0:i], facecolor='green', alpha=0.3)

            for coll in (ax62.collections):
                ax62.collections.remove(coll)            
            a=gy.BrLitD; b = gy.BrLitL 
            line621.set_ydata(spline(gy.age, a,xage)[0:i])
            line621.set_xdata(xage[0:i])
            ax62.fill_between(xage[0:i], 0, spline(gy.age, a, xage)[0:i], facecolor='gray', alpha=0.3)
            line622.set_ydata(spline(gy.age, a+b, xage)[0:i])    
            line622.set_xdata(xage[0:i])
            ax62.fill_between(xage[0:i], spline(gy.age, a, xage)[0:i], spline(gy.age, a+b,xage)[0:i], facecolor='green', alpha=0.3)

            for coll in (ax63.collections):
                ax63.collections.remove(coll)            
            a=gy.FoLitD; b = gy.FoLitL 
            line631.set_ydata(spline(gy.age, a,xage)[0:i])
            line631.set_xdata(xage[0:i])
            ax63.fill_between(xage[0:i], 0, spline(gy.age, a, xage)[0:i], facecolor='gray', alpha=0.3)
            line632.set_ydata(spline(gy.age, a+b, xage)[0:i])    
            line632.set_xdata(xage[0:i])
            ax63.fill_between(xage[0:i], spline(gy.age, a, xage)[0:i], spline(gy.age, a+b,xage)[0:i], facecolor='green', alpha=0.3)

            for coll in (ax64.collections):
                ax64.collections.remove(coll)            
            a=gy.RootLitD; b = gy.FineRLitL 
            line641.set_ydata(spline(gy.age, a,xage)[0:i])
            line641.set_xdata(xage[0:i])
            ax64.fill_between(xage[0:i], 0, spline(gy.age, a, xage)[0:i], facecolor='gray', alpha=0.3)
            line642.set_ydata(spline(gy.age, a+b, xage)[0:i])    
            line642.set_xdata(xage[0:i])
            ax64.fill_between(xage[0:i], spline(gy.age, a, xage)[0:i], spline(gy.age, a+b,xage)[0:i], facecolor='green', alpha=0.3)

            for coll in (ax65.collections):
                ax65.collections.remove(coll)            
            a=gy.CWD 
            line651.set_ydata(spline(gy.age, a,xage)[0:i])
            line651.set_xdata(xage[0:i])
            ax65.fill_between(xage[0:i], 0, spline(gy.age, a, xage)[0:i], facecolor='gray', alpha=0.3)

            for coll in (ax66.collections):
                ax66.collections.remove(coll)            
            a = gy.BaLitD + gy.BrLitD + gy.FoLitD + gy.RootLitD + gy.CWD
            b = gy.BaLitL + gy.BrLitL + gy.FoLitL + gy.FineRLitL
            line661.set_ydata(spline(gy.age, a,xage)[0:i])
            line661.set_xdata(xage[0:i])
            ax66.fill_between(xage[0:i], 0, spline(gy.age, a, xage)[0:i], facecolor='gray', alpha=0.3)
            line662.set_ydata(spline(gy.age, a+b, xage)[0:i])    
            line662.set_xdata(xage[0:i])
            ax66.fill_between(xage[0:i], spline(gy.age, a, xage)[0:i], spline(gy.age, a+b,xage)[0:i], facecolor='green', alpha=0.3)

            for coll in (ax7.collections):
                ax7.collections.remove(coll)            
            a = np.cumsum(a); b= np.cumsum(b)
            line71.set_ydata(spline(gy.age, a,xage)[0:i])
            line71.set_xdata(xage[0:i])
            ax7.fill_between(xage[0:i], 0, spline(gy.age, a, xage)[0:i], facecolor='gray', alpha=0.3)
            line72.set_ydata(spline(gy.age, a+b, xage)[0:i])    
            line72.set_xdata(xage[0:i])
            ax7.fill_between(xage[0:i], spline(gy.age, a, xage)[0:i], spline(gy.age, a+b,xage)[0:i], facecolor='green', alpha=0.3)
    
            return mplfig_to_npimage(fig) # RGB image of the figure
        
        M_ani =mpy.VideoClip(make_frame_mpl, duration=dur)
        M_ani.write_gif("/home/vagrant/.spyder2/Plantation/Hdom.gif", fps=fpers)
            
    #plt.show()



def fDrawDecomposition(gy, de, ro, pe):
    """
    Draws the key variables in decomposition 
    """
    fs = 14.   #fontsize
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)    
    dt = gy.g['dt'][0]
    time = np.arange(0, len(gy.Hdom), step = dt)/12.0
    fig= plt.figure(num = 'Plantation simulator, decomposition', facecolor=fcol, edgecolor='k',figsize=figsi)   #Figsize(w,h), tuple inches 
    ax1 = fig.add_axes([0.33, 0.50, 0.3, 0.46]) #left, bottom, width, height
    ax1.set_title('Soil storage: Organic matter, mass', fontsize=fs +1)    
    ax1.set_xlabel('Time, yrs', fontsize = fs)
    ax1.set_ylabel('ton $ha^{-1}$', fontsize = fs)
    a0 = pe.peatStorage
    a00, = ax1.plot(time, a0 * 1e-3, 'k-')
    ax1.fill_between(time, 0, a0  * 1e-3, facecolor='k', alpha=0.3)
    a =   a0 + de.cwdSto  
    aa, = ax1.plot(time, a, 'k--')
    ax1.fill_between(time, a0 * 1e-3, a * 1e-3, facecolor='0.9', alpha=0.3)
    b = a + ro.SH*10000.0
    bb, = ax1.plot(time, b * 1e-3 , '-')
    ax1.fill_between(time, a * 1e-3, b * 1e-3, facecolor='0.1', alpha=0.3)
    c = b + ro.H*10000.0
    cc, = ax1.plot(time, c * 1e-3, '.')
    ax1.fill_between(time, b * 1e-3, c * 1e-3, facecolor='0.3', alpha=0.3)
    d = c + ro.LH*10000.0
    dd, = ax1.plot(time, d * 1e-3, 'g-')
    ax1.fill_between(time, c * 1e-3, d * 1e-3, facecolor='b', alpha=0.3)
    e = d + ro.F*10000.0   
    ee, = ax1.plot(time, e * 1e-3 , 'y-')
    ax1.fill_between(time, d * 1e-3, e * 1e-3, facecolor='y', alpha=0.3)
    f = e + ro.L*10000.0
    ff, = ax1.plot(time, f * 1e-3, 'c')
    ax1.fill_between(time, e * 1e-3, f * 1e-3, facecolor='g', alpha=0.3)
    ax1.set_ylim(0, max(f)*1.1 * 1e-3)
    ax1.set_xlim(0, max(time))
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)    
    ax1.legend([ff,ee,dd,cc,bb,aa, a00], ['L', 'F', 'LH', 'H', 'SH', 'CWD', 'Peat'], loc=3, fontsize=fs-1)    
    ax11= ax1.twinx()
    subsidence = (f-f[0])/10000.0/pe.peatRhoInit*100
    ax11.plot(time, pe.surfAfterOxidation, 'k--')
    ax11.plot(time, pe.subsidence, 'r--', linewidth = 3)
    ax11.set_ylim( -f[0]/10000.0/pe.peatRhoInit*100, (max(f)*1.1-f[0])/10000.0/pe.peatRhoInit*100)
    ylab = 'Surface elevation, cm' 
    text = 'Initial density ' + str(pe.peatRhoInit) + ' kg $m^{-3}$'    
    plt.figtext(0.5,0.54,text, fontsize=fs)
    text = 'Final density  ' + str(pe.peatRhoFinal) + ' kg $m^{-3}$'    
    plt.figtext(0.5,0.52,text, fontsize=fs)

    ax11.set_ylabel(ylab, fontsize = fs)
    ax11.set_xlim(0, max(time))    
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)    
    
    
    ax2 = fig.add_axes([0.33, 0.08, 0.3, 0.33]) #left, bottom, width, height
    ax2.set_title('Carbon: cumulative release', fontsize=fs+1)        
    ax2.set_xlabel('Time, yrs', fontsize = fs)
    ax2.set_ylabel('ton $ha^{-1}$', fontsize = fs)
    a0 = np.cumsum(pe.peatCrelease) 
    a00, = ax2.plot(time, a0 * 1e-3, 'k-')
    ax2.fill_between(time, 0, a0 * 1e-3, facecolor='k', alpha=0.3)
    a = a0 + np.cumsum(de.cwdCRelease)
    aa, = ax2.plot(time, a * 1e-3, 'k--')
    ax2.fill_between(time, 0, a * 1e-3, facecolor='0.5', alpha=0.3)
    b = a + np.cumsum(ro.miner*10000.0)*ro.massToC
    bb, = ax2.plot(time, b * 1e-3, 'y-')
    ax2.fill_between(time, a * 1e-3, b * 1e-3, facecolor='y', alpha=0.3)
    ax2.set_ylim(0, max(b)*1.1 * 1e-3)
    ax2.set_xlim(0, max(time))    
    ax2.legend([bb,aa, a00], ['Stand litter', 'CWD', 'Peat'], loc=2, fontsize=fs-1)    
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    
    ax4 = fig.add_axes([0.70, 0.08, 0.28, 0.65]) #left, bottom, width, height
    ax4.set_xlabel('Time, yrs', fontsize=fs)
    ax4.set_title('Carbon balance, ton $ha^{-1}$ ', fontsize=fs+1)
    nrota = int(max(gy.rotationArr)+1) 
    standTot = (gy.MerchMass + gy.BiBark + gy.BiBranch  + gy.BiFoliage + gy.RootMass )*ro.massToC * 1e-3 
    standAndWeed = standTot + (gy.weedAbove + gy.weedBelow)*ro.massToC * 1e-3
    standRep = (gy.BiBark + gy.BiBranch  + gy.BiFoliage + gy.RootMass)*ro.massToC * 1e-3
    standRepAndWeed = standRep + (gy.weedAbove + gy.weedBelow)*ro.massToC * 1e-3
    
    """
    for r in range(nrota):    
        s = max(np.where(gy.rotationArr==r,standTot,0))
        standTot = np.where(gy.rotationArr>r,standTot+s,standTot)
        s2= max(np.where(gy.rotationArr==r,standRep,0))
        standRep = np.where(gy.rotationArr>r,standRep+s2,standRep)
    """    
    #soil = b
    soilStoChange =  (f-f[0])*ro.massToC * 1e-3
    soil = (f-f[0])*ro.massToC*-1 * 1e-3    
    balReal = standTot-soil
    balRep = standRep - soil
    ax4.plot(time, standTot , 'g--', label='Stand total')
    ax4.plot(time, standAndWeed, 'm--', label = 'Weed')
    ax4.plot(time, standRep, 'g-', label ='Stand exl stems')
    
    #ax4.plot(time, soil*-1, 'r-', label = 'Soil out')
    ax4.plot(time, balReal, color = '0.2', linestyle='--', label = 'Balance with stem', linewidth = 2)
    ax4.fill_between(time, balReal, 0, where= balReal < 0, facecolor = 'r', alpha=0.3)    
    ax4.fill_between(time, balReal, 0, where= balReal >= 0, facecolor = 'g', alpha=0.3)    
    ax4.plot(time, balRep, color = '0.2', linestyle=':', label = 'Balance exl stem', linewidth = 2)
    ax4.fill_between(time, balRep, 0, where= balRep < 0, facecolor = 'r', alpha=0.3)    
    ax4.fill_between(time, balRep, 0, where= balRep >= 0, facecolor = 'g', alpha=0.3)    
    ax4.plot(time, soilStoChange, color = 'b', linestyle=':', label = 'Soil storage change', linewidth = 2)    
    ax4.set_xlim(0, max(time))    

    ax4.plot(time, np.zeros(len(time)), 'k-')
    ax4.legend(loc='lower left', fontsize=fs-1, ncol=2)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)

    ax7= fig.add_axes([0.05, 0.61, 0.22, 0.35]) #left, bottom, width, height    
    a = gy.FoLitL  + gy.BaLitL + gy.BrLitL  + gy.FineRLitL  
    b =  gy.FoLitD +  gy.BaLitD + gy.BrLitD +  gy.RootLitD
    c = gy.CWD
    a = np.cumsum(a); b= np.cumsum(b); c = np.cumsum(c)
    cc, = ax7.plot(time, c * 1e-3, '0.9', label = 'CWD')
    ax7.fill_between(time, 0, c * 1e-3, facecolor='gray', alpha=0.4)
    bb, = ax7.plot(time, (c+b) * 1e-3, 'k--', label = 'Mort')   
    ax7.fill_between(time, c * 1e-3, (c+b) * 1e-3, facecolor='gray', alpha=0.3)
    aa, = ax7.plot(time, (c+b+a) * 1e-3, 'g', label = 'Living')   
    ax7.fill_between(time, (c+b) * 1e-3, (a+b+c) * 1e-3, facecolor='green', alpha=0.3)
    ax7.set_xlabel('Time, yrs', fontsize = fs)
    ax7.set_ylabel('ton $ha^{-1}$', fontsize = fs)
    
    ax7.set_title('Mass: input to soil as litter', fontsize = fs+1)
    ax7.set_xlim(0, time[-1])    
    ax7.set_ylim(bottom=0)    
    ax7.legend([aa,bb,cc], ['Living', 'Mort', 'CWD'], loc=2, fontsize=fs-1)    
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)

    ax8= fig.add_axes([0.05, 0.08, 0.22, 0.33]) #left, bottom, width, height    
    a = (gy.FoLitL  + gy.BaLitL + gy.BrLitL  + gy.FineRLitL)*ro.massToC  
    b =  (gy.FoLitD +  gy.BaLitD + gy.BrLitD +  gy.RootLitD)*ro.massToC
    c = gy.CWD*de.massToC
    a = np.cumsum(a); b= np.cumsum(b); c = np.cumsum(c)
    cc, = ax8.plot(time, c * 1e-3, '0.9', label = 'CWD')
    ax8.fill_between(time, 0, c * 1e-3, facecolor='gray', alpha=0.4)
    bb, = ax8.plot(time, (c+b) * 1e-3, 'k--', label = 'Mort')   
    ax8.fill_between(time, c * 1e-3, (c+b) * 1e-3, facecolor='gray', alpha=0.3)
    aa, = ax8.plot(time, (c+b+a) * 1e-3, 'g', label = 'Living')   
    ax8.fill_between(time, (c+b) * 1e-3, (a+b+c) * 1e-3, facecolor='green', alpha=0.3)
    ax8.set_xlabel('Time, yrs', fontsize = fs)
    ax8.set_ylabel('ton $ha^{-1}$', fontsize = fs)
    
    ax8.set_title('Carbon: input to soil as litter', fontsize = fs+1)
    ax8.set_xlim(0, time[-1])    
    ax8.set_ylim(bottom=0)    
    ax8.legend([aa,bb,cc], ['Living', 'Mort', 'CWD'], loc='upper left', fontsize=fs-1)    
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    #----meatadata text: Table  ---------------  
    optMetadata = False
    if optMetadata:
        le = 0.66; w = 0.045    
        row0 = 0.93; rh =0.03
        plt.figtext(le+0*w, row0+rh, 'Carbon storages and fluxes, ton $ha^{-1}$', fontsize = fs, fontweight='bold')
        plt.figtext(le+0*w, row0, 'Rotat', fontsize = fs, fontweight='bold')
        plt.figtext(le+1.0*w, row0, 'Export', fontsize = fs, fontweight='bold')        
        plt.figtext(le+2.0*w, row0, 'Resid', fontsize = fs, fontweight='bold')
        plt.figtext(le+3.0*w, row0, 'Weeds', fontsize = fs, fontweight='bold')        
        plt.figtext(le+4.0*w, row0, 'Input', fontsize = fs, fontweight='bold')
        plt.figtext(le+5.0*w, row0, 'Emiss', fontsize = fs,fontweight='bold')        
        plt.figtext(le+6.0*w, row0, 'dS', fontsize = fs,fontweight='bold')   
        plt.figtext(le+7.0*w, row0, 'Sb', fontsize = fs,fontweight='bold')   
        
        Nrotat = int(gy.g['Nrotat'][0]); rotLength = len(gy.age)/Nrotat; lit = 0 
            
        for n in range(Nrotat):
            rstart = int((n)*rotLength)
            rend = int((n+1)*rotLength-1)       
            rg = gy.rowGe[n]  
            mm = gy.MerchMass[rend]
            plt.figtext(le+0*w, row0-(n+1)*rh, str(n+1), fontsize = fs)
            plt.figtext(le+1.0*w, row0-(n+1)*rh, str(mm / 1000.0*de.massToC)[0:4], fontsize = fs) 
            r = (gy.BiBark[-1] + gy.BiBranch[rend]  + gy.BiFoliage[rend] + gy.RootMass[rend] )*ro.massToC / 1000.0         
            plt.figtext(le+2.0*w, row0-(n+1)*rh, str(r)[0:5], fontsize = fs)        
            we = (gy.weedAbove[rend] + gy.weedBelow[rend])*ro.massToC/1000.0 
            plt.figtext(le+3.0*w, row0-(n+1)*rh, str(we)[0:5], fontsize = fs)        
            inp =  a+b+c       
            plt.figtext(le+4.0*w, row0-(n+1)*rh, str(inp[rend] / 1000.0)[0:5], fontsize = fs)
            emiss = np.cumsum(pe.peatCrelease) + np.cumsum(de.cwdCRelease) +np.cumsum(ro.miner*10000.0)*ro.massToC     
            plt.figtext(le+5.0*w, row0-(n+1)*rh, str((emiss[rend]-emiss[rstart])/1000.0)[0:5], fontsize = fs)           
            soil = (de.cwdSto + (ro.SH + ro.H + ro.LH + ro.F + ro.L)*10000.0 + pe.peatStorage)*ro.massToC 
            ds = (soil[rend] -soil[rstart])
            plt.figtext(le+6.0*w, row0-(n+1)*rh, str(ds/1000.0)[0:5], fontsize = fs)    
            su = (pe.subsidence[rend] -pe.subsidence[rstart])/rotLength*12.0        
            plt.figtext(le+7.0*w, row0-(n+1)*rh, str(su)[0:5], fontsize = fs)    
        
    plt.show()

def drawCderivative(gy, de, ro, pe):

    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)    
    dt = gy.g['dt'][0]
    time = np.arange(0, len(gy.Hdom), step = dt)/12.0
    fig= plt.figure(num = 'C rates, decomposition', facecolor=fcol, edgecolor='k',figsize=figsi)   #Figsize(w,h), tuple inches 


    a0 = pe.peatStorage
    a =   a0 + de.cwdSto  
    b = a + ro.SH*10000.0
    c = b + ro.H*10000.0
    d = c + ro.LH*10000.0
    e = d + ro.F*10000.0   
    f = e + ro.L*10000.0

    ax4 = fig.add_axes([0.2, 0.2, 0.6, 0.55]) #left, bottom, width, height
    ax4.set_xlabel('Time, yrs', fontsize=14)
    ax4.set_title('Carbon balance, kg $ha^{-1} month^{-1}$ ')
    nrota = int(max(gy.rotationArr)+1)
    standTot = (gy.MerchMass + gy.BiBark + gy.BiBranch  + gy.BiFoliage + gy.RootMass )*ro.massToC 
    standAndWeed = standTot + (gy.weedAbove + gy.weedBelow)*ro.massToC
    standRep = (gy.BiBark + gy.BiBranch  + gy.BiFoliage + gy.RootMass)*ro.massToC
    standRepAndWeed = standRep + (gy.weedAbove + gy.weedBelow)*ro.massToC
    """
    for r in range(nrota):    
        s = max(np.where(gy.rotationArr==r,standTot,0))
        standTot = np.where(gy.rotationArr>r,standTot+s,standTot)
        s2= max(np.where(gy.rotationArr==r,standRep,0))
        standRep = np.where(gy.rotationArr>r,standRep+s2,standRep)
    """    
    #soil = b

    soilStoChange =  (f-f[0])*ro.massToC
    soil = (f-f[0])*ro.massToC*-1    
    balReal = standTot-soil
    balRep = standRep - soil
    ax4.plot(time, np.gradient(standTot) , 'g--', label='Stand total')
    ax4.plot(time, np.gradient(standAndWeed), 'm--', label = 'Weed')
    ax4.plot(time, np.gradient(standRep), 'g-', label ='Stand exl stems')
    
    #ax4.plot(time, soil*-1, 'r-', label = 'Soil out')
    ax4.plot(time, np.gradient(balReal), color = '0.2', linestyle='--', label = 'Balance with stem', linewidth = 2)
    ax4.fill_between(time, np.gradient(balReal), 0, where= np.gradient(balReal) < 0, facecolor = 'r', alpha=0.3)    
    ax4.fill_between(time, np.gradient(balReal), 0, where= np.gradient(balReal) >= 0, facecolor = 'g', alpha=0.3)    
    ax4.plot(time, np.gradient(balRep), color = '0.2', linestyle=':', label = 'Balance exl stem', linewidth = 2)
    ax4.fill_between(time, np.gradient(balRep), 0, where= np.gradient(balRep) < 0, facecolor = 'r', alpha=0.3)    
    ax4.fill_between(time, np.gradient(balRep), 0, where= np.gradient(balRep) >= 0, facecolor = 'g', alpha=0.3)    
    ax4.plot(time, np.gradient(soilStoChange), color = 'b', linestyle=':', label = 'Soil storage change', linewidth = 2)    
    ax4.set_xlim(0, max(time))    

    ax4.plot(time, np.zeros(len(time)), 'k-')
    ax4.legend(loc='upper left')
    plt.show()

def fDrawGy(gy):
   
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)    
        
    plt.figure(num = 'Growth and yield variables', facecolor=fcol, edgecolor='k',figsize=figsi)   #Figsize(w,h), tuple inches    
    titles =['Dominant height', 'Basal area', 'Survival', 'Mean diameter', 'Pulp yield', 'Wood consumption']
    ylab = ['m', '$m^{2}$', 'number of stems $ha^{-1}$', 'cm', 'tn $ha^{-1}$','$m{3}$ $ha^{-1}$']
    xlab = ['age, yrs']*6
    ys = [gy.Hdom, gy.BA, gy.Survival, gy.WMeanDiam, gy.PulpY, gy.WoodConsumption]
    ind = range(6)   
    for f in ind:
        plt.subplot(2, 3, f+1)
        plt.plot(gy.time, ys[f], 'r.-')
        plt.title(titles[f])
        plt.ylabel(ylab[f])
        plt.xlabel(xlab[f])
    plt.show()
    
    plt.figure(num = 'Distributions', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(25.0,12.0))   #Figsize(w,h), tuple inches    
    titles =['Diameter, cm', 'Volume, $m^{3}$ $ha^{-1}$', 'Height, m', 'Merchantable vol, $m^{3}$ $ha^{-1}$', 'Basal area, $m^{2}$ $ha^{-1}$']
    ylab = ''
    xlab = ['diameter, cm']*5
    ind = range(4)   
    xs = range(40)
    ul = int(gy.g['DiamUpperLimit'][0]); nsteps = len(gy.age)  
    DbhDist = np.reshape(gy.DbhDist, (nsteps, ul))
    Vdist= np.reshape(gy.Vdist, (nsteps, ul))
    Hdist= np.reshape(gy.Hdist, (nsteps, ul))
    MercDist= np.reshape(gy.MercDist, (nsteps, ul))
    ys = [DbhDist, Vdist, Hdist, MercDist]
    ages=[12,24,36,48,60]
    colors =['r.-', 'b.-','y.-','g.-','m.-', 'k.-']
    for f in ind:
        count = 0
        for a in ages:
            plt.subplot(2, 2, f+1)
            plt.plot(xs, ys[f][a][0:], colors[count])
            plt.title(titles[f])
            plt.ylabel(ylab)
            plt.xlabel(xlab[f])
            count +=1
    plt.show()
    
    plt.figure(num = 'Weeds', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(25.0,12.0))   #Figsize(w,h), tuple inches    
    titles =['Above ground weed', 'Below ground weed', 'Weed LAI', 'Above litter from living', 'Below ground living', 'Above litter from dead', 'Below ground dead']
    ylab = ['kg $ha^{-1}$', 'kg $ha^{-1}$', 'kg $ha-^{1} $timestep^{-1}', 'kg $ha-^{1} $timestep^{-1}', 'kg $ha-^{1} $timestep^{-1}', 'kg $ha-^{1} $timestep^{-1}', 'kg $ha-^{1} $timestep^{-1}', 'kg $ha-^{1} $timestep^{-1}']
    xlab = ['age, yrs']*7
    ys = [gy.weedAbove, gy.weedBelow, gy.weedLAI, gy.weedALitL, gy.weedBLitL, gy.weedALitD, gy.weedBLitD]
    ind = range(7)   
    for f in ind:
        plt.subplot(3, 3, f+1)
        plt.plot(gy.time, ys[f], 'r.-')
        plt.title(titles[f])
        plt.ylabel(ylab[f])
        plt.xlabel(xlab[f])
    plt.tight_layout()
    plt.show()
    
    plt.figure(num = 'Biomass', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(25.0,12.0))   #Figsize(w,h), tuple inches    
    titles =['Bark', 'Branch', 'Foliage', 'Stem', 'Root', 'Merchantable wood', 'LAI']
    ylab = ['kg $ha^{-1}$', 'kg $ha^{-1}$', 'kg $ha-^{1}$', 'kg $ha-^{1}$', 'kg $ha-^{1}$', 'kg $ha-^{1}$ ',  '$m-^{2} $ $m-^{-2} $']
    xlab = ['age, m']*7
    ys = [gy.BiBark, gy.BiBranch, gy.BiFoliage, gy.BiStem, gy.RootMass, gy.MerchMass, gy.LAI]
    ind = range(7)   
    for f in ind:
        plt.subplot(3, 3, f+1)
        plt.plot(gy.time, ys[f], 'r.-')
        plt.title(titles[f])
        plt.ylabel(ylab[f])
        plt.xlabel(xlab[f])
    plt.tight_layout()
    plt.show()
    
    plt.figure(num = 'Stand litters', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(25.0,12.0))   #Figsize(w,h), tuple inches    
    titles =['Roots, living', 'Roots, dead', 'Branch, living', 'Branch, dead', 'Bark, living', 'Bark, dead', 
             'Foliage, living', 'Foliage, dead', 'Coarse woody debris']
    ylab = 'kg $ha-^{1}$ $timestep^{-1}$'
    xlab = 'age, yrs'
    ys = [gy.FineRLitL, gy.RootLitD, gy.BrLitL, gy.BrLitD, gy.BaLitL, gy.BaLitD, gy.FoLitL, gy.FoLitD, gy.CWD]
    ind = range(9)   
    for f in ind:
        plt.subplot(3, 3, f+1)
        plt.plot(gy.time, ys[f], 'r.-')
        plt.title(titles[f])
        plt.ylabel(ylab)
        plt.xlabel(xlab)
    plt.tight_layout()
    plt.show()
 
def fDrawDec(gy, de, ro, pe):
   
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)    
        
    plt.figure(num = 'Decomposition, coarse woody debris', facecolor=fcol, edgecolor='k',figsize=figsi)   #Figsize(w,h), tuple inches    
    titles =['CWD input', 'CWD storage', 'CWD mass loss', 'C input as CWD', 'C storage in CWD', 'C emission from CWD']
    ylab = ['kg $ha^{-1}$ $timestep^{-1}$', 'kg $ha^{-1}$', 'kg $ha^{-1}$ $timestep^{-1}$', 'kg $ha^{-1}$ $timestep^{-1}$', 'kg $ha^{-1}$', 'kg $ha^{-1}$ $timestep^{-1}$']
    xlab = ['age, yrs']*6
    ys = [gy.CWD, de.cwdSto, de.cwdMloss, de.cwdCInp, de.cwdCSto, de.cwdCRelease]
    ind = range(6)   
    for f in ind:
        plt.subplot(2, 3, f+1)
        plt.plot(gy.time, ys[f], 'r.-')
        plt.title(titles[f])
        plt.ylabel(ylab[f])
        plt.xlabel(xlab[f])
    plt.tight_layout()
    plt.show()
        
    plt.figure(num = 'Decomposition, plant litter, mass', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(25.0,12.0))   #Figsize(w,h), tuple inches    
    titles =['Total input', 'L storage', 'F storage', 'LH storage', 'H storage', 'SH storage', 
             'Total storage', 'Mass loss', 'Cumulatrive input', 'Cumulative mass loss' ]
    ylab = ['kg $m^{-2}$ $timestep^{-1}$']*7 + ['kg $m^{-2}$ $timestep^{-1}$'] + ['kg $m^{-2}$']*2
    xlab = 'age, yrs'
    ys = [ro.totalLitter, ro.L, ro.F, ro.LH, ro.H, ro.H, ro.totSto[1:], ro.miner, np.cumsum(ro.totalLitter), np.cumsum(ro.miner)]
    ind = range(10)   
    for f in ind:
        plt.subplot(3, 4, f+1)
        plt.plot(gy.time, ys[f], 'r.-')
        plt.title(titles[f])
        plt.ylabel(ylab[f])
        plt.xlabel(xlab)
    plt.tight_layout()
    plt.show()
    
    
    plt.figure(num = 'Decomposition, plant litter, Nitrogen', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(25.0,12.0))   #Figsize(w,h), tuple inches    
    titles =['Total N input', 'L storage', 'F storage', 'LH storage', 'H storage', 'SH storage', 
             'Total storage', 'N release', 'Cumulative N input', 'Cumulative N release']
    ylab = ['kg N $m^{-2}$ $timestep^{-1}$']*7 + ['kg N $m^{-2}$ $timestep^{-1}$'] + ['kg N $m^{-2}$']*2
    xlab = 'age, yrs'
    ys = [ro.totalNLitter, ro.NL, ro.NF, ro.NLH, ro.NH, ro.NH, ro.totNSto[1:], ro.Nminer, np.cumsum(ro.totalNLitter), np.cumsum(ro.Nminer)]
    ind = range(10)   
    for f in ind:
        plt.subplot(3, 4, f+1)
        plt.plot(gy.time, ys[f], 'r.-')
        plt.title(titles[f])
        plt.ylabel(ylab[f])
        plt.xlabel(xlab)
    plt.tight_layout()
    plt.show()
    

    plt.figure(num = 'Decomposition, peat', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(25.0,12.0))   #Figsize(w,h), tuple inches    
    titles =['Peat storage', 'Peat C release', 'Cumulative C release from peat']
    ylab = ['kg $ha^{-1}$ ', 'kg $ha^{-1}$ $timestep^{-1}$','kg $ha^{-1}$ '] 
    xlab = 'age, yrs'
    ys = [pe.peatStorage, pe.peatCrelease, np.cumsum(pe.peatCrelease)]
    ind = range(3)   
    for f in ind:
        plt.subplot(2, 2, f+1)
        plt.plot(gy.time, ys[f], 'r.-')
        plt.title(titles[f])
        plt.ylabel(ylab[f])
        plt.xlabel(xlab)
    plt.tight_layout()
    plt.show()    
    
def fShowActiveFigs(rounds=1, tauko=1.0):
    import matplotlib.pyplot as plt
    import time
    fnames = plt.get_figlabels()
    for r in range(rounds):        
        for f in fnames:
            print (f)
            show_plot(f)        
            time.sleep(tauko)

def show_plot(figure_id=None):    
    if figure_id is not None:
        fig = plt.figure(num=figure_id)
    else:
        fig = plt.gcf()

    plt.show()
    fig.canvas.manager.window.activateWindow()
    fig.canvas.manager.window.raise_()

def fDrawNutrient(n, gy, de):
    fs = 14
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)    
    dt = gy.g['dt'][0]; dty = dt/12.0
    time = np.arange(0, len(n.stemNet))*dty
    fig= plt.figure(num = n.name, facecolor=fcol, edgecolor='k',figsize=figsi)   #Figsize(w,h), tuple inches 
    
    ax1 = fig.add_axes([0.07, 0.50, 0.25, 0.46]) #left, bottom, width, height
    ax1.set_title('Net demand', fontsize=fs)    
    #ax1.set_xlabel('Age, yrs', fontsize = 14)
    ax1.set_ylabel('Net demand kg $ha^{-1}$', fontsize=fs)
    y = n.stemNet
    ax1.plot(time, y , 'g-', label='Stem')
    ax1.fill_between(time, 0, y,facecolor='green', alpha=0.3)
    yy = y + n.barkNet
    ax1.plot(time, yy , 'b-', label='Bark')
    ax1.fill_between(time,y,yy,facecolor='blue', alpha=0.3)
    y = yy + n.branchNet
    ax1.plot(time, y , 'r-', label='Branch')
    ax1.fill_between(time,yy,y,facecolor='red', alpha=0.3)
    yy = y + n.foliageNet
    ax1.plot(time, yy , 'c-', label='Foliage')
    ax1.fill_between(time,y,yy,facecolor='cyan', alpha=0.3)
    y = yy + n.weedAboveNet
    ax1.plot(time, y, 'y-', label = 'Weed')
    ax1.fill_between(time, yy, y, facecolor = 'yellow', alpha=0.3)    
    y = n.rootNet*-1.0
    ax1.plot(time, y , 'k-', label='Root')
    ax1.fill_between(time, y ,0,facecolor='gray', alpha=0.3)
    yy = y -  n.fineRootNet
    ax1.plot(time, yy, 'k--', label = 'Fineroot' )
    ax1.fill_between(time, yy ,y,facecolor='gray', alpha=0.5)    
    y = yy - n.weedBelowNet
    ax1.plot(time, y, 'm-', label='Weed roots')
    ax1.fill_between(time, y, yy, facecolor='magenta', alpha=0.3)    
    ax1.legend(loc='upper right', fontsize=fs-1, ncol=2)
    ax1.set_xlim(0,time[-1])
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)


    ax2 = fig.add_axes([0.07, 0.07, 0.25, 0.35]) #left, bottom, width, height
    ax2.set_title('Nutrients in vegetation', fontsize=fs)    
    ax2.set_xlabel('Age, yrs', fontsize = fs)
    ax2.set_ylabel('Storage kg $ha^{-1}$', fontsize=fs)
    y = n.stemSto
    ax2.plot(time, y , 'g-', label='Stem')
    ax2.fill_between(time, 0, y,facecolor='green', alpha=0.3)
    yy = y + n.barkSto
    ax2.plot(time, yy , 'b-', label='Bark')
    ax2.fill_between(time,y,yy,facecolor='blue', alpha=0.3)
    y = yy + n.branchSto
    ax2.plot(time, y , 'r-', label='Branch')
    ax2.fill_between(time,yy,y,facecolor='red', alpha=0.3)
    yy = y + n.foliageSto
    ax2.plot(time, yy , 'c-', label='Foliage')
    ax2.fill_between(time,y,yy,facecolor='cyan', alpha=0.3)
    y = yy + n.weedAboveSto
    
    ax2.plot(time, y, 'y-', label = 'Weed')
    ax2.fill_between(time, yy, y, facecolor = 'yellow', alpha=0.3)
    y = n.rootSto*-1.0
    ax2.plot(time, y , 'k-', label='Root')
    ax2.fill_between(time, y ,0,facecolor='gray', alpha=0.3)    
    yy = y - n.fineRootSto
    ax2.plot(time, yy, 'k--', label='Fineroot')
    ax2.fill_between(time, yy ,y,facecolor='gray', alpha=0.5)    
    y = yy - n.weedBelowSto
    ax2.plot(time, y, 'm-', label='Weed roots')
    ax2.fill_between(time, y, yy, facecolor='magenta', alpha=0.3)
    ax2.legend(loc='upper right', fontsize=fs, ncol=2)
    ax2.set_xlim(0,time[-1])
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)


    ax3 = fig.add_axes([0.40, 0.50, 0.25, 0.46]) #left, bottom, width, height
    ax3.set_title('Gross demand', fontsize=fs)    
    #ax3.set_xlabel('Age, yrs', fontsize = 14)
    ax3.set_ylabel('Gross demand kg $ha^{-1}$', fontsize=fs)
    y = n.stemGro
    ax3.plot(time, y , 'g-', label='Stem')
    ax3.fill_between(time, 0, y,facecolor='green', alpha=0.3)
    yy = y + n.barkGro
    ax3.plot(time, yy , 'b-', label='Bark')
    ax3.fill_between(time,y,yy,facecolor='blue', alpha=0.3)
    y = yy + n.branchGro
    ax3.plot(time, y , 'r-', label='Branch')
    ax3.fill_between(time,yy,y,facecolor='red', alpha=0.3)
    yy = y + n.foliageGro
    ax3.plot(time, yy , 'c-', label='Foliage')
    ax3.fill_between(time,y,yy,facecolor='cyan', alpha=0.3)
    y = yy + n.weedAboveGro
    ax3.plot(time, y, 'y-', label = 'Weed')
    ax3.fill_between(time, yy, y, facecolor = 'yellow', alpha=0.3)    
    y = n.rootGro*-1.0
    ax3.plot(time, y , 'k-', label='Root')
    ax3.fill_between(time, y ,0,facecolor='gray', alpha=0.3)
    yy = y -  n.fineRootGro
    ax3.plot(time, yy, 'k--', label = 'Fineroot' )
    ax3.fill_between(time, yy ,y,facecolor='gray', alpha=0.5)    
    y = yy - n.weedBelowGro
    ax3.plot(time, y, 'm-', label='Weed roots')
    ax3.fill_between(time, y, yy, facecolor='magenta', alpha=0.3)    
    ax3.legend(loc='upper right', fontsize=fs-1, ncol=2)
    ax3.set_xlim(0,time[-1]) 
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
 

    ax4 = fig.add_axes([0.40, 0.07, 0.25, 0.35]) #left, bottom, width, height
    ax4.set_title('Cumulative gross demand', fontsize=fs)    
    ax4.set_xlabel('Age, yrs', fontsize = fs)
    ax4.set_ylabel('Gross demand kg $ha^{-1}$', fontsize=fs)
    y = np.cumsum(n.stemGro)
    ax4.plot(time, y , 'g-', label='Stem')
    ax4.fill_between(time, 0, y,facecolor='green', alpha=0.3)
    yy = y + np.cumsum(n.barkGro)
    ax4.plot(time, yy , 'b-', label='Bark')
    ax4.fill_between(time,y,yy,facecolor='blue', alpha=0.3)
    y = yy + np.cumsum(n.branchGro)
    ax4.plot(time, y , 'r-', label='Branch')
    ax4.fill_between(time,yy,y,facecolor='red', alpha=0.3)
    yy = y + np.cumsum(n.foliageGro)
    ax4.plot(time, yy , 'c-', label='Foliage')
    ax4.fill_between(time,y,yy,facecolor='cyan', alpha=0.3)
    y = yy + np.cumsum(n.weedAboveGro)
    ax4.plot(time, y, 'y-', label = 'Weed')
    ax4.fill_between(time, yy, y, facecolor = 'yellow', alpha=0.3)    
    y = np.cumsum(n.rootGro)*-1.0
    ax4.plot(time, y , 'k-', label='Root')
    ax4.fill_between(time, y ,0,facecolor='gray', alpha=0.3)
    yy = y -  np.cumsum(n.fineRootGro)
    ax4.plot(time, yy, 'k--', label = 'Fineroot' )
    ax4.fill_between(time, yy ,y,facecolor='gray', alpha=0.5)    
    y = yy - np.cumsum(n.weedBelowGro)
    ax4.plot(time, y, 'm-', label='Weed roots')
    ax4.fill_between(time, y, yy, facecolor='magenta', alpha=0.3)    
    ax4.legend(loc='upper left', fontsize=fs, ncol=2)
    ax4.set_xlim(0,time[-1]) 
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    
    totDemand = np.cumsum(n.stemGro) + np.cumsum(n.barkGro) + \
        np.cumsum(n.branchGro) +np.cumsum(n.foliageGro) + \
        np.cumsum(n.weedAboveGro) + np.cumsum(n.rootGro) + \
        np.cumsum(n.fineRootGro) +  np.cumsum(n.weedBelowGro)
        
        
    ax5 = fig.add_axes([0.72, 0.50, 0.25, 0.46]) #left, bottom, width, height
    ax5.set_title('Nutrient supply', fontsize=fs)    
    #ax5.set_xlabel('Age, yrs', fontsize = 14)
    ax5.set_ylabel('Supply kg $ha^{-1}$ $timestep^{-1}$', fontsize=fs)
    y = n.weathering
    ax5.plot(time, y , 'g-', label='Weathering')
    ax5.fill_between(time, 0, y,facecolor='green', alpha=0.3)
    yy = y + n.deposition
    ax5.plot(time, yy , 'b-', label='Deposition')
    ax5.fill_between(time,y,yy,facecolor='blue', alpha=0.3)
    y = yy + n.fixing
    ax5.plot(time, y , 'r-', label='Micr fixing')
    ax5.fill_between(time,yy,y,facecolor='red', alpha=0.3)
    yy = y + n.retranslocation
    ax5.plot(time, yy , 'c-', label='Retranslocation')
    ax5.fill_between(time,y,yy,facecolor='cyan', alpha=0.3)        

    yyy=yy+n.retranslocationWeeds
    ax5.plot(time, yyy ,  linestyle='solid', color='#F0FFFF', label='Weed retranslocation')
    ax5.fill_between(time,yy,yyy,facecolor= '#F0FFFF', alpha=0.3)        

    y = yyy + n.releasePeat
    ax5.plot(time, y , 'k-', label='Peat')
    ax5.fill_between(time, yyy ,y,facecolor='gray', alpha=0.3) 
    yy = y +  n.releaseCwd
    ax5.plot(time, yy, 'k--', label = 'CWD' )
    ax5.fill_between(time, y ,yy,facecolor='gray', alpha=0.6)
    y = yy + n.releaseLitter
    ax5.plot(time, y, 'm-', label='Litter')
    ax5.fill_between(time, yy, y, facecolor='magenta', alpha=0.3)     
    yy = y + n.ferRel 
    ax5.plot(time, yy, 'y-', label = 'Fertilizer')
    ax5.fill_between(time, y, yy, facecolor = 'yellow', alpha=0.3)    
    ax5.legend(loc='upper left')
    ax5.set_xlim(0,time[-1]) 
    
    ax5.plot(time, n.totalDemand, 'r--', linewidth = 3, label ='Total demand')    
    ax5.legend(loc='upper right', fontsize=fs-2, ncol=2)
    ax5.set_xlim(0,time[-1]) 
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    

    ax6 = fig.add_axes([0.72, 0.07, 0.25, 0.35]) #left, bottom, width, height
    ax6.set_title('Cumulative nutrient supply', fontsize=fs)    
    ax6.set_xlabel('Age, yrs', fontsize = fs)
    ax6.set_ylabel('Supply kg $ha^{-1}$ ', fontsize=fs)
    y = np.cumsum(n.weathering)
    ax6.plot(time, y , 'g-', label='Weathering')
    ax6.fill_between(time, 0, y,facecolor='green', alpha=0.3)
    yy = y + np.cumsum(n.deposition)
    ax6.plot(time, yy , 'b-', label='Deposition')
    ax6.fill_between(time,y,yy,facecolor='blue', alpha=0.3)
    y = yy + np.cumsum(n.fixing)
    ax6.plot(time, y , 'r-', label='Micr fixing')
    ax6.fill_between(time,yy,y,facecolor='red', alpha=0.3)
    yy = y + np.cumsum(n.retranslocation)
    ax6.plot(time, yy , 'c-', label='Retranslocation')
    ax6.fill_between(time,y,yy,facecolor='cyan', alpha=0.3)        

    yyy=yy+np.cumsum(n.retranslocationWeeds)
    ax6.plot(time, yyy ,  linestyle='solid', color='#F0FFFF', label='Weed retranslocation')
    ax6.fill_between(time,yy,yyy,facecolor= '#F0FFFF', alpha=0.3)        


    y = yyy + np.cumsum(n.releasePeat)
    ax6.plot(time, y , 'k-', label='Peat')
    ax6.fill_between(time, yyy ,y,facecolor='gray', alpha=0.3) 
    yy = y +  np.cumsum(n.releaseCwd)
    ax6.plot(time, yy, 'k--', label = 'CWD' )
    ax6.fill_between(time, y ,yy,facecolor='gray', alpha=0.6)
    y = yy + np.cumsum(n.releaseLitter)
    ax6.plot(time, y, 'm-', label='Litter')
    ax6.fill_between(time, yy, y, facecolor='magenta', alpha=0.3)     
    yy = y + np.cumsum(n.ferRel) 
    ax6.plot(time, yy, 'y-', label = 'Fertilizer')
    ax6.fill_between(time, y, yy, facecolor = 'yellow', alpha=0.3)    
   

    ax6.plot(time, np.cumsum(n.totalDemand), 'r--', linewidth = 3, label ='Total demand')    
    ax6.legend(loc='upper left', fontsize=fs-2, ncol=2)
    ax6.set_xlim(0,time[-1]) 
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)

    plt.show()

def fDrawValid(pe, de, ro, gy):
    import matplotlib.pyplot as plt
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)    
    fig= plt.figure(num = 'Plantation simulator test', facecolor=fcol, edgecolor='k',figsize=figsi)   #Figsize(w,h), tuple inches 
    
    ax1 = fig.add_axes([0.03, 0.52, 0.3, 0.4]) #left, bottom, width, height
    ax1.set_title('Dwt')    
    ax1.set_xlabel('time, months', fontsize = 14)
    ax1.set_ylabel('dwt, cm')
    ax1.plot(pe.time, pe.dwt*-1 , 'g-', label='dwt')
     
    text = 'Mean dwt ' + str(np.mean(pe.dwt)*-1)[:5] + ' cm'    
    plt.figtext(0.08,0.9,text, fontsize=14)
    
    
    ax2 = fig.add_axes([0.03, 0.05, 0.3, 0.4]) #left, bottom, width, height
    ax2.set_title('Soil surface elevation')    
    ax2.set_xlabel('time, months', fontsize = 14)
    ax2.set_ylabel('Elevation cm')

    f = pe.peatStorage + de.cwdSto + ro.SH*10000.0 + \
        ro.H*10000.0 + ro.LH*10000.0 + ro.F*10000.0   + \
        ro.L*10000.0  

    subsidence = (f-f[0])/10000.0/pe.peatRhoInit*100
    ax2.plot(pe.time, pe.surfAfterOxidation, 'k--')
    ax2.plot(pe.time, pe.subsidence, 'r--', linewidth = 3)
    ax2.set_ylim( -f[0]/10000.0/pe.peatRhoInit*100, (max(f)*1.1-f[0])/10000.0/pe.peatRhoInit*100)
    ylab = 'Surface elevation, cm' 
    text = 'Initial density ' + str(pe.peatRhoInit) + ' kg $m^{-3}$'    
    plt.figtext(0.08,0.2,text, fontsize=14)
    text = 'Final density  ' + str(pe.peatRhoFinal) + ' kg $m^{-3}$'    
    plt.figtext(0.08,0.17,text, fontsize=14)
    text = 'Total subsidence  ' + str(pe.subsidence[-1])[:6] + ' cm'    
    plt.figtext(0.08,0.14,text, fontsize=14)

    ax2.set_ylabel(ylab, fontsize = 14)
    ax2.set_xlim(0, max(pe.time))    

    import seaborn as sns
    ax3 = fig.add_axes([0.38, 0.05, 0.3, 0.4]) #left, bottom, width, height
    ax3.set_title('Subsidence rate')    
    ax3.set_xlabel('rate, cm $yr^{-1}$', fontsize = 14)
    ax3.set_ylabel('frequency')
    rates =  np.diff(pe.subsidence)*12.0  
    ax3 = sns.kdeplot(rates, shade=True, color="r")
    ax3.set_xlim(-20,10)
    text = 'mean  ' + str(np.average(rates))[:5] + ' cm $yr^{-1}$' 
    plt.figtext(0.41,0.40,text, fontsize=14)

    text = 'sd ' +str(np.std(rates))[:5] + ' cm $yr^{-1}$'
    plt.figtext(0.41,0.37,text, fontsize=14)
   
    ax4 = fig.add_axes([0.38, 0.52, 0.3, 0.4]) #left, bottom, width, height
    ax4.set_title('Stand')    
    ax4.set_xlabel('time, yrs', fontsize = 14)
    ax4.set_ylabel('dwt, cm')
    ax4.plot(gy.time, gy.Hdom, 'r-', label='$H_{dom}$, m')
    ax4.plot(gy.time, gy.BA, 'b-', label='BA, $m^{2}$')
    ax4.plot(gy.time[1:], gy.WMeanDiam[1:], 'c-', label='$\overline{dbh}$, cm')
    ax4.legend(loc='upper left')
    ax11=ax4.twinx()
    ax11.plot(gy.time, gy.Survival, 'g-', label='Stocking')
    ax11.legend(loc='lower right')
    ax11.set_ylim(0,max(gy.Survival)*1.2)
    
    ul = int(gy.g['DiamUpperLimit'][0]); nsteps = len(gy.Hdom)
    MercVol = np.reshape(gy.MercDist, (nsteps, ul))
    ax4 = fig.add_axes([0.73, 0.52, 0.25, 0.4]) #left, bottom, width, height
    ax4.set_title('Merchantable volume')    
    ax4.set_xlabel('time, yrs', fontsize = 14)
    ax4.set_ylabel('Merc vol, $m^{3}$ $ha^{-1}$')
    ax4.plot(gy.time, np.sum(MercVol, axis=1), 'r-', label='vol')
    ax4.legend(loc='lower right')

    reft = -26
    V = np.sum(MercVol, axis=1)[reft]
    text = 'Merchantable volume  ' +str(V)[:5] + ' $m^{3}$ $ha^{-1}$ at age ' + str(len(gy.Hdom)+reft) + ' months'
    plt.figtext(0.73,0.225,text, fontsize=14)        
    hdom  = gy.Hdom[reft]
    text = '$H_{dom}$  ' +str(hdom)[:4] + ' m at age ' + str(len(gy.Hdom)+reft) + ' months'
    plt.figtext(0.73,0.2,text, fontsize=14)    
    stocking = gy.Survival[reft]
    text = 'Survival  ' +str(int(stocking))[:4] + ' stems at age ' + str(len(gy.Hdom)+reft) + ' months'
    plt.figtext(0.73,0.175,text, fontsize=14)    
    BA = gy.BA[reft]
    text = 'Basal area  ' +str(BA)[:4] + ' $m^{2}$ $ha^{-1}$ at age ' + str(len(gy.Hdom)+reft) + ' months'
    plt.figtext(0.73,0.15,text, fontsize=14)    
    dbh = gy.WMeanDiam[reft]
    text = 'Weigthed mean diameter  ' +str(dbh)[:4] + ' cm at age ' + str(len(gy.Hdom)+reft) + ' months'
    plt.figtext(0.73,0.125,text, fontsize=14)    

    plt.show()    


def fDrawTreeNutBal(N,P,K, gy, de):
    fs=14
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)    
    dt = gy.g['dt'][0]; dty = dt/12.0
    time = np.arange(0, len(N.stemNet))*dty
    fig= plt.figure(num = 'Tree nutrient balance', facecolor=fcol, edgecolor='k',figsize=figsi)   #Figsize(w,h), tuple inches 
    
    ax1 = fig.add_axes([0.07, 0.50, 0.25, 0.46]) #left, bottom, width, height
    ax1.set_title('Tree N balance', fontsize=fs)    
    #ax1.set_xlabel('Age, yrs', fontsize = 14)
    ax1.set_ylabel('kg $ha^{-1}$', fontsize=fs)
    ax1.plot(time, N.treeDemand, 'r-', label= 'demand'); ax1.plot(time, N.treeSupply, 'g-', label='supply')
    ax1.fill_between(time, N.treeSupply, N.treeDemand, where=  N.treeSupply< N.treeDemand, facecolor = 'r', alpha=0.3)    
    ax1.fill_between(time,N.treeDemand, N.treeSupply,  where=  N.treeSupply> N.treeDemand, facecolor = 'g', alpha=0.3)    
    ax1.legend(loc='lower right', fontsize=fs-1, ncol=1)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    


    ax2 = fig.add_axes([0.07, 0.08, 0.25, 0.35]) #left, bottom, width, height
    ax2.set_title('Cumulative tree N balance', fontsize=fs)    
    ax2.set_xlabel('Age, yrs', fontsize = fs)
    ax2.set_ylabel('kg $ha^{-1}$', fontsize=fs)
    ax2.plot(time, np.cumsum(N.treeDemand), 'r-', label='demand'); ax2.plot(time, np.cumsum(N.treeSupply), 'g-', label='supply')
    ax2.fill_between(time, np.cumsum(N.treeSupply), np.cumsum(N.treeDemand), \
            where=  np.cumsum(N.treeSupply)< np.cumsum(N.treeDemand), facecolor = 'r', alpha=0.3)    
    ax2.fill_between(time, np.cumsum(N.treeDemand), np.cumsum(N.treeSupply), \
            where=  np.cumsum(N.treeSupply)> np.cumsum(N.treeDemand), facecolor = 'g', alpha=0.3)    
    ax2.legend(loc='lower right', fontsize=fs-1, ncol=1)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    
    
    ax3 = fig.add_axes([0.40, 0.50, 0.25, 0.46]) #left, bottom, width, height
    ax3.set_title('Tree P balance', fontsize=fs)    
    #ax3.set_xlabel('Age, yrs', fontsize = 14)
    ax3.set_ylabel('kg $ha^{-1}$', fontsize=fs)
    ax3.plot(time, P.treeDemand, 'r-', label='demand'); ax3.plot(time, P.treeSupply, 'g-', label ='supply')
    ax3.fill_between(time, P.treeSupply, P.treeDemand, where=  P.treeSupply< P.treeDemand, facecolor = 'r', alpha=0.3)    
    ax3.fill_between(time, P.treeDemand, P.treeSupply, where=  P.treeSupply> P.treeDemand, facecolor = 'g', alpha=0.3)    
    ax3.legend(loc='lower right', fontsize=fs-1)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    
 
    ax4 = fig.add_axes([0.40, 0.08, 0.25, 0.35]) #left, bottom, width, height
    ax4.set_title('Cumulative tree P balance', fontsize=fs)    
    ax4.set_xlabel('Age, yrs', fontsize = fs)
    ax4.set_ylabel('kg $ha^{-1}$', fontsize=fs)
    ax4.plot(time, np.cumsum(P.treeDemand), 'r-', label='demand'); ax4.plot(time, np.cumsum(P.treeSupply), 'g-', label='supply')
    ax4.fill_between(time, np.cumsum(P.treeSupply), np.cumsum(P.treeDemand), \
            where=  np.cumsum(P.treeSupply)< np.cumsum(P.treeDemand), facecolor = 'r', alpha=0.3)    
    ax4.fill_between(time, np.cumsum(P.treeDemand), np.cumsum(P.treeSupply), \
            where=  np.cumsum(P.treeSupply)> np.cumsum(P.treeDemand), facecolor = 'g', alpha=0.3)    
    ax4.legend(loc='lower right', fontsize=fs-1, ncol=1)    
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
        
    ax5 = fig.add_axes([0.72, 0.50, 0.25, 0.46]) #left, bottom, width, height
    ax5.set_title('Tree K balance', fontsize=fs)    
    #ax5.set_xlabel('Age, yrs', fontsize = 14)
    ax5.set_ylabel('kg $ha^{-1}$ $timestep^{-1}$', fontsize=fs)
    ax5.plot(time, K.treeDemand, 'r-',label='demand'); ax5.plot(time, K.treeSupply, 'g-', label='supply')
    ax5.fill_between(time, K.treeSupply, K.treeDemand, where=  K.treeSupply< K.treeDemand, facecolor = 'r', alpha=0.3)    
    ax5.fill_between(time, K.treeDemand, K.treeSupply, where=  K.treeSupply> K.treeDemand, facecolor = 'g', alpha=0.3)    
    ax5.legend(loc='lower right', fontsize=fs-1, ncol=1)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)

    ax6 = fig.add_axes([0.72, 0.08, 0.25, 0.35]) #left, bottom, width, height
    ax6.set_title('Cumulative tree K balance', fontsize=fs)    
    ax6.set_xlabel('Age, yrs', fontsize = fs)
    ax6.set_ylabel('kg $ha^{-1}$ ', fontsize=fs)
    ax6.plot(time, np.cumsum(K.treeDemand), 'r-', label='demand');ax6.plot(time, np.cumsum(K.treeSupply), 'g-', label='supply')
    ax6.fill_between(time, np.cumsum(K.treeSupply), np.cumsum(K.treeDemand), \
            where=  np.cumsum(K.treeSupply)< np.cumsum(K.treeDemand), facecolor = 'r', alpha=0.3)    
    ax6.fill_between(time, np.cumsum(K.treeDemand), np.cumsum(K.treeSupply), \
            where=  np.cumsum(K.treeSupply)> np.cumsum(K.treeDemand), facecolor = 'g', alpha=0.3)    
    ax6.legend(loc='lower right', fontsize=fs-1, ncol=1)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)

    plt.show()
