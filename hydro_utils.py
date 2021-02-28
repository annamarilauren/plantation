# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 16:48:25 2017

@author: lauren
"""
import pandas as pd
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as interS
from scipy.interpolate import interp1d


def getParams(ParamFile=None, shname='mesic', lrs=60):    
    
    #-----Profile parameters --------------
    #self.nLyrs = 40 
    profPara={'nLyrs':lrs,  'profname':'testi' }
    nLyrs = profPara['nLyrs']        
    #------Parameters from file----------
    pFdf = pd.read_excel(ParamFile, sheet_name='Retention', skiprows=0)          #get database
    pFdf=pFdf.set_index('Number')                                               #index by soil number
    allLyrs =pd.read_excel(ParamFile, sheet_name=shname, skiprows=0)          #get the soil codes in the profie
    dz = allLyrs['dz']; lyrs = allLyrs['Lyrs']       
    pF = pFdf.iloc[list(lyrs)].copy()                                             #filter with the list
    pF['lyr']=range(len(pF))                                                    #add layer number to the df
    pF = pF.set_index('lyr')                                                    #set it to index
    pF['dz']=dz                                                                 #thickness of layer (m)  
    z = np.cumsum(dz) - dz/2.0                                                  #depth of node center (m)     
    pF['z']= z
    pF = pF[:nLyrs].to_dict()                                                   #cut the unnecessary layer away
    d={}
    d['profPara']=profPara; d['pF']=pF;     
    return d      

def getRainfall(rainFile='C:\Apps\WinPython-64bit-2.7.10.3\IPEWG\hydro\\rainfall.csv'):
    df=pd.read_csv(rainFile, names=['Date', 'P mm'], skiprows=1)
    df['Date']= pd.to_datetime(df['Date'])
    df.index= df['Date']
    del df['Date']     
    return df    

def getET(start_doy, length):
    fre=0.0163218562365; amp=0.136089550465; pha=-1.58636102519; offs=3.46707833369
    x = np.array(range(start_doy, length))   
    return np.sin(x * fre + pha) * amp + offs


def wrc(pF, x=None, var=None):
    """
    vanGenuchten-Mualem soil water retention curve\n
    IN:
        pF - dict['ThetaS': ,'ThetaR': ,'alpha':, 'n':,] OR
           - list [ThetaS, ThetaR, alpha, n]
        x  - soil water tension [m H2O = 0.1 kPa]
           - volumetric water content [vol/vol]
        var-'Th' is x=vol. wat. cont.
    OUT:
        res - Theta(Psii) or Psii(Theta)
    NOTE:\n
        sole input 'pF' draws water retention curve and returns 'None'. For drawing give only one pF-parameter set. 
        if several pF-curves are given, x can be scalar or len(x)=len(pF). In former case var is pF(x), in latter var[i]=pf[i,x[i]]
               
    Samuli Launiainen, Luke 2/2016
    """

    if type(pF) is dict: #dict input
        #Ts, Tr, alfa, n =pF['ThetaS'], pF['ThetaR'], pF['alpha'], pF['n']
        Ts=np.array(list(pF['ThetaS'].values())); 
        Tr=np.array(list(pF['ThetaR'].values())); 
        alfa=np.array(list(pF['alpha'].values())); 
        n=np.array(list( pF['n'].values()))
        m= 1.0 -np.divide(1.0,n)
       
    else: #list input
        pF=np.array(pF, ndmin=1) #ndmin=1 needed for indexing to work for 0-dim arrays
        Ts=pF[0]; Tr=pF[1]; alfa=pF[2]; n=pF[3] 
        
        m=1.0 - np.divide(1.0,n)
    
        
    def theta_psi(x): #'Theta-->Psi'
        x=np.minimum(x,Ts) 
        x=np.maximum(x,Tr) #checks limits
        s= ((Ts - Tr) / (x - Tr))#**(1/m)
        Psi=-1e-2/ alfa*(s**(1/m)-1)**(1/n) # in m
        return Psi
        
    def psi_theta(x): # 'Psi-->Theta'
        x=100*np.minimum(x,0) #cm
        Th = Tr + (Ts-Tr)/(1+abs(alfa*x)**n)**m
        return Th     
         
    # This does all the work    
    if x is None and np.size(Ts)==1:  #draws pf-curve
        import matplotlib.pylab as plt
        xx=-np.logspace(-4,5,100) #cm
        yy=psi_theta(xx)
        #field capacity and wilting point
        fc=psi_theta(-1.0) 
        wp=psi_theta(-150.0)
        fig=plt.figure()
        fig.suptitle('vanGenuchten-Mualem WRC', fontsize=16)
        #ttext=str(pF).translate(None,"{}'")
        ttext= r'$\theta_s=$'+ str(Ts)+ r', $\theta_r=$' +str(Tr)+ r', $\alpha=$'+str(alfa)+ ',n='+str(n)

        plt.title(ttext, fontsize=14)
        plt.semilogx(-xx,yy,'g-')
        plt.semilogx(1, fc,'ro', 150, wp, 'ro')#fc, wp
        plt.text(1,1.1*fc,'FC'), plt.text(150, 1.2*wp,'WP')
        plt.ylabel(r'$\theta$  $(m^3m^{-3})$', fontsize=14) 
        plt.xlabel('$\psi$ $(cm)$', fontsize=14)
        plt.ylim(0.8*Tr, min(1,1.1*Ts))
        del xx, yy
        return None
    elif x is None: 
        print ('wrc: To draw water-retention curve give only one pF -parameter set')
        return None
        
    if var is 'Th': y=theta_psi(x) #'Theta-->Psi'
                
    else: y=psi_theta(x) # 'Psi-->Theta'
               
    return y


def CWTr(profPara, pF, direction='positive'):
    #-------Parameters ---------------------
    nLyrs = profPara['nLyrs']
    z = np.array(pF['z'].values())   
    dz =np.array(pF['dz'].values())
    
    #--------- Connection between gwl and water storage------------
    d = 6 if direction == 'positive' else -6   
    gwl=np.linspace(0,d,150)

    if direction == 'positive':
        sto = [sum(wrc(pF, x = np.minimum(z-g, 0.0))*dz) for g in gwl]     #equilibrium head m
    else:
        sto = [sum(wrc(pF, x = np.minimum(z+g, 0.0))*dz) for g in gwl]     #equilibrium head m
    gwlToSto = interp1d(np.array(gwl), np.array(sto), fill_value='extrapolate')
    sto = list(sto); gwl= list(gwl)        
    sto.reverse(); gwl.reverse()
    stoToGwl =interp1d(np.array(sto), np.array(gwl), fill_value='extrapolate')
    del gwl, sto
    
    
    #----------Transmissivity-------------------
    K=np.array(pF['KsatHorizontal'].values())/100.0*24.0   #from cm/h to m/day
    tr =[sum(K[t:]*dz[t:]) for t in range(nLyrs)]        
    if direction=='positive':        
        gwlToTra = interS(z, np.array(tr))            
    else:
        z= list(z);  z.reverse(); tr.reverse()
        gwlToTra = interS(-np.array(z), np.array(tr))            
        
    del tr
            
    return gwlToSto, stoToGwl, gwlToTra
        
def adjacent_mean(a, ny, nx, drmask):
    """
    computes mean phi (response variable) in four adjacent cells (north, south, east, west), excluding ditch nodes.        
    Input: \n
        a in flattened array \n 
        ny, nx dimensions of the original 2D array \n
        drmask drain mask in flattened array \n       
    Output:  \n
        mean of variable 'a' surrounding the drain mask cells, excluding the drain cells
    """    
    aa=np.empty((ny+2,nx+2))
    aa[:,:]=np.nan; aa[1:-1,1:-1]=np.reshape(a, (ny,nx))
    ya = np.array(drmask)/nx + 1
    xa=np.mod(np.array(drmask),nx) +1
    aa[ya,xa] = np.nan
    w=aa[ya,xa-1]
    e=aa[ya,xa+1]
    n=aa[ya-1,xa]
    s=aa[ya+1,xa]
    return np.nanmean([w,e,n,s], axis=0)
        


#getRaifall()
        
def writeRainToP():
    from calendar import monthrange
    from calendar import isleap
    folder ='C:\Apps\WinPython-64bit-2.7.10.3\IPEWG\hydro\\'
    xf = folder + 'PKU_rainfall.xls'
        
    p = pd.date_range('1994-01-01', '2016-12-31')
    prec = []
    for yr in range(1994,2017):
        if yr==1998 or yr==2005 or yr==2006:
            dlen = 365 if isleap(yr)==False else 366
            prec = prec+ dlen*list([np.nan])
        else:
            df0 = pd.read_excel(xf, skiprows=9, sheet_name=str(yr)[2:], index_col=None, parse_cols='A:M', header=None)
            df0=df0[0:31]
            df0=df0.replace('TTU', np.nan)
            df0=df0.replace('-', np.nan)
            for m in range(1,13):
                _,d = monthrange(yr,m)
                print (yr, m,d)
                df00 = df0[:d][m]
                mp=np.nanmean(df00)
                df00=df00.replace(np.nan, mp)
                prec = prec + list(df00.values)
    parr=np.array(prec)
    idx = np.where(np.isnan(parr))
    print (idx)
    indday= np.mod(idx,365)
    dfp = pd.DataFrame(np.array(prec), index= p)
    d00=dfp.groupby([dfp.index.dayofyear]).mean()
    a00= np.ravel(d00.values)    
    parr[idx]= a00[indday]            
    dfp = pd.DataFrame(parr, index= p)
    dfp.to_csv(folder + 'rainfall.csv', sep=',')
    ysum=dfp.resample('A', how='sum')    
    print (ysum)
#writeRainToP()

def compET():
    from calendar import monthrange
    from calendar import isleap
    folder ='C:\Apps\WinPython-64bit-2.7.10.3\IPEWG\hydro\\'
    xf = folder + 'PenmannMonteithFAO.xlsm'        
    p = pd.date_range('2007-02-01', '2012-08-31')
    df0 = pd.read_excel(xf, skiprows=0, sheet_name='ET', index_col=0, parse_cols='A:B')
    mET= df0.resample('M', how='mean')        
    mET.to_csv(folder + 'monthlyET.csv', separator=',')     

def my_sin(x, freq, amplitude, phase, offset):
    # create the function we want to fit
    return np.sin(x * freq + phase) * amplitude + offset

def fit_sinET():
    from scipy.optimize import curve_fit
    import matplotlib.pylab as plt
    #Curve fit here
    optFit = True
    if optFit == True:        
        et = [3.35, 3.38, 3.32, 3.47, 3.56, 3.58, 3.53, 3.68, 3.44, 3.62, 3.4,3.31]; et = np.array(et)
        months = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 310, 365]; months= np.array(months)
        guess_freq = 1/365.0
        guess_amplitude = 3*np.std(et)/(2**0.5)
        guess_phase = 12/(2.0*np.pi)
        guess_offset = np.mean(et)
        
        p0=[guess_freq, guess_amplitude,
            guess_phase, guess_offset]
        
        # now do the fit
        #my_sin = lambda time, *p0: np.sin(time * p0[0] + p0[2]) * p0[1] + p0[3]
              
        fit = curve_fit(my_sin, months, et, p0=p0)   
        
        fre= fit[0][0]
        amp = fit[0][1]
        pha = fit[0][2]
        offs = fit[0][3] 
        print (fre, amp, pha, offs)        
    

    x=np.array(range(365));
    y=np.sin(x * fre + pha) * amp + offs
    fig= plt.figure(num = 'Peat', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(25.0,12.0))   #Figsize(w,h), tuple inches 
    ax1 = fig.add_axes([0.33, 0.47, 0.3, 0.46]) #left, bottom, width, height
    ax1.set_title('ET')    
    ax1.set_xlabel('x', fontsize = 14)
    ax1.plot(months, et, 'ko')

    ax1.plot(x, y, 'r-')
    plt.show()

    return fre, amp, pha,offs

def getET(start_doy, length):
    fre=0.0163218562365; amp=0.136089550465; pha=-1.58636102519; offs=3.46707833369
    x = np.array(range(start_doy, length))   
    return np.sin(x * fre + pha) * amp + offs

#fit_sinET()

def gmeanTr(Tr):
    """
    Input: 
        Transmissivity vector, tr in node center point
    Output:
        Transmissivity, tr in west surface sqrt(Tr(i-1)*Tr(i)) and east sqrt(Tr(i)*Tr(i+1)) 
    """
    n = len(Tr)
    Trwest = np.sqrt(Tr[:n-1]*Tr[1:]); Trwest=np.append(Trwest, 0.0)
    Treast = np.sqrt(Tr[1:]*Tr[:n-1]); Treast=np.insert(Treast, 0, 0.0)                
    return Trwest, Treast

def Hadjacent(H):
    """
    Input:
        H vector, H in each node
    Output:
        Hwest H(i-1), Heast H(i+1)
    """
    n=len(H)
    Hwest = H[0:n-1]; Hwest=np.append(Hwest, 0.0)
    Heast = H[1:]; Heast=np.insert(Heast, 0, 0.0)
    return Hwest, Heast  

def Amatrix(A, n, implic, Trwest, Treast, alfa):
    """
    Construction of tridiagonal matrix
    """     
    i,j = np.indices(A.shape)
    A[i==j]= implic*(Trwest+Treast)+alfa                                 #diagonal element
    A[i==j+1]=-implic*Trwest[:n-1]                                       #West element
    A[i==j-1]=-implic*Treast[1:]                                         #East element     

    return A

def boundConst(A, n):
    """
    Diriclet (constant head boundary conditions)
    """    
    A[0,0]=1; A[0,1]=0.                                             #Dirichlet, west boundary
    A[n-1,n-1]=1.; A[n-1, n-2]=0.                                   #Dirichlet, east boundary
    return A

def boundNoFlow(A, n, implic, Trwest, Treast, alfa):
    """
    Diriclet (constant head boundary conditions)
    """    
                                                #Dirichlet, west boundary
    A[0,0]= 2.*implic*(Treast[0])+alfa[0]                                 #diagonal element
    A[0,1]=-2*implic*Trwest[0]                                         #East element     
    A[n-1,n-1]=2.*implic*(Trwest[n-1])+alfa[0]
    A[n-1, n-2]=-2*implic*Treast[n-1]                                       #West element
    return A

def boundConstW(A, n):
    """
    Diriclet (constant head boundary conditions)
    """    
    A[0,0]=1; A[0,1]=0.                                             #Dirichlet, west boundary
    A[n-1,n-1]=1.; A[n-1, n-2]=0.                                   #Dirichlet, east boundary
    return A

def boundNoFlowE(A, n, implic, Trwest, Treast, alfa):
    """
    Neumann 
    """    
                                                                            
    A[0,0]= 2.*implic*(Treast[0])+alfa[0]                                 #diagonal element
    A[0,1]=-2*implic*Trwest[0]                                         #East element     
    A[n-1,n-1]=2.*implic*(Trwest[n-1])+alfa[0]
    A[n-1, n-2]=-2*implic*Treast[n-1]                                       #West element
    return A

def rightSide(S,dt,dy, implic,alfa, H, Trminus0,Hminus,Trplus0,Hplus, DrIrr, Htmp1, ele, h0):
    hs = S*dt*dy**2 + alfa*H + (1-implic)*(Trminus0*Hminus) - (1-implic)*(Trminus0 + Trplus0)*H  + (1-implic)*(Trplus0*Hplus)
    n=len(Htmp1)
            
    if DrIrr==False:                
        hs[0]=Htmp1[1] if Htmp1[0]>Htmp1[1] else min(ele[0]+h0, Htmp1[1])
        hs[n-1]=Htmp1[n-2] if Htmp1[n-1]>Htmp1[n-2] else min(ele[n-1]+h0, Htmp1[n-2])    #if wt below canal water level, lower the canal wl to prevent water inflow to compartment
    else:
        hs[0]=ele[0]+h0
        hs[n-1]=ele[n-1]+h0 
    return hs


def rightSide2(S,dt,dy, implic,alfa, H, Trminus0,Hminus,Trplus0,Hplus, DrIrr, Htmp1, ele, h0, hup):
    hs = S*dt*dy**2 + alfa*H + (1-implic)*(Trminus0*Hminus) - (1-implic)*(Trminus0 + Trplus0)*H  + (1-implic)*(Trplus0*Hplus)
    n=len(Htmp1)
            
    if DrIrr==False:                
        hs[0]=Htmp1[1] if Htmp1[0]>Htmp1[1] else min(ele[0]+h0, Htmp1[1])
        hs[n-1]=Htmp1[n-2] if Htmp1[n-1]>Htmp1[n-2] else min(ele[n-1]+hup, Htmp1[n-2])    #if wt below canal water level, lower the canal wl to prevent water inflow to compartment
    else:
        hs[0]=ele[0]+h0
        hs[n-1]=ele[n-1]+hup 
    return hs

def getOnsiteData():
    #get TPK 10 results for testing, rainfall + water tables        
    sfolder ='C:\Apps\WinPython-64bit-2.7.10.3\IPEWG\hydro\\'    
    dfP = pd.read_csv(sfolder+'onsite.csv') #, index_col=0)
    P = dfP['Rainfall'].values
    dates=pd.to_datetime(dfP['DateTime'],format='%Y-%m-%d %H:%M:%S')
    dfP['dates'] = pd.Series(dates); dfP.index= dfP['dates']; dfP=dfP.drop('DateTime',1)    
    dfP=dfP.drop('dates',1)
    dfG = pd.read_csv(sfolder+'onsiteGate.csv') #, index_col=0)
    return P, dates, dfG
        
