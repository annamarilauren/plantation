# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 18:25:51 2016

@author: Ari Lauren
Computes decomposition of ogranic matter:
Litter, Coarse woody debris and Peat
Litter decomposition, CO2 release, N dynamics:\n
    -Chertov, O. G., Komarov, A. S., Nadporozhskaya, M., Bykhovets, S. S., Zudin, S. L., 2001. ROMUL – a model of forest soil organic matter dynamics as a substantial tool for forest ecosystem modeling. Ecol. Model. 138, 289-308.
CWD decomposition and release of CO2: \n
    -MacKenzen, J, Bauhus, J & Webber, E. 2003. Decomposition rates of coarse woody debis - A review with particular emphasis on Australian tree species. Australian Journal of Botany 51: 27-37.
CO2 efflux from tropical peat and conscutive decompositin: \n
    -Hooijer, A., Page, S., Canadell, J. G., Silvius, M., Kwadijk, J., Wösten, H., and Jauhiainen, J. 2010. Current and future CO2 emissions from drained peatlands in
        Southeast Asia, Biogeosciences 7, 1505-1514.
    -Jauhiainen, J., Hooijer, A., and PageS. E. 2012.Carbon dioxide emissions from an Acacia plantation on peatland in
        Sumatra, Indonesia, Biogeosciences, 9, 617–630
Nutrient storages in cwd and the 
"""
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

class DecomCWD():
    def __init__(self, gen, dpara, gy, N, P, K):                               # general parameters, decomposition parameters, growth and yield instance
        self.g = gen                                                           # General parameters
        self.cwdSto = np.empty(0)                                              # Storage of coarse woody debris kg ha-1
        self.cwdMloss = np.empty(0)                                            # Mass loss of CWD, kg ha-1 timestep-1
        self.cwdCInp = np.empty(0)                                             # Carbon input as CWD into soil, kg ha-1 timestep-1
        self.cwdCSto = np.empty(0)                                             # Carbon storage in CWD kg ha-1
        self.cwdCRelease = np.empty(0)                                         # Carbon emission due to CWD decay kg ha-1
        self.massToC = 0.5                                                     # C content of CWD kg kg-1
        self.dpara=dpara; rg = 0                                               # Initial values and row in input file
        rg = 0 ; rd = 0; rn = 0                                                # Row in general parameters, in decomposition parameters, in nutrient parameters
        to = dpara['soilTopen'][rd]                                            # Soil temerature open area in deg C       
        tc = dpara['soilTcanopy'][rd]                                          # Soil temerature under canopy in deg C       
        mgr = self.g['maxGreenMass'][rg]                                       # maximum green mass includes above ground weed mass and tree stand foliage mass              
        gr = gy.BiFoliage + gy.weedAbove                                       # existing green biomass kg ha-1        
        self.soilT = to - (gr / mgr)*(to - tc)

        self.N = N['concWood'][rn]                                             # N content of wood (%)
        self.P = P['concWood'][rn]                                             # P content of wood (%)
        self.K = K['concWood'][rn]                                             # K content of wood (%)       
        self.cwdNutSto ={}        
        self.cwdNutRelease ={}        
        
    def testSet(self, cwdIn, rowGe, rot=None, dbh = 15.0):
        """
        Run CWD model with example set.
        """
        length= len(cwdIn)
        cwdIn = np.zeros(length); cwdIn[0]=100.0; tm=self.soilT; dbh = np.zeros(length); dbh[0]=8.0        
        #print cwdIn
        #print dbh
        self.decomCWD(cwdIn, rowGe, rot = rot, dbh = dbh, tm=32)
        x = range(length)
        fig= plt.figure(num = 'CWD', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(25.0,12.0))   #Figsize(w,h), tuple inches 
        ax1 = fig.add_axes([0.33, 0.47, 0.3, 0.46]) #left, bottom, width, height
        ax1.set_title('CWD model')    
        ax1.set_xlabel('time, months', fontsize = 14)
        ax1.set_ylabel('CWD mass remaining', fontsize = 14)
        ax1.plot(x, self.cwdSto, 'r-')
        plt.show()

    def decomCWD(self, cwdIn, rowGe, rot=None, dbh = None, tm = None):
        """
        Computes storage and decompostion on coarse woody debris. 
        Ref: MacKenzen, J, Bauhus, J & Webber, E. 2003. Decomposition rates of coarse woody debis - A review with particular emphasis on Australian tree species. Australian Journal of Botany 51: 27-37.
        Includes two models: \n
            a. f(temperature), Table 2, Eq. 1  \n
            b. f(temperature, diameter, wood density), Eq 7 \n   
            c. f(temperature given as array, computed T(stand green mass)) \n
            d. f(temperature as array, diameter, wood density)  \n
        Arguments: \n
            cwdIn time series of cwd input to soil (kg/ha/timestep) \n
        Kwargs: \n
            dbh - time series of stem mean diameter (in cm), if given computes with model b, if None uses model a \n 
            rot - array telling in which rotation the stand is, used to pick correct parameters for the computation \n
            tm - annual mean temperature (deg C), array  \n
        Output: \n
            Storage of CWD in kg/ha (mass) \n
            Storage of CWD in C kg/ha \n
        """
        
        from scipy.misc import derivative
        from scipy.interpolate import InterpolatedUnivariateSpline as interS
        
        #--------parameters--------------
        dt = self.g['dt'][0]; length = len(cwdIn); dty = dt/12.0        
        time = np.arange(0, length)*dty; 
        N = self.N; P = self.P; K = self.K               
        
        cwdSto = np.zeros(length); cwdMloss = np.zeros(length) 
        cwdNsto = np.zeros(length); cwdNrelease = np.zeros(length)
        cwdPsto = np.zeros(length); cwdPrelease = np.zeros(length)
        cwdKsto = np.zeros(length); cwdKrelease = np.zeros(length)
        
        kAsConst = bool
        if tm is None: 
            tm = self.soilT
            kAsConst = False                
        #------data points for nutrient release model (Palviainen & Finér 2015. Decomposition and nutrient release from Norway spruce coarse roots
                                        #and stumps – A 40-year chronosequence study. Forest Ecology and Management 358 (2015) 1–11)                               
        #-----------N, P, K storages and release as a function of remaining mass
        m_arr = np.array([0, 0.1, 0.2, 0.25, 1.0])
        n_arr = np.array([0, 0.56, 0.9, 1.0, 1.0])
        fn = interS(m_arr,n_arr,k=1) 
        m_arr = np.array([0, 0.06, 0.2, 0.25, 0.81, 1.0])
        p_arr = np.array([0, 0.20, 0.44, 0.48, 0.60, 1.0])
        fp = interS(m_arr,p_arr,k=1)
        m_arr = np.array([0, 0.2, 0.40, 0.60, 0.81, 1.0])
        k_arr = np.array([0, 0.18, 0.30, 0.35, 0.40, 1.0])
        fk = interS(m_arr,k_arr,k=1)
        
        #-------Model computation-----------
        #kAsConst = True; dbh = None; tm = 32. 
        if kAsConst==True: 
            # -----constant temperature is given
            if dbh is None: 
                #---- only mean temperature as explanatory
                k = 0.0166 * np.exp(0.093 * tm) 
                fSto = lambda tim: np.exp(-k*tim)
                fRelease = lambda tim: (derivative(fSto, tim))
                for n in range(length):
                    ts = time[n:]-time[n]
                    s = cwdIn[n]* fSto(ts)
                    r = cwdIn[n]* fRelease(ts)*dty                              # storage and release from this cohort        
                    cwdSto[n:] = cwdSto[n:] + s
                    cwdMloss[n:] = cwdMloss[n:] + r                             #Mass loss in time step
                    #-------------nutrients--------------------------------                    
                    cwdNsto[n:] = cwdNsto[n:] + N/100.0 * cwdIn[n] * fn(fSto(ts))  #Nitrogen storage, release
                    cwdNrelease[n+1:] = cwdNrelease[n+1:] + np.diff( N/100.0 * cwdIn[n] * fn(fSto(ts)))*-1.0             #Mass loss in time step
                    cwdPsto[n:] = cwdPsto[n:] + P/100.0 * cwdIn[n] * fp(fSto(ts))  #Phosphorus storage, release
                    cwdPrelease[n+1:] = cwdPrelease[n+1:] + np.diff( P/100.0 * cwdIn[n] * fp(fSto(ts)))*-1.0             #Mass loss in time step
                    cwdKsto[n:] = cwdKsto[n:] + K/100.0 * cwdIn[n] * fk(fSto(ts))  #Potassium storage, release
                    cwdKrelease[n+1:] = cwdKrelease[n+1:] + np.diff( K/100.0 * cwdIn[n] * fk(fSto(ts)))*-1.0             #Mass loss in time step
                
            else:
                # -----CWD properties: mean diameter and wood density,
                rg = [rowGe[int(r)] for r in rot]                                           #picks the correct row according to the ongoing rotation
                X0 = np.array([self.g['WoodDens'][r] for r in rg])
                k = (0.07816 + 0.010413 * tm - 0.002012 * dbh - 0.16749 * X0/1000.0)*2.0   #calibrated to measurements in April
                for n in range(length):
                    ts = time[n:]-time[n]; kstar = k[n:]
                    cwdSto[n:] = cwdSto[n:] + cwdIn[n]*np.exp(-kstar*ts)
                    cwdMloss[n+1:] = cwdMloss[n+1:] + np.diff(cwdIn[n]*np.exp(-kstar*ts))             #Mass loss in time step
                    #---------nutrients------------------------------
                    cwdNsto[n:] = cwdNsto[n:] + N/100.0 * cwdIn[n] * fn(np.exp(-kstar*ts))  #Nitrogen storage, release
                    cwdNrelease[n+1:] = cwdNrelease[n+1:] + np.diff( N/100.0 * cwdIn[n] * fn(np.exp(-kstar*ts)))*-1.0             #Mass loss in time step
                    cwdPsto[n:] = cwdPsto[n:] + P/100.0 * cwdIn[n] * fp(np.exp(-kstar*ts))  #Phosphorus storage, release
                    cwdPrelease[n+1:] = cwdPrelease[n+1:] + np.diff( P/100.0 * cwdIn[n] * fp(np.exp(-kstar*ts)))*-1.0             #Mass loss in time step
                    cwdKsto[n:] = cwdKsto[n:] + K/100.0 * cwdIn[n] * fk(np.exp(-kstar*ts))  #Potassium storage, release
                    cwdKrelease[n+1:] = cwdKrelease[n+1:] + np.diff( K/100.0 * cwdIn[n] * fk(np.exp(-kstar*ts)))*-1.0             #Mass loss in time step
        
        else:                                                                   #if k is an array with length of number of time steps
            #------temperature is given as array
            if dbh is None:            
                #-------Only temperature array is given as explanatory------                
                k = 0.0166 * np.exp(0.093 * tm)
                for n in range(length):
                    s=cwdIn[n]
                    cwdSto[n] = cwdSto[n] + s
                    nn = s*N/100.0; pp = s*P/100.0; kk = s*K/100.0
                    for m in range(n+1,length):
                        ds = s*(1.0 - np.exp(-k[m]*dty))                        
                        s =  s * np.exp(-k[m]*dty)
                        cwdSto[m] = cwdSto[m] + s                    
                        cwdMloss[m] = cwdMloss[m] + ds
                        """
                        #------- nutrients------------                       
                        dn = N/100.0 * cwdIn[n] * (fn(np.exp(-k[m]*m*dty))- fn(np.exp(-k[m]*(m+1)*dty)))
                        nn = nn-dn                        
                        cwdNsto[m] = cwdNsto[m] + nn  #Nitrogen storage, release
                        cwdNrelease[m] = cwdNrelease[m] + dn
                        dp = P/100.0 * cwdIn[n] * (fp(np.exp(-k[m]*m*dty))- fp(np.exp(-k[m]*(m+1)*dty)))
                        pp = pp-dp
                        cwdPsto[m] = cwdPsto[m] + pp  #Nitrogen storage, release
                        cwdPrelease[m] = cwdPrelease[m] + dp
                        dk = P/100.0 * cwdIn[n] * (fk(np.exp(-k[m]*m*dty))- fk(np.exp(-k[m]*(m+1)*dty)))
                        kk = kk - dk 
                        cwdKsto[m] = cwdKsto[m] + kk  #Nitrogen storage, release
                        cwdKrelease[m] = cwdKrelease[m] + dk
                        """
            else:
                # -----CWD properties: mean diameter as array and wood density, temperature as array                 
                rg = [rowGe[int(r)] for r in rot]                                           #picks the correct row according to the ongoing rotation
                dbh[0] = dbh[-1]                
                X0 = np.array([self.g['WoodDens'][r] for r in rg])
                k = (0.07816 + 0.010413 * tm - 0.002012 * dbh - 0.16749 * X0/1000.0)*2.0   #calibrated to measurements in April
                for n in range(length):
                    s=cwdIn[n]; s0=s 
                    cwdSto[n] = cwdSto[n] + s
                    nn = s*N/100.0; pp = s*P/100.0; kk = s*K/100.0
                    for m in range(n+1,length):  #
                        ds = s*(1.0 - np.exp(-k[m]*dty))                        
                        s =  s * np.exp(-k[m]*dty)
                        cwdSto[m] = cwdSto[m] + s                    
                        cwdMloss[m] = cwdMloss[m] + ds
                        #------- nutrients------------   
                        pro = 0 if cwdIn[n] == 0 else s/cwdIn[n]
                        pro0 = 0 if cwdIn[n] == 0 else s0/cwdIn[n]                    
                        dn = N/100.0 * cwdIn[n] * (fn(pro0) - fn(pro))
                        nn = nn-dn                        
                        cwdNsto[m] = cwdNsto[m] + nn  #Nitrogen storage, release
                        cwdNrelease[m] = cwdNrelease[m] + dn
                        dp = P/100.0 * cwdIn[n] * (fp(pro0) - fp(pro))
                        pp = pp-dp
                        cwdPsto[m] = cwdPsto[m] + pp  #Nitrogen storage, release
                        cwdPrelease[m] = cwdPrelease[m] + dp
                        dk = K/100.0 * cwdIn[n] * (fk(pro0) - fk(pro))
                        kk = kk - dk 
                        cwdKsto[m] = cwdKsto[m] + kk  #Nitrogen storage, release
                        cwdKrelease[m] = cwdKrelease[m] + dk
                        s0=s
        
        #------locate to instace variables
        self.cwdCInp = np.append(self.cwdCInp, cwdIn*self.massToC); self.cwdSto = np.append(self.cwdSto, cwdSto) 
        self.cwdMloss = np.append(self.cwdMloss, cwdMloss);self.cwdCSto = np.append(self.cwdCSto, cwdSto*self.massToC)
        #self.cwdCRelease = np.append(self.cwdCRelease,cwdMloss*self.massToC*-1.0)
        self.cwdCRelease = np.append(self.cwdCRelease,cwdMloss*self.massToC)
        
        self.cwdNutSto['Nitrogen'] = cwdNsto; self.cwdNutRelease['Nitrogen'] = cwdNrelease 
        self.cwdNutSto['Phosphorus'] = cwdPsto; self.cwdNutRelease['Phosphorus'] = cwdPrelease 
        self.cwdNutSto['Potassium'] = cwdKsto; self.cwdNutRelease['Potassium'] = cwdKrelease 

        del cwdIn, cwdSto, cwdMloss, time, cwdNsto, cwdNrelease, cwdPsto, cwdPrelease, cwdKsto, cwdKrelease

class Romul():
    """
    Computes mass, carbon, N, P, and K storages and fluxes in above ground and below ground organic material. \n 
    Chertov, O. G., Komarov, A. S., Nadporozhskaya, M., Bykhovets, S. S., Zudin, S. L., 2001. ROMUL – a model of forest soil organic matter dynamics as a substantial tool for forest ecosystem modeling. Ecol. Model. 138, 289-308.
    New moisture and temperature restriction models: Romul manual (Romul description, www address ), where soil \n
    moisture is given as relative water content W/Wfc where Wfc is field papacity moisture either in volumetric or gravimetric water content. \n
    Decomposition rate is controlled by: \n
        - litter quality (N  content, ash content, lignin content) \n
        - environmental conditions (soil temperature, soil water content) \n
        - soil parameters (pH) \n
        - litter location (above, below ground)
    Decomposition processes in successional stages: \n
        - L rapid decomposition, low release on N \n
        - F decomposition slows, more N release \n
        - H humus na s stabile humus \n
    Output: \n
        - storage of soil organic matter \n
        - mineralization of organic matter \n
    """
    def __init__(self, dpara, gy, de, Npara, Ppara, Kpara, test = False):
        #------------- Setting the parameters and initial values -----------------
        self.p=dpara; rg = 0                                                   # Initial values and row in input file
        rn = 0; rp=0; rk = 0                                                   # input row for N, P, K parameters        
        self.gy = gy                                                           # Litter input from growth and yield
        self.length =  len(self.gy.age)                                        # Number of timesteps
        self.dt = gy.g['dt'][0]                                                # Time step in months
        self.soilT = de.soilT #self.p['soilTopen'][rg]                                       # Soil temerature in deg C       
        self.soilW = self.p['soilW'][rg]                                       # Soil moisture, relative to field capacity       
        self.soilpH = self.p['soilpH'][rg]                                     # Soil pH
        self.massToC = 0.5                                                     # C content of litter kg kg-1        
        self.test = test                                                       # Runs the test set

        #-------Declare litter parameters --------------        
        if self.test==True:
            TestInp = np.zeros(self.length); TestInp[0]=1; te = 0
            self.litters_a = {'TestLit':{'Desc': 'Test litter, not in computation', 'NreTransL':0.5, 'Ncont': 2.4, 'Pcont' : 0.1, 'Kcont':0.4,
                               'ashCont': 2.0, 'ligninCont':0.0, 'above': True,
                               'Linput': TestInp, 'Mstorage': [], 'Nstorage': []}}
            self.litters_b = {} #{'TestLit':{'Desc': 'Test root litter, not in computation', 'NreTrans':0.5, 'Ncont': 3.0, 'ashCont': 2.0, 'ligninCont':0.0, 'above': False,
                           #'Linput': TestInp, 'Mstorage': [], 'Nstorage': []}}
        else:
            te=1
            self.litters_a={'FoLitL':{'Desc': 'Foliar litter from living trees', 
                               'NreTrans':Npara['retrans'][rn], 'PreTrans':Ppara['retrans'][rp], 'KreTrans':Kpara['retrans'][rk], 
                               'Ncont': Npara['concLeaf'][rn], 'Pcont': Ppara['concLeaf'][rp],'Kcont': Kpara['concLeaf'][rk],
                               'ashCont': 2.0, 'ligninCont':0.0, 'above': True,
                               'Linput': self.gy.FoLitL/10000.0, 'Mstorage': [], 'Nstorage': [], 'Pstorage': [],'Kstorage': []},
                     'FoLitD':{'Desc': 'Foliar litter from dead trees', 'NreTrans':0.0, 'PreTrans':0.0, 'KreTrans':0.0,
                               'Ncont': Npara['concLeaf'][rn], 'Pcont': Ppara['concLeaf'][rp],'Kcont': Kpara['concLeaf'][rk], 
                               'ashCont': 2.0, 'ligninCont':0, 'above': True,
                               'Linput': self.gy.FoLitD/10000.0, 'Mstorage': [], 'Nstorage': [],'Pstorage': [],'Kstorage': []},
                     'BaLitL':{'Desc': 'Bark litter from living trees', 
                               'NreTrans':Npara['retrans'][rn], 'PreTrans':Ppara['retrans'][rp],'KreTrans':Kpara['retrans'][rk],
                               'Ncont': Npara['concBark'][rn], 'Pcont': Ppara['concBark'][rp],'Kcont': Kpara['concBark'][rk], 
                               'ashCont': 2.0, 'ligninCont':0, 'above': True, 
                               'Linput': self.gy.BaLitL/10000.0, 'Mstorage': [], 'Nstorage': [],'Pstorage': [],'Kstorage': []},                               
                     'BaLitD':{'Desc': 'Bark litter from dead trees', 'NreTrans':0.0, 'PreTrans':0.0, 'KreTrans':0.0,
                               'Ncont': Npara['concBark'][rn], 'Pcont': Ppara['concBark'][rp],'Kcont': Kpara['concBark'][rk], 
                               'ashCont': 2.0, 'ligninCont':0, 'above': True,
                               'Linput': self.gy.BaLitD/10000.0, 'Mstorage': [], 'Nstorage': [],'Pstorage': [],'Kstorage': []},                               
                     'BrLitL':{'Desc': 'Branch litter from living trees', 
                               'NreTrans':Npara['retrans'][rn], 'PreTrans':Ppara['retrans'][rp],'KreTrans':Kpara['retrans'][rk],
                               'Ncont': Npara['concBranch'][rn], 'Pcont': Ppara['concBranch'][rp],'Kcont': Kpara['concBranch'][rk], 
                               'ashCont': 2.0, 'ligninCont':0, 'above': True, 
                               'Linput': self.gy.BrLitL/10000.0, 'Mstorage': [], 'Nstorage': [],'Pstorage': [],'Kstorage': []},                               
                     'BrLitD':{'Desc': 'Branch litter from dead trees', 'NreTrans':0.0, 'PreTrans':0.0, 'KreTrans':0.0,
                               'Ncont': Npara['concBranch'][rn], 'Pcont': Ppara['concBranch'][rp],'Kcont': Kpara['concBranch'][rk], 
                               'ashCont': 2.0, 'ligninCont':0, 'above': True, 
                               'Linput': self.gy.BrLitD/10000.0, 'Mstorage': [], 'Nstorage': [],'Pstorage': [],'Kstorage': []},
                     'weedALitL':{'Desc': 'Weed litter from above ground living plants', 
                                  'NreTrans':Npara['retransWeeds'][rn], 'PreTrans':Ppara['retransWeeds'][rp],'KreTrans':Kpara['retransWeeds'][rk],
                                  'Ncont': Npara['concWeeds'][rn], 'Pcont': Ppara['concWeeds'][rp],'Kcont': Kpara['concWeeds'][rk], 
                                  'ashCont': 2.0, 'ligninCont':0, 'above': True, 
                               'Linput': self.gy.weedALitL/10000.0, 'Mstorage': [], 'Nstorage': [],'Pstorage': [],'Kstorage': []},                               
                     'weedALitD':{'Desc': 'Weed litter from above ground dead plants', 'NreTrans':0.0, 'PreTrans':0.0,'KreTrans':0.0,
                                  'Ncont': Npara['concWeeds'][rn], 'Pcont': Ppara['concWeeds'][rp],'Kcont': Kpara['concWeeds'][rk], 
                                  'ashCont': 2.0, 'ligninCont':0, 'above': True, 
                                  'Linput': self.gy.weedALitD/10000.0, 'Mstorage': [], 'Nstorage': [],'Pstorage': [],'Kstorage': []}}

            self.litters_b ={'FineRLitL':{'Desc': 'Fine root litter from living trees', 
                                'NreTrans':Npara['retrans'][rn], 'PreTrans':Ppara['retrans'][rp], 'KreTrans':Kpara['retrans'][rk], 
                                'Ncont': Npara['concRoot'][rn], 'Pcont': Ppara['concRoot'][rp],'Kcont': Kpara['concRoot'][rk], 
                                'ashCont': 2.0, 'ligninCont':0, 'above': False,
                                'Linput': self.gy.FineRLitL/10000.0, 'Mstorage': [], 'Nstorage': [],'Pstorage': [],'Kstorage': []},
                     'weedBLitL':{'Desc': 'Weed litter from beloww ground living plants', 
                                'NreTrans':Npara['retransWeeds'][rn],'PreTrans':Ppara['retransWeeds'][rp], 'KreTrans':Kpara['retransWeeds'][rk],  
                                'Ncont': Npara['concWeeds'][rn], 'Pcont': Ppara['concWeeds'][rp],'Kcont': Kpara['concWeeds'][rk], 
                                'ashCont': 2.0, 'ligninCont':0, 'above': False,
                                'Linput': self.gy.weedBLitL/10000.0, 'Mstorage': [], 'Nstorage': [],'Pstorage': [],'Kstorage': []},
                     'weedBLitD':{'Desc': 'Weed litter from below ground dead plants', 'NreTrans':0.0,'PreTrans':0.0,'KreTrans':0.0,
                                  'Ncont': Npara['concWeeds'][rn], 'Pcont': Ppara['concWeeds'][rp],'Kcont': Kpara['concWeeds'][rk], 
                                 'ashCont': 2.0, 'ligninCont':0, 'above': False,
                                 'Linput': self.gy.weedBLitD/10000.0, 'Mstorage': [], 'Nstorage': [],'Pstorage': [],'Kstorage': [] }}
        
        
        # ---------Inputs of litter for mass balance computation -----------------------
        self.TotInp = 0.0; self.TotInpN = 0.0; self.TotInpP = 0.0; self.TotInpK = 0.0                         # Total litter and N,P,K input as scalar kg m-2
        self.totalLitter = 0.0; self.totalNLitter = 0.0; self.totalPLitter = 0.0 ;self.totalKLitter = 0.0     # Total litter in time series kg m-2 timestep-1
        for c in self.litters_a.keys():
            self.TotInp = self.TotInp + sum(self.litters_a[c]['Linput'])
            self.TotInpN = self.TotInpN + sum((self.litters_a[c]['Linput'])*self.litters_a[c]['Ncont']/100.0 * \
                    (1.0-self.litters_a[c]['NreTrans']))
            self.TotInpP = self.TotInpP + sum((self.litters_a[c]['Linput'])*self.litters_a[c]['Pcont']/100.0 * \
                    (1.0-self.litters_a[c]['PreTrans']))
            self.TotInpK = self.TotInpK + sum((self.litters_a[c]['Linput'])*self.litters_a[c]['Kcont']/100.0 * \
                    (1.0-self.litters_a[c]['KreTrans']))            
            self.totalLitter = self.totalLitter + self.litters_a[c]['Linput']
            self.totalNLitter =  self.totalNLitter + self.litters_a[c]['Linput']*self.litters_a[c]['Ncont']/100.0 * \
                (1.0-self.litters_a[c]['NreTrans'])
            self.totalPLitter =  self.totalPLitter + self.litters_a[c]['Linput']*self.litters_a[c]['Pcont']/100.0 * \
                (1.0-self.litters_a[c]['PreTrans'])
            self.totalKLitter =  self.totalKLitter + self.litters_a[c]['Linput']*self.litters_a[c]['Kcont']/100.0 * \
                (1.0-self.litters_a[c]['KreTrans'])
        
        for c in self.litters_b.keys():
            self.TotInp = self.TotInp + sum(self.litters_b[c]['Linput'])
            self.TotInpN = self.TotInpN + sum((self.litters_b[c]['Linput'])*self.litters_b[c]['Ncont']/100.0 * \
                (1.0-self.litters_b[c]['NreTrans']))
            self.TotInpP = self.TotInpP + sum((self.litters_b[c]['Linput'])*self.litters_b[c]['Pcont']/100.0 * \
                (1.0-self.litters_b[c]['PreTrans']))
            self.TotInpK = self.TotInpK + sum((self.litters_b[c]['Linput'])*self.litters_b[c]['Kcont']/100.0 * \
                (1.0-self.litters_b[c]['KreTrans']))
            self.totalLitter = self.totalLitter + self.litters_b[c]['Linput']
            self.totalNLitter =  self.totalNLitter + self.litters_b[c]['Linput']*self.litters_b[c]['Ncont']/100.0 * \
                (1.0-self.litters_b[c]['NreTrans'])           
            self.totalPLitter =  self.totalPLitter + self.litters_b[c]['Linput']*self.litters_b[c]['Pcont']/100.0 * \
                (1.0-self.litters_b[c]['PreTrans'])           
            self.totalKLitter =  self.totalKLitter + self.litters_b[c]['Linput']*self.litters_b[c]['Kcont']/100.0 * \
                (1.0-self.litters_b[c]['KreTrans'])           
        
        
        # ---------Initial storage for mass balance computation ------------------
        self.InitMass = (self.p['aboveMass'][rg] + self.p['belowMass'][rg])/10000.0*te  #Initial soil organic matter mass in kg m-2 
        self.InitN = (self.p['aboveN'][rg] + self.p['belowN'][rg])/10000.0*te           #Ititial nitrogen content in soil organic matter in kg m-2
        self.InitP = (self.p['aboveP'][rg] + self.p['belowP'][rg])/10000.0*te           #Ititial nitrogen content in soil organic matter in kg m-2
        self.InitK = (self.p['aboveK'][rg] + self.p['belowK'][rg])/10000.0*te           #Ititial nitrogen content in soil organic matter in kg m-2

        # ----------Initialization of the soil organic matter succession stages ------------------------- 
        Ffrac = 0.9; Hfrac = 1.0-Ffrac                                                  #Dividing initial mass into F and H pool
        self.L = np.zeros(self.length)
        self.F = np.zeros(self.length); self.F[0] = self.p['aboveMass'][rg]/10000.0 * Ffrac * te
        self.LH = np.zeros(self.length)
        self.H = np.zeros(self.length); self.H[0] = self.p['aboveMass'][rg]/10000.0 * Hfrac * te
        self.SH = np.zeros(self.length); self.SH[0] = self.p['belowMass'][rg]/10000.0 * Hfrac * te        
        
        self.NL = np.zeros(self.length)
        self.NF = np.zeros(self.length); self.NF[0] = self.p['aboveN'][rg]/10000.0 * Ffrac * te
        self.NLH = np.zeros(self.length); self.NLH[0] = self.p['belowN'][rg]/10000.0 * Ffrac * te        
        self.NH = np.zeros(self.length); self.NH[0] = self.p['aboveN'][rg]/10000.0 * Hfrac * te
        self.NSH = np.zeros(self.length); self.NSH[0] = self.p['belowN'][rg]/10000.0 * Hfrac * te    

        self.PL = np.zeros(self.length)
        self.PF = np.zeros(self.length); self.PF[0] = self.p['aboveP'][rg]/10000.0 * Ffrac * te
        self.PLH = np.zeros(self.length); self.PLH[0] = self.p['belowP'][rg]/10000.0 * Ffrac * te        
        self.PH = np.zeros(self.length); self.PH[0] = self.p['aboveP'][rg]/10000.0 * Hfrac * te
        self.PSH = np.zeros(self.length); self.PSH[0] = self.p['belowP'][rg]/10000.0 * Hfrac * te    

        self.KL = np.zeros(self.length)
        self.KF = np.zeros(self.length); self.KF[0] = self.p['aboveK'][rg]/10000.0 * Ffrac * te
        self.KLH = np.zeros(self.length); self.KLH[0] = self.p['belowK'][rg]/10000.0 * Ffrac * te        
        self.KH = np.zeros(self.length); self.KH[0] = self.p['aboveK'][rg]/10000.0 * Hfrac * te
        self.KSH = np.zeros(self.length); self.KSH[0] = self.p['belowP'][rg]/10000.0 * Hfrac * te    

        self.totSto = np.zeros(self.length);self.totNSto = np.zeros(self.length) 
        self.totPsto = np.zeros(self.length); self.totKsto = np.zeros(self.length)
        self.Pminer = np.zeros(self.length); self.Kminer=np.zeros(self.length)        
        
    def decomRomul(self):
        """
        Run Romul model, computation timestep is one day, which is here integrated over the timestep of Plantation simulator.         
        """
        def fk1(Ncont, ashCont, ligninCont, soilpH, soilT, soilW, aboveGr=True): 
            #Description of Romul, Table 1
            ligninCont = np.array(ligninCont)
            psii = 0.0005+0.0054*Ncont  if aboveGr==True else np.maximum(0.0136+0.0006*ashCont, 0.0) #ashContent range 0...20%
            mu = 1.0 if ligninCont < 10.0 else  0.092*(ligninCont/Ncont)**-0.7396            
            ita = 0.701*soilpH-1.618-0.038*soilpH**2.0
            fT = np.maximum(0.0, np.where(soilT<=35.0,np.where(soilT<=1.0, 0.1595 + 0.0319*soilT, 
                                                         0.1754*np.exp(0.0871*soilT)), 8.791-0.1465*soilT))
            fW = np.where(soilW <=1.333, np.minimum(1.0, np.where(soilW < 0.023, 0.0, 9.297*soilW**2.5493)), 29.53*0.07889**soilW)
            return psii * mu * ita * fT * fW
            
        def fk2(Ncont, ligninCont, soilpH, soilT, soilW, aboveGr = True):
            ligninCont = np.array(ligninCont)
            psii = 0.000496 if aboveGr == True else 0.00126
            mu = 1.0 if ligninCont < 10.0 else 0.0027*(ligninCont/Ncont)**-0.3917            
            ita = 0.701*soilpH-1.618-0.038*soilpH**2.0
            fT= np.maximum(0.0, np.where(soilT <=35, np.where(soilT<25.0, np.where(soilT<=1.0, 0.1595+0.0319*soilT, 0.1754*np.exp(0.0871*soilT)), 1.534),
                     3.69-0.0615*soilT))
            fW = np.where(soilW <=1.333, np.minimum(1.0, np.where(soilW < 0.023, 0.0, 9.297*soilW**2.5493)), 29.53*0.07889**soilW)
            return psii * mu * ita * fT * fW
        
        def fk3(Ncont, ashCont, ligninCont, soilT, soilW, aboveGr=True):
            ligninCont = np.array(ligninCont)
            psii = 0.0089+0.0078*Ncont if aboveGr==True else np.maximum(0.0394-0.0021*ashCont, 0.0) 
            mu = 1.0 if ligninCont <10.0 else 0.0622*(ligninCont/Ncont)**-0.397
            ita = 1.0
            fT = np.maximum(0.0, np.where(soilT < 7.0, np.where(soilT <= 0.0, 1.3/(1.0+1.97*soilT**2.0),1.3), 
                                                     (78.0-1.3*soilT)/53.0))
            fW = np.where(soilW <=1.333, np.minimum(1.0, np.where(soilW < 0.023, 0.0, 9.297*soilW**2.5493)), 29.53*0.07889**soilW)
            return psii * mu * ita * fT * fW
    
        def fk4(Ncont, soilT, soilW):
            psii = np.where(Ncont<2.0, 0.0005*Ncont, 0.001)
            mu = 1.0
            ita = 1.0
            fT = np.maximum(0.0, np.where(soilT <=40.0, np.where(soilT<=20.0, np.where(soilT<=1.0, 0.1595+0.0319*soilT, 0.1754*np.exp(0.0871*soilT)), 1.0),
                     2.0-0.025*soilT))
            fW = np.maximum(0.0, np.where(soilW <=1.333, np.minimum(1.0, 7.5*soilW),2.333-soilW))
            return psii * mu * ita * fT * fW
        
        def fk5(Ncont, soilT, soilW):
            psii = np.maximum(np.where(Ncont<=2.0, 0.007*(2.0*Ncont-1)/3.0, 0.007), 0.0)
            mu = 1.0
            ita = 1.0
            fT = np.maximum(0.0, np.where(soilT <=25.0, np.where(soilT<=13.0, np.where(soilT<=1.0, 0.078+0.0156*soilT, 0.0675*np.exp(0.2088*soilT)), 1.0),
                     2.0-0.04*soilT))
            fW = np.maximum(0.0, np.where(soilW<=2.333, np.minimum(np.where(soilW < 0.067, 0.0, 2.307*soilW-0.1538), 1.0), 2.4-0.6*soilW))
            return psii * mu * ita * fT * fW
        
        def fk6(soilT, soilW):
            psii = 0.0006     #sandy soils 0.0006 -> heavy clay 0.00006 (Raomul description)       
            fT = np.maximum(0.0, np.where(soilT <=35.0, np.where(soilT<=27.5, np.where(soilT<=1.0, 0.1595+0.0319*soilT, 0.1754*np.exp(0.0871*soilT)), 1.95),
                     4.68-0.078*soilT))
            fW = np.where(soilW <=1.333, np.minimum(1.0, np.where(soilW < 0.023, 0.0, 9.297*soilW**2.5493)), 29.53*0.07889**soilW)
            return psii * fT * fW
            
        def fGetRestrictions(Ncont, ashCont, ligninCont, soilpH, soilT, soilW, aboveGr):
            Ncont = 0.01 if Ncont < 0.01 else Ncont
            k1 = fk1(Ncont, ashCont, ligninCont, soilpH, soilT, soilW, aboveGr=aboveGr)
            k2 = fk2(Ncont, ligninCont, soilpH, soilT, soilW, aboveGr = aboveGr)       
            k3 = fk3(Ncont, ashCont, ligninCont, soilT, soilW, aboveGr=aboveGr)       
            k4 = fk4(Ncont, soilT, soilW)
            k5 = fk5(Ncont, soilT, soilW)
            k6 = fk6(soilT, soilW)       
            return k1, k2, k3, k4, k5, k6

        def fMF(nlit, lmass, nf, fmass):
            crit = nf/fmass*100.0-1.16*nlit/lmass*100.0
            MF = 0.1 if crit <0.44 else (0.5 if crit < 1.5 else 1.0)
            return MF           

        def fMSH(sh, nsh):
            msh = 0.5 if nsh==0.0 else (1.0 if sh*0.5/nsh < 8.0 else 0.5)
            return msh
            
       #---------Romul parameters & site parameters ------------------
        time = np.arange(0, self.length, step = self.dt); rdt = 30              #time months from beginning, Romul integration step in days
        delta_a = 24.0; delta_l = 12.8; gamma = 0.8                             #Parameters from the documentation  
        ML = 0.1; MH = 1.0; d = 0.5
        soilT = self.soilT; soilW = self.soilW; soilpH = self.soilpH                          
     
        #--------Local state variables --------------------------        
        L=np.zeros(self.length); cL=np.zeros(self.length); L3a=np.zeros(self.length)
        L3b=np.zeros(self.length);cL3a=np.zeros(self.length); cL3b=np.zeros(self.length)
        NL=np.zeros(self.length); cNL=np.zeros(self.length); NL3a=np.zeros(self.length)
        NL3b=np.zeros(self.length); cNL3a=np.zeros(self.length);cNL3b=np.zeros(self.length)
        miner = np.zeros(self.length); Nminer = np.zeros(self.length)

        PL=np.zeros(self.length); cPL=np.zeros(self.length); PL3a=np.zeros(self.length)
        PL3b=np.zeros(self.length); cPL3a=np.zeros(self.length);cPL3b=np.zeros(self.length)
        Pminer = np.zeros(self.length)

        KL=np.zeros(self.length); cKL=np.zeros(self.length); KL3a=np.zeros(self.length)
        KL3b=np.zeros(self.length); cKL3a=np.zeros(self.length);cKL3b=np.zeros(self.length)
        Kminer = np.zeros(self.length)

        #-----Computation of litter cohorts separately; first above then belowground cohorts -------------------
        for c in self.litters_a.keys():
            NreTrans = self.litters_a[c]['NreTrans']; Ncont = self.litters_a[c]['Ncont']*(1-NreTrans) 
            PreTrans = self.litters_a[c]['PreTrans']; Pcont = self.litters_a[c]['Pcont']*(1-PreTrans) 
            KreTrans = self.litters_a[c]['KreTrans']; Kcont = self.litters_a[c]['Kcont']*(1-KreTrans)            
            ashCont = self.litters_a[c]['ashCont']; ligninCont =self.litters_a[c]['ligninCont']; aboveGr =self.litters_a[c]['above']
            Linput =self.litters_a[c]['Linput']
            #k1,k2,k3,k4,k5,k6 = fGetRestrictions(Ncont, ashCont, ligninCont, soilpH, soilT, soilW, aboveGr=aboveGr )
            cL = cL*0.0; cNL = cNL*0.0; cL3a = cL3a*0.0; cNL3a = cNL3a*0.0            
            cPL = cPL*0.0; cPL3a = cPL3a*0.0; cKL = cKL*0.0; cKL3a = cKL3a*0.0             
            m = 0.0; nm = 0.0; pm = 0.0; km = 0.0            
            for n in range(self.length):
                k1,k2,k3,k4,k5,k6 = fGetRestrictions(Ncont, ashCont, ligninCont, soilpH, soilT[n], soilW, aboveGr=aboveGr )
                cL[n] = Linput[n] + m * (np.exp(-(k1+k3)*rdt))
                cNL[n] = Linput[n]*(Ncont/100.0) + nm*np.exp(-(ML*k1+k3)*rdt)
                cPL[n]= Linput[n]*(Pcont/100.0) + pm*np.exp(-(k1+k3)*rdt)
                cKL[n]= Linput[n]*(Kcont/100.0) + km*np.exp(-(k1+k3)*rdt)  #Tähän suuremmat kertoimet
                cL3a[n] = m *(1.0 - np.exp(-k3*rdt))
                cNL3a[n] = nm *(1.0 - np.exp(-k3*rdt))  
                cPL3a[n] = pm *(1.0 - np.exp(-k3*rdt))  
                cKL3a[n] = km *(1.0 - np.exp(-k3*rdt))  
                miner[n] = miner[n] + m*(1.0 - np.exp(-k1*rdt))
                Nminer[n] = Nminer[n] + nm*(1.0 - np.exp(-ML*k1*rdt))
                Pminer[n] = Pminer[n] + pm*(1.0 - np.exp(-k1*rdt))
                Kminer[n] = Kminer[n] + km*(1.0 - np.exp(-k1*rdt))
                m = cL[n];nm = cNL[n]; pm = cPL[n]; km = cKL[n]      
            self.litters_a[c]['Mstorage']=cL.copy(); self.litters_a[c]['Nstorage']=cNL.copy()
            self.litters_a[c]['Pstorage']=cPL.copy();self.litters_a[c]['Kstorage']=cKL.copy()
            L = L + cL; NL = NL + cNL; PL = PL + cPL; KL = KL + cKL;
            L3a = L3a + cL3a; NL3a = NL3a + cNL3a ; PL3a = PL3a + cPL3a ; KL3a = KL3a + cKL3a
        
        #-------- Below ground litter cohorts here ----------------------            
        for c in self.litters_b.keys():
            NreTrans = self.litters_b[c]['NreTrans']; Ncont = self.litters_b[c]['Ncont']*(1-NreTrans) 
            Pcont = self.litters_b[c]['Pcont']*(1-PreTrans) ;  Kcont = self.litters_b[c]['Kcont']*(1-KreTrans) 
            ashCont = self.litters_b[c]['ashCont']; ligninCont =self.litters_b[c]['ligninCont']; aboveGr =self.litters_b[c]['above']
            Linput =self.litters_b[c]['Linput']

            #k1,k2,k3,k4,k5,k6 = fGetRestrictions(Ncont, ashCont, ligninCont, soilpH, soilT, soilW, aboveGr=aboveGr )
            cL = cL*0.0; cNL = cNL*0.0; cL3b = cL3b*0.0; cNL3b = cNL3b*0.0            
            cPL = cPL*0.0; cPL3a = cPL3a*0.0; cKL = cKL*0.0; cKL3a = cKL3a*0.0             
            m = 0.0; nm = 0.0; pm = 0.0; km = 0.0            
            for n in range(self.length):
                k1,k2,k3,k4,k5,k6 = fGetRestrictions(Ncont, ashCont, ligninCont, soilpH, soilT[n], soilW, aboveGr=aboveGr )
                cL[n] = Linput[n] + m * (np.exp(-(k1+k3)*rdt))
                cNL[n] = Linput[n]*(Ncont/100.0) + nm*np.exp(-(ML*k1+k3)*rdt)
                cPL[n]= Linput[n]*(Pcont/100.0) + pm*np.exp(-(k1+k3)*rdt)
                cKL[n]= Linput[n]*(Kcont/100.0) + km*np.exp(-(k1+k3)*2*rdt)  #Tähän suuremmat kertoimet
                cL3b[n] = m *(1.0 - np.exp(-k3*rdt))
                cNL3b[n] = nm*(1.0 - np.exp(-k3*rdt))  
                cPL3a[n] = pm *(1.0 - np.exp(-k3*rdt))  
                cKL3a[n] = km *(1.0 - np.exp(-k3*rdt))  
                miner[n] = miner[n] + m*(1.0 - np.exp(-k1*rdt))
                Nminer[n] = Nminer[n] + nm*(1.0 - np.exp(-ML*k1*rdt))
                Pminer[n] = Pminer[n] + pm*(1.0 - np.exp(-k1*rdt))
                Kminer[n] = Kminer[n] + km*(1.0 - np.exp(-k1*rdt))
                m = cL[n];nm = cNL[n]; pm = cPL[n]; km = cKL[n]      
            self.litters_b[c]['Mstorage']=cL.copy(); self.litters_b[c]['Nstorage']=cNL.copy()
            self.litters_b[c]['Pstorage']=cPL.copy();self.litters_b[c]['Kstorage']=cKL.copy()
            L = L + cL; NL = NL + cNL; PL = PL + cPL; KL = KL + cKL;
            L3a = L3a + cL3a; NL3a = NL3a + cNL3a ; PL3a = PL3a + cPL3a ; KL3a = KL3a + cKL3a

        #-----Computation of F and H, all cohorts are lumped together --------------
        self.L = L; self.NL = NL; self.PL = PL; self.KL = KL
        f = self.F[0]; lh= self.LH[0]; h = self.H[0]; sh = self.SH[0]  
        nf = self.NF[0]; nlh = self.NLH[0]; nh = self.NH[0]; nsh = self.NSH[0]
        pf = self.PF[0]; plh = self.PLH[0]; ph = self.PH[0]; psh = self.PSH[0]
        kf = self.KF[0]; klh = self.KLH[0]; kh = self.KH[0]; ksh = self.KSH[0]        

        #----Computed thru all Plantation simulator time steps, in each integrated over all the days in the timestep        
        for m in range(self.length):
            #------parameters--------
            Ncont = max(self.NL[m]/self.L[m]*100.0, 0.01); ashCont=2; ligninCont = 0.0
            MF = fMF(self.NL[m], self.L[m], nf, f)     
            MSH = fMSH(sh, nsh) 
            k1,k2,k3,k4,k5,k6 = fGetRestrictions(Ncont, ashCont, ligninCont, soilpH, soilT[m], soilW, aboveGr=False )
            
            #-----storages------------------
            self.F[m] = L3a[m] + f*np.exp(-(k2+k4+k5)*rdt)
            self.LH[m] = L3b[m] + lh*np.exp(-(k2+k4+k5)*rdt)

            self.NF[m] = NL3a[m] + nf*np.exp(-(MF*k2+k4+k5) *rdt)            
            self.NLH[m] = NL3b[m] + nlh*np.exp(-(MF*k2+k4+k5) *rdt)            

            self.PF[m] = PL3a[m] + pf*np.exp(-(k2+k4+k5) *rdt)            
            self.PLH[m] = PL3b[m] + plh*np.exp(-(k2+k4+k5) *rdt)            
                   
            self.KF[m] = KL3a[m] + kf*np.exp(-(k2+k4+k5) *rdt)            
            self.KLH[m] = KL3b[m] + klh*np.exp(-(k2+k4+k5) *rdt)            

            self.H[m] = nf*(1.0 -np.exp(-k4*delta_a*d*rdt)) + \
                    h*np.exp(-(k5+k6)*rdt)
            self.NH[m] = nf*(1.0-np.exp(-k4*gamma*d*rdt)) + \
                    nh*np.exp(-(k5+k6)*rdt)
            
            # P and K storages in H are  computed using mass

            cP = 0 if self.F[m] == 0 else self.PF[m]/self.F[m] 
            cK = 0 if self.F[m] == 0 else self.KF[m]/self.F[m] 

            self.PH[m] = cP*(nf*(1.0 -np.exp(-k4*delta_a*d*rdt)) + \
                    h*np.exp(-(k5+k6)*rdt))
            self.KH[m] = cK*(nf*(1.0 -np.exp(-k4*delta_a*d*rdt)) + \
                    h*np.exp(-(k5+k6)*rdt))
            
            self.SH[m] = nf*(1.0-np.exp(-((delta_a*(1-d)*k4) + delta_l*k5)*rdt)) + \
                    nlh*(1.0-np.exp(-(delta_a*k4 + delta_l*k5)*rdt)) + \
                    nh*(1.0-np.exp(-k5*delta_l*rdt)) + \
                    sh*(np.exp(-k6*rdt)) 

            self.NSH[m] = nf*(1.0-np.exp(-gamma*(k4*(1-d)+k5)*rdt))  + \
                    nlh*(1.0-np.exp(-gamma*(k4+k5)*rdt)) +\
                    nh*(1.0-np.exp(-gamma*k5*rdt)) + \
                    nsh*(np.exp(-MSH*k6*rdt))

            cPb = 0 if self.LH[m]==0 else self.PLH[m]/self.LH[m]
            cKb = 0 if self.LH[m]==0 else self.KLH[m]/self.LH[m]                         
            self.PSH[m] = cPb*(nf*(1.0-np.exp(-((delta_a*(1-d)*k4) + delta_l*k5)*rdt)) + \
                    nlh*(1.0-np.exp(-(delta_a*k4 + delta_l*k5)*rdt)) + \
                    nh*(1.0-np.exp(-k5*delta_l*rdt)) + \
                    sh*(np.exp(-k6*rdt))) 

            self.KSH[m] = cKb*(nf*(1.0-np.exp(-((delta_a*(1-d)*k4) + delta_l*k5)*rdt)) + \
                    nlh*(1.0-np.exp(-(delta_a*k4 + delta_l*k5)*rdt)) + \
                    nh*(1.0-np.exp(-k5*delta_l*rdt)) + \
                    sh*(np.exp(-k6*rdt))) 

            
            miner[m] = miner[m] + \
                    f*(1.0-np.exp(-(k2+k4+k5)*rdt)) - nf*(1.0 - np.exp(-(delta_a*k4+delta_l*k5)*rdt))+ \
                    lh*(1.0-np.exp(-(k2+k4+k5)*rdt)) - nlh*(1.0 - np.exp(-(delta_a*k4+delta_l*k5)*rdt))+ \
                    h*(1.0-np.exp(-(k5+k6)*rdt)) - nh*(1.0-np.exp(delta_l*k5*rdt)) +\
                    sh*(1.0-np.exp(-k6*rdt))
            
            Nminer[m] = Nminer[m] + \
                    nf*(1.0-np.exp(-(MF*k2+(1-gamma)*(k4+k5))*rdt)) + \
                    nlh*(1.0-np.exp(-(MF*k2+(1-gamma)*(k4+k5))*rdt)) + \
                    nh*(1.0 - np.exp(-(k6+(1-gamma)*k5)*rdt)) +\
                    nsh*(1.0-np.exp(-MF*k6*rdt))
            
            Pminer[m] = Pminer[m] + \
                    pf*(1.0-np.exp(-(k2+k4+k5)*rdt)) - cP*nf*(1.0 - np.exp(-(delta_a*k4+delta_l*k5)*rdt))+ \
                    cPb*lh*(1.0-np.exp(-(k2+k4+k5)*rdt)) - cPb*nlh*(1.0 - np.exp(-(delta_a*k4+delta_l*k5)*rdt))+ \
                    cP*h*(1.0-np.exp(-(k5+k6)*rdt)) - cPb*nh*(1.0-np.exp(delta_l*k5*rdt)) +\
                    cPb*sh*(1.0-np.exp(-k6*rdt))
            
            Kminer[m] = Kminer[m] + \
                    kf*(1.0-np.exp(-(k2+k4+k5)*rdt)) - cK*nf*(1.0 - np.exp(-(delta_a*k4+delta_l*k5)*rdt))+ \
                    cKb*lh*(1.0-np.exp(-(k2+k4+k5)*rdt)) - cKb*nlh*(1.0 - np.exp(-(delta_a*k4+delta_l*k5)*rdt))+ \
                    cK*h*(1.0-np.exp(-(k5+k6)*rdt)) - cKb*nh*(1.0-np.exp(delta_l*k5*rdt)) +\
                    cKb*sh*(1.0-np.exp(-k6*rdt))

            f=self.F[m]; lh=self.LH[m]; nf=self.NF[m]; nlh = self.NLH[m]; h = self.H[m]; nh = self.NH[m]
            sh = self.SH[m]; nsh = self.NSH[m]; pf = self.PF[m]; plh = self.PLH[m]; ph = self.PH[m]; psh = self.PSH[m]; 
            kf = self.KF[m]; klh = self.KLH[m]; kh = self.KH[m]; ksh = self.KSH[m]
        
        totSto = self.L + self.F + self.LH + self.H + self.SH
        totSto = np.insert(totSto,0, self.InitMass)        
        deltaSto = np.diff(totSto)        
        self.mbal = self.totalLitter - miner - deltaSto
        
        totNSto = self.NL + self.NF + self.NLH + self.NH + self.NSH
        totNSto = np.insert(totNSto,0, self.InitN)        
        deltaNSto = np.diff(totNSto)        
        self.mbalN = self.totalNLitter - Nminer - deltaNSto

        self.miner = miner+self.mbal ; self.Nminer = Nminer+self.mbalN
        self.Pminer = Pminer; self.Kminer = Kminer        
        self.totSto = totSto; self.totNSto = totNSto
        
        popt = False
        if popt ==True:
            #---computing mass balance ----
            print ('---------Romul Mass balance, biomass ---------')        
            In = self.TotInp*10000.0; Out = sum(miner)*10000.0
            print ('In, kg/ha', In, ' Out, kg/ha', Out, ' In-out', In-Out)
            Sto0 = self.InitMass*10000.0
            StoEnd = (self.L[-1] + self.F[-1] + self.LH[-1] + self.H[-1] + self.SH[-1])*10000.0
            print ('Inital storage, kg/ha', Sto0, 'End storage, hg/ha', StoEnd, ' DeltaS ', StoEnd-Sto0)
            MB = In - Out - (StoEnd-Sto0)    
            print ('Mass balance, ', MB)
            
            print  ('Romul Mass balance, Nitrogen ---------------')
            In = self.TotInpN*10000.0; Out = sum(Nminer)*10000.0
            print ('In, kg/ha', In, ' Out, kg/ha', Out, ' In-out', In-Out)
            Sto0 = self.InitN*10000.0
            StoEnd = (self.NL[-1] + self.NF[-1] + self.NLH[-1] + self.NH[-1] + self.NSH[-1])*10000.0
            print ('Inital storage, kg/ha', Sto0, 'End storage, hg/ha', StoEnd, ' DeltaS ', StoEnd-Sto0)
            MB = In - Out - (StoEnd-Sto0)    
            print ('Nitrogen mass balance, ', MB)

            print  ('Romul Mass balance, Potassium ---------------')
            In = self.TotInpK*10000.0; Out = sum(Kminer)*10000.0
            print ('In, kg/ha', In, ' Out, kg/ha', Out, ' In-out', In-Out)
            Sto0 = self.InitK*10000.0
            StoEnd = (self.KL[-1] + self.KF[-1] + self.KLH[-1] + self.KH[-1] + self.KSH[-1])*10000.0
            print ('Inital storage, kg/ha', Sto0, 'End storage, hg/ha', StoEnd, ' DeltaS ', StoEnd-Sto0)
            MB = In - Out - (StoEnd-Sto0)    
            print ('Potassium mass balance, ', MB)

        if self.test==True:
            #----väliaikainen piirto        
            x = range(self.length)
            fig= plt.figure(num = 'Romul', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(25.0,12.0))   #Figsize(w,h), tuple inches 
            ax1 = fig.add_axes([0.33, 0.47, 0.3, 0.46]) #left, bottom, width, height
            ax1.set_title('Romul')    
            ax1.set_xlabel('x', fontsize = 14)
            ax1.plot(x, self.L, 'r-')
            ax1.plot(x, self.F, 'b-')
            ax1.plot(x, self.LH, 'r--')
            ax1.plot(x, self.H, 'g-')
            ax1.plot(x, self.SH, 'g--')
            ax1.plot(x, np.cumsum(self.mbal), '0.75')
            y = self.L+self.F+ self.LH+self.H + self.SH
            ax1.plot(x, y ,'k-')
            ax1.plot(x, np.cumsum(self.miner), 'co')
            ax2 = fig.add_axes([0.33, 0.05, 0.3, 0.35]) #left, bottom, width, height
            ax2.plot(x, self.NL, 'r-')
            ax2.plot(x, self.NF, 'b-')
            ax2.plot(x, self.LH, 'b--')
            ax2.plot(x, self.NH, 'g-')
            ax2.plot(x, self.NSH, 'g--')
            y = self.NL+self.NF+self.NLH+ self.NH + self.NSH
            ax2.plot(x, y,'k-')
            ax2.plot(x, np.cumsum(self.Nminer), 'co')
            ax2.plot(x, np.cumsum(self.mbalN), '0.75')

            ax3 = fig.add_axes([0.7, 0.05, 0.25, 0.35]) #left, bottom, width, height
            y = (self.NL + self.NF + self.NLH + self.NH + self.NSH)/(self.L+self.F+self.LH + self.H + self.SH)       
            ax3.plot(x, y, 'k-')
            
            plt.show()
            
class DecomPeat():
    """
    Computes CO2 flux from tropical peat as a function of depth of water table. \n
    Hoojier, A., Page, S., Canadell, J.G., Silvius, M., Kwadijk, J., Wösten, H., Jauhiainen, J. 2010. Current and future
    CO2 emissions from drained peatlands in Southeast Asia. Biogeosciences, 7: 1505-1514. \n
    Jauhiainen, J., Hoojier, A., Page, S. 2012. Carbon dioxide emissions from an Acacia plantation on peatland in Sumatra, indonesia. Biogeosciences, 9: 617-630.
    Input: \n
        -length - simulation length in number of time steps \n
        -dt - time step length in months \n
        -gwl - mean ground water table, if float, time series for water table (dwt) is computed as sine curve, \n
            if array, dwt is interpolated between observations \n
        -de - CWD instance \n
        -ro - Romul instance \n
        -optH - use Hoojier model \n
        -optSine - compose dwt time series using sine curve \n
        -optCWD - if True, use CWD carbon flux in computation of peat mass loss \n
    Output: \n
    
        
    """
    def __init__(self, length, dt, gwl=80.0, dwtAmp = 8.0, dwtArr=np.empty, peatdatafile=None):
        #------generate DWT dynamics, monthly mean values, fit sine curve
        self.time = np.arange(0, length, step = dt); self.dt = dt; self.length = length
        self.gwl = gwl
        self.dwtAmp = dwtAmp
        self.peatdatafile = peatdatafile
        if dwtArr is not np.empty:
            self.dwt=dwtArr
        
        
    def decomposePeat(self, dpara, gy, de, ro, optH = True, optSine = True, optCWD = False):
        if optSine == True:        
            fre = 0.5; amp = self.dwtAmp; pha = 0.29; offs=self.gwl             # DWT to offset, amplitude dwt fluctuation
            self.dwt = self.my_sin(self.time, fre, amp,pha, offs )        
        else:
            if self.dwt is np.empty:
                self.dwt=np.ones(self.length)*self.gwl
        #------- set parameters --------------------  
        rd = 0                                                                  # row index in decompose page        
        soilT = de.soilT                                                        # Temperature time series created in CWD instance
        pdt = self.dt * 30.0 / 365.0                                            # time step length here, self.dt in months 
        self.peatRhoInit = dpara['peatRhoInit'][rd]                             # initial bulk density, kg m-3
        self.peatRhoFinal = dpara['peatRhoFinal'][rd]                           # peat bulk density in the end of computation, kg m-3 
        peatDepth = 1.                                                         # m 
        CO2toC = 12.0 / 44.0                                                    # conversion factor from CO2 to C
        CtoMass = 2.0                                                           # conversion factor from C to mass
        refT = 28.0                                                             # reference temperature for efflux calculation, Q10 used to scale it further
        q10 = 2.0; a_q10 = q10**((soilT-refT)/10.0)                             # temperature correction to emission (Jauhiainen et al. 2012)
        fRho = interp1d([0, self.length], [self.peatRhoInit, self.peatRhoFinal])        
        self.peatRho = fRho(range(self.length))
        
        #********* Shrinking and swelling of peat *************
        #Camporese et al. 2006, Hydrological modelling in swelling/shrinking peat soils. WRR 42
        shrink = np.zeros(self.length)        
        optShrink = True
        if optShrink==True:        
            from hydro_utils import getParams, wrc
            #paramFile=r"C:/Users/alauren/OneDrive - University of Eastern Finland/codes/plantation_simulator/peat_data.xlsx"        
            #paramFile = self.peatdatafile
            peatset=['mesic', 'fibric', 'TPmean', 'TPQ1', 'TPQ3', 'woody']        
            soilpara = getParams(self.peatdatafile, shname = peatset[2], lrs=176)                                         # soil layer hydraulic parameters 
            pF=soilpara['pF']        
            nlyrs = 30        
            dz =np.array(list(pF['dz'].values()))
            z = np.array(list(pF['z'].values()))
            poro = np.array(list(pF['ThetaS'].values()))
            poro= poro[:nlyrs]
            vw = 1./1000.                                                       #specific volume of water m3/kg
            vs = 1./1600.                                                       #specific volume of solids m3/kg
            void_r = np.zeros(self.length)                                      #void ratio (vol pores / vol solids, or poro/(1-poro))           
            #moist_r = np.zeros(self.length)                                     #moisture ratio (vol water / volsolids) Vw/Vs            
            alfa = 0.2
            delta =1.0/(2*alfa+1)
            limit_v0 = 8.8                                                      #limit moisture ratio
            
            for gwl, rh, d in zip(self.dwt[:], self.peatRho, range(self.length)):
                #drymass = rh*dz[:nlyrs]
                theta_ini =wrc(pF, x = np.minimum(z-0/100.0, 0.0))
                theta_ini = theta_ini[:nlyrs]
                #watmass_ini = theta_ini*dz[:nlyrs]*1000.0
                theta = wrc(pF, x = np.minimum(z-gwl/100.0, 0.0))    
                theta = theta[:nlyrs]
                
                #watmass = theta*dz[:nlyrs]*1000.0
                #grav_w=watmass/drymass            
                #grav_w_ini = watmass_ini/drymass
                #lchange = sum(dz[:nlyrs]) - sum(((grav_w*vw + vs)/(grav_w_ini*vw + vs))**0.33*dz[:nlyrs])
                #shrink[d]= lchange*1.0              
                
                moist_r_ini = theta_ini/(1.-poro)
                ee_ini = np.where(moist_r_ini<limit_v0, (limit_v0 + 1.)**(1.-delta)*(moist_r_ini + 1.)**delta -1., moist_r_ini)
                moist_r = theta / ((1.-poro))
                ee = np.where(moist_r<limit_v0, (limit_v0 + 1.)**(1.-delta)*(moist_r + 1.)**delta -1., moist_r)
                deltah = sum(dz[:nlyrs]) - sum(((ee/(1.0+ee)) / (ee_ini/(1.0 + ee_ini))) *dz[:nlyrs])                
                shrink[d] = deltah
            #print shrink
        optStandMassRho = False                                                 # if true, compute top peat rho from stand mass, else constant increase of rho
        if optStandMassRho==True:
            """
            Duraisamy et al. 2007. Engineering Properties and Compressibility Behavior of Tropical Peat Soil. 
            American Journal of Applied Sciences 4 (10): 768-773, 2007
            Yamamoto et al. 2003.Moisture Distribution in Stems of Acacia mangium, A. auriculiformis and Hybrid Acacia Trees
            """
            excavatorGroundPressure = 90.0                                     # kPa, 15 ton excavator with 40...50 kPa, heavy 90 kPa
            cutsIn = [int((gy.g['Lrotat'][0]*12 / gy.g['dt'][0] ) *r) for r in range(int(gy.g['Nrotat'][0]))]
            cutsIn = cutsIn[-1]            
            standMass = gy.BiBark + gy.BiBranch + gy.BiFoliage + gy.BiStem + gy.RootMass +gy.weedAbove+ gy.weedBelow    #kg ha-1
            litterMass = (ro.L+ro.F+ro.H) * 10000.0                             # kg ha-1
            gms = 140.0                                                         # gravimetric moisture of stems %, Yamamoto et al             
            gml = 500.0                                                         # gravimetric moisure of litter %, guess           
            mass = standMass*(gms+100.0)/100.0 + litterMass*(gml+100.)/100.0 + de.cwdSto*(gml+100.)/100.0    # fresh mass 
            g = 9.81                                                            # acceleration of gravity m s-2
            Fd = mass /10000.0*g/1000.0                                         # downward force kPa (mass kg/ha->kg/m2, to N, to N/m2 (Pa), and to kPa)                    
            Fd[cutsIn]=Fd[cutsIn] + excavatorGroundPressure            
            logdeltaef = np.log10(Fd) - np.log10(Fd[0])                         # change in the downward force
            fii0 = 1.0 - self.peatRho/1000./1.6                                 # porosity [-] without compression
            e0 = fii0/(1.0-fii0)                                                # void ratio without compression
            cc={'fibric':2.752, 'hemic':2.165, 'sapric':1.935}                  # compressibility indices for tropical peats from Duraisamy et al.
            deltae = logdeltaef * cc['hemic']                                   # change in void ratio
            e1 = e0 - deltae                                                    # new void ratio
            fii1 = e1/(e1+1.)                                                   # new porosity
            rho1 = 1.6*(1.0-fii1)                                               # new bulk density
            
            #*************Linear delay difference equation *********************(see Kolari et al 2007. Temperature... Tellus 59B, 3)
            tau = 4.0                                                           # time delay, in global dt (months)
            tauS = 2.0                                                           # time delay for shrinking and swelling
            dSdt = np.zeros(len(rho1))                                          # delayed change in rho1
            Srho = np.zeros(len(rho1)); Srho[0]=rho1[0]                         # delayed rho1
            length = len(rho1)
            dShrink = np.zeros(len(rho1));Sshrink = np.zeros(len(rho1)); Sshrink[0]=shrink[0] 
            for d in range(1,length):
                dSdt[d] = (rho1[d]-Srho[d-1])/tau
                Srho[d]=dSdt[d]+Srho[d-1]
                dShrink[d] = (shrink[d]-Sshrink[d-1])/tauS
                Sshrink[d] = dShrink[d] + Sshrink[d-1]
            shrink= Sshrink
            """
            import matplotlib.pylab as plt
            fig=plt.subplot(211); plt.plot(range(length), rho1, 'b')
            plt.plot(range(length), Srho, 'g-');plt.show()
            plt.subplot(212); plt.plot(range(length), Fd, 'k--')            
            """            
            #for r, s in zip(rho1,Srho):
            #    print r,s
            self.peatRho = Srho*1000.0    

        #initPeatStorage = 10000.0 * peatDepth * self.peatRhoInit                # kg/ha
        initPeatStorage = 10000.0 * peatDepth * self.peatRho[0]                # kg/ha    
        # ------Efflux models ---------------
        if optH == True:        
            # Hoojier 2010:
            CO2efflux = (91.0*self.dwt/100.0*1000.0*pdt)*a_q10                  # kg CO2 /ha/timestep 
            print ('Hoojier 2010')
        else:
            # Jauhiainen 2012:
            #self.dwt=np.append(self.dwt,np.mean(self.dwt))   # TEMP!!!!!! REMOVE THIS
            CO2efflux = ((71.1*self.dwt/100.0 + 23.15)*1000.0*pdt)*a_q10
            print ('Jauhiainen 2012')

        #------Divide mass loss between peat, cwd and stand litter --------
        cwdMloss = de.cwdMloss*-1.0                                             # kg / ha / timestep
        romulMloss = ro.miner*10000.0                                           # kg / ha / timestep
        cwdShareCaught = 1.0 #0.5        
        if optCWD == True:            
            stoChange = CO2efflux*CO2toC*CtoMass  - romulMloss - cwdMloss * cwdShareCaught       # kg / ha / timestep
        else:
            #-----Assumption: outflux from CWD is not caught with chambers         
            stoChange = CO2efflux*CO2toC*CtoMass - romulMloss                   # kg / ha / timestep        
        stoChange[stoChange<0]=0.0
        peatStorage = initPeatStorage - np.cumsum(stoChange)
        
        # -------compute peat subsidence - separate oxidation and consolidation 
        soilSto = peatStorage + de.cwdSto + (ro.SH + ro.H + ro.LH + ro.F + ro.L)*10000.0
        surfAfterOxidation = (soilSto-soilSto[0])/10000.0/self.peatRhoInit*100
        
        #--------combined peat oxidation and consolidation----------
        subsidence = soilSto / 10000.0/self.peatRho*100.0
        subsidence = subsidence - subsidence[0] - shrink*100.
        #print subsidence                        
        annualSubsidence = np.diff(subsidence)/pdt
        
        #--------Share of oxidation from the total subsidence (75-92% Hoojier et al. 2012, Biogeosciences 9: 1053-1071)         
        oxidShare = (surfAfterOxidation - surfAfterOxidation[0])/(subsidence)
        oxidShare = np.zeros(len(subsidence))        
        stmp = surfAfterOxidation - surfAfterOxidation[0]
        oxidShare = np.where(subsidence < 0.0, stmp/subsidence, 0.0)

        #--------Locate to instance variables------------------
        self.CO2efflux = CO2efflux; self.peatCrelease = stoChange / CtoMass 
        self.peatStorage = peatStorage; self.surfAfterOxidation = surfAfterOxidation
        self.subsidence = subsidence; self.oxidShare = oxidShare        
        self.annualSubsidence=annualSubsidence
        self.CtoMass = CtoMass 
    
    def fitDWT(self):   
        #Curve fit here
        optFit = False
        if optFit == True:        
            dwt = [50.0, 60.0, 50.0, 40.0, 45.0, 50.0 ]; dwt = np.array(dwt)
            months = [0, 2, 4, 6, 8, 10]; months= np.array(months)
            guess_freq = 1/12.0
            guess_amplitude = 3*np.std(dwt)/(2**0.5)
            guess_phase = 12/(2.0*np.pi)
            guess_offset = np.mean(dwt)
            
            p0=[guess_freq, guess_amplitude,
                guess_phase, guess_offset]
            
            # now do the fit
            #my_sin = lambda time, *p0: np.sin(time * p0[0] + p0[2]) * p0[1] + p0[3]
                  
            fit = curve_fit(self.my_sin, months, dwt, p0=p0)   
            
            fre= fit[0][0]
            amp = fit[0][1]
            pha = fit[0][2]
            offs = fit[0][3] 
            print (fre, amp, pha, offs)        
        
        
        fig= plt.figure(num = 'Peat', facecolor=(232/255.0, 243/255.0, 245.0/255), edgecolor='k',figsize=(25.0,12.0))   #Figsize(w,h), tuple inches 
        ax1 = fig.add_axes([0.33, 0.47, 0.3, 0.46]) #left, bottom, width, height
        ax1.set_title('Peat')    
        ax1.set_xlabel('x', fontsize = 14)
        #ax1.scatter(months, dwt, color ='r', marker= 'o')
        ax1.plot(self.time, self.dwt, 'k-')
        plt.show()

        return fre, amp, pha,offs

    def my_sin(self, x, freq, amplitude, phase, offset):
        # create the function we want to fit
        return np.sin(x * freq + phase) * amplitude + offset
