# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 17:08:23 2016

@author: Ari Lauren
"""
import numpy as np
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline as interS
from scipy.interpolate import interp1d
from scipy.ndimage.interpolation import shift

class cGrowthAndYield():
    def __init__(self, para, gen, bi, thin):
        #----parameters---------        
        self.p = para  #Growth and yield model parameters  
        self.g = gen   #General parameters
        self.b = bi    #Biomass model parameters
        self.thin = thin #THinning parameters

        #-----initialize instance variables ----------------        
        self.nsteps = int()
        self.age = np.empty(0)                                                  # biological age of the stand, yrs
        self.time = np.empty(0)                                                 # time, yrs, from the beginning of the simulation until the end (including multiple rotations)
        self.appAge = np.empty(0)                                               # apparent age, yrs modified by Type 1 growth response
        self.appSI = np.empty(0)                                                # apparanet site index, m modified by Type 2 growth response
        self.Hdom = np.empty(0)                                                 # dominant height, m
        self.Survival = np.empty(0)                                             # number of stems, N ha-1
        self.BA = np.empty(0)                                                   # basal area m2 ha-1
        self.WMeanDiam = np.empty(0)                                            # basal area weighed mean diameter, cm
        self.PulpY = np.empty(0)                                                # pulp yield, tn ha-1
        self.WoodConsumption = np.empty(0)                                      # wood consumtion to ton of pulp, m3 ton-1

        self.DbhDist = np.empty(0)                                              # dbh distribution 2D: age, diam class
        self.Vdist = np.empty(0)                                                # volume distribution 2D: age, diam class
        self.Hdist = np.empty(0)                                                # height distribution 2D: age, diam class
        self.MercDist = np.empty(0)                                             # merchantable volume distribution 2D: age, diam class
        self.BAdist = np.empty(0)                                               # basal area distribution 2D: age, diam class
        self.rotationArr = np.empty(0)                                          # number of current rotation
        self.V =np.empty(0)                                                     # total volume of the stand 1D: age
        self.MV =np.empty(0)                                                    #Merchantable volume of the stand 1D: age
        
        self.BiBark = np.empty(0)                                               # bark biomass, kg ha-1
        self.BiBranch = np.empty(0)                                             # branch biomass, kg ha-1
        self.BiFoliage = np.empty(0)                                            # foliage biomass, kg ha-1
        self.BiStem = np.empty(0)                                               # stem biomass, kg ha-1
        self.RootMass = np.empty(0)                                             # root biomass, kg ha-1
        self.MerchMass = np.empty(0)                                            # merchantable biomass, kg ha-1
        self.LAI = np.empty(0)                                                  # leaf area index, m2 m-2

        self.weedLAI = np.empty(0)                                              # weed leaf area, m2 m-2
        self.weedAbove = np.empty(0)                                            # above ground weed mass, kg ha-1
        self.weedBelow = np.empty(0)                                            # below ground weed mass, kg ha-1                     
        self.weedings = np.empty(0)
        
        self.FineRLitL = np.empty(0)                                            # fineroot litter from living trees, kg ha-1 timestep-1
        self.RootLitD = np.empty(0)                                             # root litter from dead trees, kg ha-1 timestep-1
        self.BrLitD = np.empty(0)                                               # branch litter from dead trees, kg ha-1 timestep-1
        self.BrLitL = np.empty(0)                                               # branch litter from living trees, kg ha-1 timestep-1
        self.BaLitD = np.empty(0)                                               # bark litter from dead trees, kg ha-1 timestep-1
        self.BaLitL = np.empty(0)                                               # bark litter from living trees, kg ha-1 timestep-1
        self.FoLitD = np.empty(0)                                               # foliage litter from dead trees, kg ha-1 timestep-1
        self.FoLitL = np.empty(0)                                               # foliage litter from living trees, kg ha-1 timestep-1
        self.CWD = np.empty(0)                                                  # coarse woody debris from dead trees, kg ha-1 timestep-1
        self.weedALitL = np.empty(0)                                            # above ground weed litter from living weeds, kg ha-1 timestep-1         
        self.weedBLitL = np.empty(0)                                            # below ground weed litter from living weeds, kg ha-1 timestep-1
        self.weedALitD = np.empty(0)                                            # above ground weed litter from dead weeds, kg ha-1 timestep-1         
        self.weedBLitD = np.empty(0)                                            # below ground weed litter from dead weeds, kg ha-1 timestep-1         
        
        self.currentRowGe =  np.empty(0)                                        # array of current row codes in General parameters   
        self.currentRowGr = np.empty(0)                                         # array of current row codes in Growth and yield parameters
        
        self.calibrateBA = 1.0                                                  #adjusts tha basal area according to measurements
        
        self.FoGrTree = np.empty(0)                                             # foliage growth per tree, kg tree-1 timestep-1 
        self.BaGrTree = np.empty(0)                                             # bark growth per tree, kg tree-1 timestep-1 
        self.BrGrTree = np.empty(0)                                             # branch growth per tree, kg tree-1 timestep-1 
        self.RootGrTree = np.empty(0)                                           # root growth per tree, kg tree-1 timestep-1 
        self.FineRootGrTree = np.empty(0)                                       # fine root growth per tree, kg tree-1 timestep-1 
        self.StemGrTree = np.empty(0)                                           # stem growth per tree
        
    def fComputeGrowthAndYield(self, dwtArr=None, surTs =None, rowGr =[0,0,0,0,0,0], rowGe=[0,0,0,0,0,0]):
        Nrotat = int(self.g['Nrotat'][0]) 
        Nmts = int(self.g['Lrotat'][0]*12+1)        
        self.rowGr =rowGr; self.rowGe =rowGe                    # these set which row to input in consequtive rotations
        for nr, r, rg in zip(range(Nrotat), self.rowGr, self.rowGe):
            if dwtArr is not None:
                tmpa = dwtArr[nr*Nmts:(nr+1)*Nmts]    
            else:
                tmpa=None
            self.rotation = nr            
            self.fHdom(rowGrowth = r, rowGen = rg, dwtArr=tmpa, surTs=surTs)
            self.fDbhDistribution(rowGrowth = r, rowGen = rg)
            self.fBiomass(rowGen = rg)
            self.fLitter(rowGen = rg)
            self.fWeeds(rowGrowth = r, rowGen = rg)
        
        
    def fHdom(self, rowGrowth=0, rowGen = 0, dwtArr=None, surTs =None):
        """
        Computes Stand variables and appends them into instance variables: 
            Dominant height (Hdom, m) \n
            Basal area (BA, m2) \n
            Survival (Number of stems/ha) \n
            BA-weighed mean breast height diameter \n
            Consumption of wood for air-dry ton \n
        Uses models from \n 
             (Forss, E., Gadow, K.v., Saborowski, J. 1996. Growth models for unthinned Acacia mangium plantations in South Kalimantan, Indonesia. Journal of Tropical Forest Research, 8: 449-462)
             Growth response variables (Snowdon, P. 2002. Modeling Type 1 and Type 2 growth responses in plantations after application of fertilizer or other silvicultural treatments. Forest Ecology and Management 163:229-244.) \n
        Input, kwargs: \n
            rowGrowth - row index in the parameter file Growth and yield paras \n
            rowGen - row index in parameter file General parameters
            dwtArr - array of simulated/measured ground water depths, in cm, positive downwards, used in adjusting mortality rate
        """
        print ('*********************')
        print ('GROWTH AND YIELD')
        print ('    + Stand')
        #print dwtArr
        #-------Setting parameters -----------------------
        r = rowGrowth; rg = rowGen; dt = self.g['dt'][0]                       # row indices, time step in months
        Lrotat = self.g['Lrotat'][0]                                           # length of rotation, in yrs
        NMonths = 12 * Lrotat                                                  # months in simulation
        SI = self.g['SI'][rg]; Iage = self.g['Iage'][0]                        # site index, m and index age, yrs
        b1 = self.p['Eq2Beta1'][r]; b2 =self.p['Eq2Beta2'][r]
        age = np.arange(0, NMonths + dt, step = dt) / 12.0 +0.001              # biological age, yrs, add 0.001 to avoid zero diviasion
        rotationArr = np.ones(len(age))*self.rotation
        
        #---------Growth response parameters, Type1 first then Type2-----------                
        tit = self.p.keys()
        apara = sorted([x for x in tit if 'T1a' in x])
        upara = sorted([x for x in tit if 'T1u' in x])
        Wpara = sorted([x for x in tit if 'T1W' in x])
        appAge = np.zeros(len(age))
        for n in range(len(apara)):
            a =self.p[apara[n]][r]; u =self.p[upara[n]][r]; W =self.p[Wpara[n]][r]                
            if a != 0:            
                grResponse = interp1d([0.0,a,a+W,age[-1]],[0.0,0.0,u,u], fill_value='extrapolate')
                appAge = appAge + grResponse(age)
        appAge = appAge + age
        del apara, upara, Wpara
        
        apara = sorted([x for x in tit if 'T2a' in x])
        upara = sorted([x for x in tit if 'T2u' in x])
        Wpara = sorted([x for x in tit if 'T2W' in x])
        appSI = np.zeros(len(age))
        for n in range(len(apara)):
            a =self.p[apara[n]][r]; u =self.p[upara[n]][r]; W =self.p[Wpara[n]][r]                
            if a != 0:            
                grResponse = interp1d([0.0,a,a+W,age[-1]],[0.0,0.0,u,u], fill_value='extrapolate')
                appSI = appSI + grResponse(age)     
        appSI = appSI + SI
        
        #---------Eq2, dominant height, m---------
        Hdom = appSI * ((1.0 - np.exp(-1.0 * b1 * appAge)) / (1.0 - np.exp(-1.0 * b1 * Iage))) ** b2

        #---------Eq4, survival, unit number of stems/ha        
        APla = 1.0 / 12.0                                                           #Age in planting 
        NPla = self.g['N'][rg]; b0 = self.p['Eq4Beta0'][r]; b1 = self.p['Eq4Beta1'][r] 
        b2 = self.p['Eq4Beta2'][r]; b3 = self.p['Eq4Beta3'][r]
        
        optThin = False
        if optThin:
            Mrate = self.g['Mrate'][rg]
            thin = self.thin
            rot = age[-1]*12./dt
            thin = thin[thin.index < rot]
            removal = np.ones(len(age))*Mrate
            removal[thin.index] = removal[thin.index] + thin['removal']
            Survival = NPla - np.cumsum(removal)
        
        else:
            indSurRate = True
            if surTs is None:
                Mrate = self.g['Mrate'][rg]
                if indSurRate==False:
                    Survival = (NPla ** b0 + (b1 + b2 / SI) * ((age / 10.0) ** b3 - (APla / 10.0) ** b3)) ** (1.0 / b0)
                else:
                    Survival = NPla - np.cumsum(np.ones(len(age))*Mrate)        
                    if dwtArr is not None:
                        dwtArr = dwtArr*-100
                        a=15. ; b=-0.3                                                # adjust the mortality rate according to water high table      
                        dwtA = np.where(dwtArr<50., a+b*dwtArr,0.0)
                        adj = np.cumsum(dwtA) 
                        Survival = np.maximum(Survival - adj, np.zeros(len(Survival)))
            else:
                Survival = surTs['Nstems'].values

        
        #--------Eq5, basal area, m2/ha
        b0 = self.p['Eq5Beta0'][r]; b1=self.p['Eq5Beta1'][r]; b2=self.p['Eq5Beta2'][r]; b3=self.p['Eq5Beta3'][r]
        BA = np.exp(b0 + b1 * 1.0 / (appAge) + b2 * np.log(Hdom) + b3 * np.log(Survival))*self.calibrateBA
        if optThin:
            BAremoval = np.zeros(len(age))             
            BAremoval[thin.index] = BAremoval[thin.index] + thin['BAremoval']
            BA = BA - np.cumsum(BAremoval)
            
        #-------Mean diameter, cm        
        WMeanDiam = (BA * 10000.0 / (3.14 * Survival)) ** 0.5 * 2.0   
        
        #--------WoodConsumption, m3/ton
        b0 = self.p['WoodBeta0'][r]; b1 = self.p['WoodBeta1'][r]
        WoodConsumption = b0 * np.log(age * 12.0) + b1 
        
        #------append results into instance variables--------
        self.age=np.append(self.age, age);
        self.appAge = np.append(self.appAge, appAge) 
        self.appSI = np.append(self.appSI, appSI); self.Hdom = np.append(self.Hdom, Hdom)
        self.Survival = np.append(self.Survival, Survival); self.BA = np.append(self.BA, BA)
        self.WMeanDiam = np.append(self.WMeanDiam, WMeanDiam); self.WoodConsumption = np.append(self.WoodConsumption, WoodConsumption)
        self.rotationArr = np.append(self.rotationArr, rotationArr);  self.time = np.append(self.time, age+(Lrotat*self.rotation))        
        self.currentRowGe = np.append(self.currentRowGe,rowGen)
        self.currentRowGr = np.append(self.currentRowGr,rowGrowth)        
        del age, appAge, appSI, Hdom, Survival, BA, WMeanDiam, WoodConsumption, rotationArr
        
    def fDbhDistribution(self, rowGrowth =0, rowGen = 0):
        """
        Computes diameter distribution:  \n
            Kwargs: \n 
            - rowGrowth - row index in parameter file growth and yield page \n
            - rowGen - row index in parameter file general parameters \n            
            Uses models from: \n
            Individual tree characteristics \n
            - describing diameter (dbh, in cm) distribution in time (Weibull distribution) \n
            - describing individual tree height (h, in m) as a function of dbh \n
            - describing individual tree volume (v, in m3) as a function of dbh and h \n
            - describing individual tree merchantable (mv, in m3/ha) volume as a function of dbh and v \n
            - describing pulp yield (m3/air dry ton) \n
            Returns variables in 2-dimensional arrays [age][dbh]: DbhDist, Hdist, Vdist, MercDist. \n
            Returns PulpYield in 1-dimensional array [age]. \n
         (Forss, E., Maltamo, M. & Saramäki, J. 1998. Static stand and tree characteristics model for Acacia mangium plantations. Journal of Tropical Forest Science 10, 318-336 ) \n

        """
        #--------Reshape variables -----------
        self.nsteps = int(len(self.age)/(self.rotation +1))
        nsteps = self.nsteps 
        appAge = self.appAge[-nsteps:]; WMeanDiam = self.WMeanDiam[-nsteps:]; Hdom = self.Hdom[-nsteps:]; Survival = self.Survival[-nsteps:]
        
        #--------Parameters---------------
        r=rowGrowth; rg = rowGen          
        b00 = self.p['Eq4Beta00'][r];b01 = self.p['Eq4Beta01'][r];b02 = self.p['Eq4Beta02'][r]
        b10 = self.p['Eq7Beta10'][r]; b11 = self.p['Eq7Beta11'][r]                 

        WeibullB = b00 + b01 * WMeanDiam + b02 * appAge   
        
        WeibullC = np.exp(b10 + b11 * appAge)  

        #------ Individual tree height  parameters Eqs 3, 8, 9
        b00 = self.p['Eq8Beta00'][r]; b01 = self.p['Eq8Beta01'][r]; b02 = self.p['Eq8Beta02'][r]
        b10 = self.p['Eq9Beta10'][r];b11 = self.p['Eq9Beta11'][r];b12 = self.p['Eq9Beta12'][r];b13 = self.p['Eq9Beta13'][r]
        Eq3Beta0 = b00 + b01 * WMeanDiam + b02 * Hdom
        Eq3Beta1 = b10 + b11 * WMeanDiam + b12 * Hdom + b13 * np.log(Survival)

        #--------Tree volume parameters, eq10, Merchantable volume, eq11
        b0 = self.p['Eq10Beta0'][r]; b1 = self.p['Eq10Beta1'][r]; b2 = self.p['Eq10Beta2'][r]
        Mlimit = self.g['Tdia'][rg]; bM0 = self.p['Eq11Beta0'][r]; bM1 = self.p['Eq11Beta1'][r]               

        #-----------Generate height, diameter, volume and merch volume distributions
        ul = int(self.g['DiamUpperLimit'][0])
        DbhDist = np.reshape(np.zeros(ul*nsteps), (nsteps, ul)); Hdist = np.reshape(np.zeros(ul*nsteps), (nsteps, ul))
        Vdist = np.reshape(np.zeros(ul*nsteps), (nsteps, ul)); MercDist = np.reshape(np.zeros(ul*nsteps), (nsteps, ul))
        WoodConsumption = self.WoodConsumption[-nsteps:]
        
        """
        for a in range(nsteps):
            for d in range(ul):
                if Hdom[a] > 1.3:
                    fld = float(d) if d>0.01 else 0.01                    
                    DbhDist[a][d]= (WeibullC[a] / WeibullB[a] * (fld / WeibullB[a]) ** (WeibullC[a] - 1.0) \
                        * np.exp(-1.0 * (fld / WeibullB[a]) ** WeibullC[a]))*Survival[a]        
                    Hdist[a][d] = 1.3 + Eq3Beta0[a] * np.exp(-1.0 * (Eq3Beta1[a] / fld)) 
                    Vdist[a][d] = (np.exp(b0 + b1 * np.log(fld) + b2 * np.log(Hdist[a][d])))* DbhDist[a][d]
                    MercDist[a][d] = Vdist[a][d] * np.exp(bM0 * (Mlimit / fld) ** bM1) 
            #if sum(DbhDist[a][:]) > Survival[a]:
            #    DbhDist[a][:] = Survival[a] /sum(DbhDist[a][:]) * DbhDist[a][:]
            print sum(DbhDist[a][:]), Survival[a], WeibullB[a], WeibullC[a]                
        #for a in range(nsteps):    
        #    print WMeanDiam[a]
        """
        seedling_limit = 1.0
        a_ind=np.where(WMeanDiam>seedling_limit)[0][0]    #index where diameter exceeds 1.0 cm
        for a in range(nsteps):
            for d in range(ul):
                if WMeanDiam[a] > seedling_limit:
                    fld = float(d) if d>0.01 else 0.01                                        
                    DbhDist[a][d]= (WeibullC[a] / WeibullB[a] * (fld / WeibullB[a]) ** (WeibullC[a] - 1.0) \
                        * np.exp(-1.0 * (fld / WeibullB[a]) ** WeibullC[a]))*Survival[a]        
                    Hdist[a][d] = 1.3 + Eq3Beta0[a] * np.exp(-1.0 * (Eq3Beta1[a] / fld)) 
                    Vdist[a][d] = (np.exp(b0 + b1 * np.log(fld) + b2 * np.log(Hdist[a][d])))* DbhDist[a][d]
                    MercDist[a][d] = Vdist[a][d] * np.exp(bM0 * (Mlimit / fld) ** bM1) 
                
            if WMeanDiam[a] < seedling_limit: DbhDist[a][0]=float(a)/a_ind * Survival[a]
            if sum(DbhDist[a][:]) > Survival[a]:
                DbhDist[a][:] = Survival[a] /sum(DbhDist[a][:]) * DbhDist[a][:]
                

        #compute total standvolume and merchantable volume
        V=np.sum(Vdist, axis=1)
        MV=np.sum(MercDist, axis=1)                                

        #------------Basal area (m2/ha) computed from dbh distribution
        BAdist = np.zeros(nsteps); PulpY = np.zeros(nsteps)
        for a in range(nsteps):
            BAdist[a] = sum((DbhDist[a][0:] * 3.14 * (np.arange(ul)/200.0 ) ** 2.0))
        PulpY = np.sum(MercDist, axis=1) / WoodConsumption
        
        #------Append results to instance variables --------------------
        self.DbhDist  = np.append(self.DbhDist, DbhDist);   self.BAdist = np.append(self.BAdist, BAdist)
        self.Vdist    = np. append(self.Vdist, Vdist);      self.Hdist  = np.append(self.Hdist, Hdist)
        self.MercDist = np.append(self.MercDist, MercDist); self.PulpY  = np.append(self.PulpY, PulpY)
        self.V = np.append(self.V, V); self.MV=np.append(self.MV, MV)
        del DbhDist, BAdist, Vdist, Hdist, MercDist, PulpY, WoodConsumption
        
    def fBiomass(self, rowGen = 0):
        """
        Computes biomass of fractions of individual trees (kg) as a function of breast height diameter: \n
            - stem \n
            - bark \n
            - branch\n
            - foliage \n
        Returns forest biomass (kg/ha) by multiplying the individual tree biomass by dbh distribution. \n 
            - 1 dimensional arrays [age]: BiStem, BiBark, BiBranch, BiFoliage. \n
            - Merchantable mass as a product of merchantable volume and wood density, 1 d array [age] \n
            - Leaf area index as (m2/m2) as a product of BiFoliage and specific leaf area. \n
        Key word arguments: \n
            - rowGen row index in parameter file General parameters
        Uses: \n
        Bi,H., Turner, J. & Lambert, M. 2004. Additive biomass functions for native eucalypt forest trees of temperate Australia.	
        Trees 18: 467-479.	 \n
        Input parameter table includes:        
        Estimated coefficients of the system equations (Eq. 6) and fit indices (R2) for the 15 tree species where D is the only independent
        variable in the equation.
        """
        print ('    + Biomass')
        #--------Parameters----------------------
        rg = rowGen; s = self.g['BiomSpe'][rg]   #species, refers to row in parameter table
        self.species = str(self.b['species'][s])[0:7]        
        b10 = self.b['b10'][s]; b11 = self.b['b11'][s]   #stem
        b20 = self.b['b20'][s]; b21 = self.b['b21'][s]   #bark
        b30 = self.b['b30'][s]; b31 = self.b['b31'][s]   #branch
        b40 = self.b['b40'][s]; b41 = self.b['b41'][s]   #foliage        
        
        rot = self.rotation; ul = int(self.g['DiamUpperLimit'][0]); nsteps = int(len(self.age)/(rot+1)) 
        rho = self.g['WoodDens'][rg]; BgAgRatio = self.g['BgAgRatio'][rg]; SLAI = self.g['SLAI'][rg]       #rho kg/m3; Belowground / Aboveground biomass [-]; Specific leaf area [m2/kg]
        BiStem = np.zeros(nsteps);BiBark = np.zeros(nsteps);BiBranch = np.zeros(nsteps)
        BiFoliage = np.zeros(nsteps); MerchMass = np.zeros(nsteps); RootMass = np.zeros(nsteps)
        LAI = np.zeros(nsteps)

        rg = rowGen
        BarkLon = self.g['BarkLon'][rg]; BranchLon = self.g['BranchLon'][rg]; LeafLon = self.g['LeafLon'][rg]
        FineRLon = self.g['FineRLon'][rg]; FineRootPr = self.g['FineRootPr'][rg]; #nsteps = len(self.age)/(self.rotation +1)
        dt = self.g['dt'][0]/12.0; Survival = self.Survival[-nsteps:]
        FoGrTree=np.zeros(nsteps); BaGrTree=np.zeros(nsteps); BrGrTree=np.zeros(nsteps);                                             # branch growth per tree, kg tree-1 timestep-1 
        RootGrTree=np.zeros(nsteps); FineRootGrTree=np.zeros(nsteps); StemGrTree=np.zeros(nsteps)            
        age = self.age[-nsteps:]
        
        #-------Time series of above ground biomass development, kg/ha
        DbhDist = self.DbhDist[nsteps*ul*rot:]
        DbhDist = np.reshape(DbhDist, (nsteps, ul))
        MercDist = np.reshape(self.MercDist[nsteps*ul*rot:], (nsteps, ul))        
        for a in range(nsteps):
            BiStem[a] = sum([np.exp(b10)*float(d)**b11*DbhDist[a][d] for d in range(ul)])
            BiBark[a] = sum([np.exp(b20)*float(d)**b21*DbhDist[a][d] for d in range(ul)])
            BiBranch[a] = sum([np.exp(b30)*float(d)**b31*DbhDist[a][d] for d in range(ul)])
            BiFoliage[a] = sum([np.exp(b40)*float(d)**b41*DbhDist[a][d] for d in range(ul)])
        RootMass = (BiStem + BiBark + BiBranch + BiFoliage) * BgAgRatio
        MerchMass = np.sum(MercDist, axis=1) * rho                            
        LAI = BiFoliage / 10000.0 * SLAI

        def smooth_bm(age, lon, dt, biomass):
            idx_juvenile=np.ravel(np.where(age/dt < lon))
            #idx_juvenile=np.ravel(np.where(self.WMeanDiam < 2.0))
            end = min(idx_juvenile[-1]+1, len(age))
            SI=biomass[end]; Iage=age[end]
            #print lon, end, SI, Iage            
            #from scipy.interpolate import interp1d
            #f=interp1d(np.array([0, end]), np.array([0, SI]))
            #biomass[idx_juvenile] = f(idx_juvenile)            
            p0 = 0.16; p1 = 3.3
            biom2 = lambda a, b1,b2: SI * ((1.0 - np.exp(-1.0 * b1 * a)) / (1.0 - np.exp(-1.0 * b1 * Iage))) ** b2
            biomass[idx_juvenile] = biom2(age[idx_juvenile],p0, p1)
            return biomass
        
        BiFoliage = smooth_bm(age,int(LeafLon/dt), dt, BiFoliage )    
        BiBark = smooth_bm(age,int(BarkLon/dt), dt, BiBark )    
        BiBranch = smooth_bm(age,int(BranchLon/dt), dt, BiBranch )    

        #---------Indicidual tree biomass growth taking into account the litterfall kg/tree/timestep (lonegevity of biomass componen)     
        def solve_Gr(n, lon, biom=None):
            """
            IN:             
                n - number of time steps, 
                lon - longevity of biomass component in timesteps,
                biom - amount of biomass component in one tree
            OUT:
                actual growth rate of biomass component in kg/tree/timestep
            """
            A = np.zeros((n,n))
            i,j = np.indices(np.shape(A))
            for k in range(lon):
                A[i==j+k]=1                               
            return np.linalg.multi_dot([np.linalg.inv(A),biom])                        # solve the true rate of biomass growth
        nsteps=self.nsteps
        FoGrTree = solve_Gr(nsteps, int(LeafLon/dt), BiFoliage/Survival)
        BaGrTree = solve_Gr(nsteps, int(BarkLon/dt), BiBark/Survival)
        BrGrTree = solve_Gr(nsteps, int(BranchLon/dt), BiBranch/Survival)
        RootGrTree = solve_Gr(nsteps, nsteps, RootMass/Survival)
        FineRootGrTree = solve_Gr(nsteps, int(FineRLon/dt), (RootMass*FineRootPr)/Survival)
        StemGrTree = solve_Gr(nsteps, nsteps, BiStem/Survival)
        
        #------Append results to instance variables -----------------
        self.BiStem   = np.append(self.BiStem, BiStem);     self.BiBark    = np.append(self.BiBark, BiBark)
        self.BiBranch = np.append(self.BiBranch, BiBranch); self.BiFoliage = np.append(self.BiFoliage,BiFoliage)
        self.RootMass = np.append(self.RootMass, RootMass); self.MerchMass = np.append(self.MerchMass, MerchMass)
        self.LAI = np.append(self.LAI, LAI)
        self.FoGrTree = np.append(self.FoGrTree, FoGrTree); self.BaGrTree = np.append(self.BaGrTree, BaGrTree)
        self.BrGrTree = np.append(self.BrGrTree, BrGrTree); self.RootGrTree = np.append(self.RootGrTree, RootGrTree)
        self.FineRootGrTree = np.append(self.FineRootGrTree, FineRootGrTree); 
        self.StemGrTree = np.append(self.StemGrTree, StemGrTree)
        
        del BiStem, BiBark, BiBranch, BiFoliage, RootMass, MerchMass, LAI 
        del FoGrTree, BaGrTree, BrGrTree, RootGrTree, FineRootGrTree, StemGrTree 

    def fLitter(self, rowGen = 0):
        """
        Computes litterfall from longevity (life span, years) and quantity of biomass components.  
        Makes distinction between litter generated from living trees and that from dead stems.
        Composes continous function from acculumation of biomass fractions, computes its derivative, and uses it with 
        biomass fraction life span to compute the litter output. \n
        Output: \n
            FoLitL - foliage litter from living trees (kg/ha/time step) \n
            BaLitL - bark litter from living trees (kg/ha/timestep) \n
            BrLitL - branch litter from living trees (kg/ha/timestep) \n
            FineRLitL - fine root litter from living trees (kg/ha/timestep) \n
            \n
            FoLitD - foliage litter from dead trees (kg/ha/timestep) \n
            BaLitD - bark litter from dead trees (kg/ha/timestep) \n 
            BrLitD - branch litter from dead trees (kg/ha/timestep) \n
            RootLitD - root litter from dead trees (kg/ha/timestep) \n
            CWD - stems (coarse woody debris, kg/ha/timestep)
        Kwargs: \n
            rowGen  - row index in parameterfile General parameters
        """
        #-----Parameters------------------------- 
        print ('    + Litter')        
        rg = rowGen
        BarkLon = self.g['BarkLon'][rg]; BranchLon = self.g['BranchLon'][rg]; LeafLon = self.g['LeafLon'][rg]
        FineRLon = self.g['FineRLon'][rg]; FineRootPr = self.g['FineRootPr'][rg]; nsteps = int(len(self.age)/(self.rotation +1))
        dt = self.g['dt'][0]/12.0; Nrotat = int(self.g['Nrotat'][0])   
        age = self.age[-nsteps:]; Survival = self.Survival[-nsteps:]; BiFoliage = self.BiFoliage[-nsteps:]
        BiBark = self.BiBark[-nsteps:]; BiBranch = self.BiBranch[-nsteps:]; RootMass = self.RootMass[-nsteps:]
        MerchMass = self.MerchMass[-nsteps:]
        
        #-----compute litter production from living trees
        FoGrTree=self.FoGrTree[-nsteps:]; BaGrTree=self.BaGrTree[-nsteps:] 
        BrGrTree=self.BrGrTree[-nsteps:]; RootGrTree = self.RootGrTree[-nsteps:]
        FineRootGrTree = self.FineRootGrTree[-nsteps:]
        
        FoLitL = shift(FoGrTree,int(LeafLon/dt),cval=0.0)*Survival              # litterfall kg/ha/timestep
        BaLitL = shift(BaGrTree,int(BarkLon/dt),cval=0.0)*Survival
        BrLitL = shift(BrGrTree,int(BranchLon/dt),cval=0.0)*Survival              
        FineRLitL = shift(FineRootGrTree,int(FineRLon/dt),cval=0.0)*Survival

        #-------compute litter production caused by mortality
        fsur = interS(age, Survival, k=1)
        dS = fsur.derivative()    #rate of mortailty
        M = dS(age)*dt*-1.0  #Number of dead trees in time step
        FoLitD = BiFoliage * M / Survival
        BaLitD = BiBark * M / Survival
        BrLitD = BiBranch * M / Survival
        RootLitD = RootMass * M / Survival
        CWD = MerchMass * M / Survival

        #-------add logging residues to litter production, initialization and then at the end of rotation ------------
        optLResid = True   #in Ireland, no previous stand    
        if optLResid:
            if self.rotation ==0: 
                FoLitD[0] = BiFoliage[-1] 
                BaLitD[0] = BiBark[-1] * 0.3
                BrLitD[0] = BiBranch[-1]
                RootLitD[0] = RootMass[-1] * 0.5                                    #stumps & roots
                CWD[0] = sum(CWD)*0.5  + RootMass[-1] * 0.5
            if self.rotation < Nrotat-1:
                FoLitD[nsteps-1] = BiFoliage[-1] 
                BaLitD[nsteps-1] = BiBark[-1] * 0.3 
                BrLitD[nsteps-1] = BiBranch[-1]
                RootLitD[nsteps-1] = RootMass[-1] * 0.5
                CWD[nsteps-1] = RootMass[-1] * 0.5
                    
        #-------append results into instance variables --------------
        self.FoLitL = np.append(self.FoLitL, FoLitL); self.BaLitL    = np.append(self.BaLitL, BaLitL) 
        self.BrLitL = np.append(self.BrLitL, BrLitL); self.FineRLitL = np.append(self.FineRLitL, FineRLitL)
        self.FoLitD = np.append(self.FoLitD, FoLitD); self.BaLitD    = np.append(self.BaLitD, BaLitD) 
        self.BrLitD = np.append(self.BrLitD, BrLitD); self.RootLitD  = np.append(self.RootLitD, RootLitD)
        self.CWD    = np.append(self.CWD, CWD)

        del FoLitL, FoLitD, BaLitL, BaLitD, BrLitL, BrLitD, FineRLitL, RootLitD, CWD
        
    def fWeeds(self, rowGrowth=0, rowGen=0):
        """
        Computes weed biomass, LAI, and litterfall. Allows setting the maximum green mass (stand foliage + weeds) \n
        limit; if exceeded, weed mass is reduced. Weeding is given as time list (time in months); in weeding \n 
        all weed biomass is located to litter from dead plants (no retranslocation). Litter produced  \n
        because of weed longevity is located to litter from living plants (retranslocation allowed):  \n
            Kwargs: \n 
            - rowGrowth - row index in parameter file growth and yield page \n
            - rowGen - row index in parameter file general parameters \n            
            Input: \n
            - maxGreenMass - maximum allowed sum of tree stand foliage and weed above ground mass, kg ha-1  \n
            - weedings - time of weeding, as list, month at when weeding is applied \n
            - appAge - apparent age in yrs, np.array (adjusted age because of Type 1 growth response) \n
            - BiFoliage  - mass of tree stand foliage, kg ha-1 \n
            Output: \n
            - weedAbove - above ground weed biomass, kg ha-1
            - weedBelow - below ground weed biomass, kg ha-1
            - weedALitL - weed litter from living above ground plant parts, allows retranslocation, kg ha-1 timestep-1
            - weedALitD - weed litter from dead above ground plant parts, no retranslocation, kg ha-1 timestep-1
            - weedBLitL - weed litter from living below ground plant parts, allows retranslocation, kg ha-1 timestep-1
            - weedBLitD - weed litter from dead below ground plant parts, no retranslocation, kg ha-1 timestep-1
            - weedLAI - weed leaf area index, m2 m-2
                     
        """
        print ('    + Weeds')
        #-----These to arguments--------------
        if len(self.weedings)==0:        
            self.weedings =np.array([0,3,6,9,18])                               #give the age when weeding is done, in months         
        weedings = self.weedings
        age = self.age; nsteps = int(len(age)/(self.rotation +1))
        appAge = self.appAge[-nsteps:]; BiFoliage = self.BiFoliage[-nsteps:]
        weedALitD = np.zeros(nsteps); weedBLitD = np.zeros(nsteps)
        
        #--------Parameters---------------
        r=rowGrowth; rg = rowGen          
        maxGreenMass = self.g['maxGreenMass'][rg]                               #maximum green mass includes above ground weed mass and tree stand foliage mass         
        self.maxGreenMass = maxGreenMass
        b0 = self.p['Wbeta0'][r];b1 = self.p['Wbeta1'][r]                       #Weed biomass growth parameters, result in kg/ha
        sa = self.g['WeedLeafA'][rg]; rsr = self.g['WeedRSRatio'][rg]           #Specific leaf area weeds kg/m2, root to shoot-ratio
        wl = self.g['WeedsLon'][rg]                                             #weed logevity
        dt = self.g['dt'][0]/12.0; Nrotat = int(self.g['Nrotat'][0])   
        a_age = np.arange(dt,dt*nsteps+dt, dt)                                  #age between 0...end of rotation yrs
        fW = lambda a: b0 * np.exp(-b1/a)       

        #----Construction of weed age series, in weeding set to 0. Apparent age is used in the growth computation, allows consuderation of fertilization,
        a1 = np.diff(appAge) ; a1 = np.insert(a1, 0, a1[0])        
        n=0; aa=0; tmpAges=[]; k=0
        for a in range(nsteps):
            if n<len(weedings): 
                if a==weedings[n]:
                    aa=0; n+=1
                else:
                    aa+=a1[k]
                tmpAges.append(aa)
            else:
                aa+=a1[k]
                tmpAges.append(aa)                            
        tmpAges = np.array(tmpAges) + 0.001  #to avoid zero division       
  
        #---weed biomass computation ----------------------
        WLAI = fW(tmpAges)*sa/10000.0                                           # Leaf area index m2/m2
        weedAbove = fW(tmpAges)                                                 # Above ground weed mass kg/ha

        weedBelow = weedAbove * rsr                                             # Below ground weed mass kg/ha
        potTot = weedAbove + BiFoliage                                          # potential total green biomass in site
        cutTot = np.minimum(potTot, maxGreenMass)                               # cut the fraction exceeding the maximum green mass
        cutWeed = np.where(cutTot-potTot < 0, cutTot-potTot, 0.0)               # the exceeding part is subtracted from weed mass
        weedAbove = np.maximum(weedAbove + cutWeed, 0.0)                        # adjust the weed mass
        weedBelow = weedAbove * rsr

        t = interS(a_age, weedAbove, k=1)
        #t = interS(tmpAges, weedAbove, k=1)        
        dW = t.derivative()
       #-------Litter production by age and longevity (retranslocation allowed) ---------------
        weedALitL = self.computeWeedLitter(age, nsteps, wl, dt, weedAbove)      # above ground litter production, kg ha-1 timestep-1
        weedBLitL = weedALitL * rsr                                             # below ground litter production, kg ha-a timestep-1
        weedDecline = np.where(dW(a_age)<0, dW(a_age), 0.0)*dt*-1.0
        
        weedALitL = weedALitL + weedDecline        
        weedBLitL = weedBLitL + weedDecline * rsr

        
        #-------In weeding, locate weed biomass to litter from dead plants (no retranslocation)
        n=0
        for a in range(nsteps):
            if a==weedings[n]:
                aa= tmpAges[a-1]+dt if a!=0  else 0.001
                weedALitL[a] = 0.0; weedBLitL[a] = 0.0 
                weedALitD[a] = fW(aa); weedBLitD[a] = weedALitD[a] * rsr
                n+=1
                if n==len(weedings):
                    break

        #-------Initialize litter with final biomass
        if self.rotation ==0: 
            weedALitD[0] = weedAbove[-1]
            weedBLitD[0] = weedBelow[-1]
        if self.rotation < Nrotat-1:
            weedALitD[nsteps-1] = weedAbove[-1] 
            weedBLitD[nsteps-1] = weedBelow[-1]
        

        #-------append results into instance variables --------------
        self.weedAbove =  np.append(self.weedAbove, weedAbove); self.weedBelow = np.append(self.weedBelow, weedBelow) 
        self.weedALitL = np.append(self.weedALitL, weedALitL); self.weedBLitL = np.append(self.weedBLitL, weedBLitL)
        self.weedALitD = np.append(self.weedALitD, weedALitD); self.weedBLitD = np.append(self.weedBLitD, weedBLitD)
        self.weedLAI = np.append(self.weedLAI, WLAI)
      

    def computeWeedLitter(self, age, nsteps, lon, dt, biomass):
        L = np.zeros(nsteps); lspan = int(lon/dt); lcohort=np.zeros(lspan)
        for a in range(nsteps):
            L[a] = lcohort[0]
            if biomass[a]==0.0:
                lcohort= lcohort*0.0
            lcohort=np.delete(lcohort,0)
            new =biomass[a]-np.sum(lcohort)
            lcohort = np.append(lcohort, max([new, 0.0]))
        return L

