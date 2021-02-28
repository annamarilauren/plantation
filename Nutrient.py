# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 10:23:52 2016
Describes nutrient uptake and nutrient supply scheme for different nutrients (N, P, K) 

@author: Ari Laurén

"""
from scipy.interpolate import InterpolatedUnivariateSpline as interS
import numpy as np

class cNutrient():
    def __init__(self, name, para, gen, fer, ro):
        self.name = name    #N, P, K
        self.para = para    # nutrient parameters from file
        self.gen = gen      # general parameters
        self.fer = fer      # fertilization time, dose and quality 
        self.ro = ro        # romul variables
        print ('    + ' +self.name)
    
    def nutDemand(self, gy):
        #------General parameters ------------        
        row = 0
        dt = gy.g['dt'][0]/12.0; Nrotat = int(gy.g['Nrotat'][0])
        lr = gy.g['Lrotat'][0]; dtInRot = int(lr/dt  + 1)

        #------Nutrient parameters - --------
        cWo =  self.para['concWood'][row]/100.0;	cLe= self.para['concLeaf'][row]/100.0;
        cBa = self.para['concBark'][row]/100.0; cRo = self.para['concRoot'][row]/100.0;
        cBr = self.para['concBranch'][row]/100.0; cWe = self.para['concWeeds'][row]/100.0
        retrans = self.para['retrans'][row]; retransWeeds = self.para['retransWeeds'][row]
 
        #--------Growth and yield variables--------------------
        time = gy.time; age = np.reshape(gy.age,(Nrotat, dtInRot)); 
        stem  = np.reshape(gy.BiStem,(Nrotat, dtInRot)) 
        bark  = np.reshape(gy.BiBark,(Nrotat, dtInRot)) 
        branch  = np.reshape(gy.BiBranch,(Nrotat, dtInRot)) 
        foliage  = np.reshape(gy.BiFoliage,(Nrotat, dtInRot)) 
        root  = np.reshape(gy.RootMass,(Nrotat, dtInRot)) #* (1.0-fineRootPr)
        survival = np.reshape(gy.Survival, (Nrotat, dtInRot))
        weedAbove = np.reshape(gy.weedAbove, (Nrotat, dtInRot))
        weedBelow = np.reshape(gy.weedBelow, (Nrotat, dtInRot))

        FoGrTree = np.reshape(gy.FoGrTree,(Nrotat, dtInRot))
        StemGrTree = np.reshape(gy.StemGrTree,(Nrotat, dtInRot))
        BrGrTree = np.reshape(gy.BrGrTree,(Nrotat, dtInRot))
        BaGrTree = np.reshape(gy.BaGrTree,(Nrotat, dtInRot))
        RootGrTree = np.reshape(gy.RootGrTree,(Nrotat, dtInRot))
        FineRootGrTree = np.reshape(gy.FineRootGrTree,(Nrotat, dtInRot))


        #-------Initialize result arrays, net uptake---------
        self.stemNet = np.empty(0); self.branchNet=np.empty(0)
        self.barkNet = np.empty(0); self.foliageNet = np.empty(0);
        self.rootNet = np.empty(0); self.weedAboveNet = np.empty(0);
        self.fineRootNet = np.empty(0); self.weedBelowNet = np.empty(0)
        
        self.stemGro = np.empty(0); self.branchGro=np.empty(0)
        self.barkGro = np.empty(0); self.foliageGro = np.empty(0);
        self.rootGro = np.empty(0); self.weedAboveGro = np.empty(0);
        self.fineRootGro = np.empty(0); self.weedBelowGro = np.empty(0)
        

        #-------Initialize storage arrays, storage---------
        self.stemSto = np.empty(0); self.branchSto=np.empty(0)
        self.barkSto = np.empty(0); self.foliageSto = np.empty(0);
        self.rootSto = np.empty(0); self.weedAboveSto = np.empty(0);
        self.fineRootSto = np.empty(0); self.weedBelowSto = np.empty(0)

        #-------Computation-----------------------------       
        
        def computeNet(b, a, co, dia=gy.WMeanDiam, survival=None, weed=False):
            #net per tree
            #b - biomass, a-age, co -concentration
            if survival is None: survival=np.ones(len(a))
            idx_seedling=[0,1,2,3,4] #np.where(dia<3.5)
            mass = interS(a, b/survival, k=1)
            dm = mass.derivative()
            rate = dm(a) * survival * dt
            if weed is False: rate[idx_seedling] = b[idx_seedling]*0.2*dt
            return np.maximum(co * rate, 0.0)
        
        
        for r in range(Nrotat):           
            frp = self.gen['FineRootPr'][int(gy.currentRowGe[r])]       
            self.stemNet = np.append(self.stemNet, computeNet(stem[r,:], age[r,:], cWo, survival=survival[r,:]))
            self.barkNet = np.append(self.barkNet, computeNet(bark[r,:], age[r,:], cBa, survival=survival[r,:]))
            self.branchNet = np.append(self.branchNet, computeNet(branch[r,:], age[r,:], cBr, survival=survival[r,:]))
            self.foliageNet = np.append(self.foliageNet, computeNet(foliage[r,:], age[r,:], cLe, survival=survival[r,:]))
            self.rootNet = np.append(self.rootNet, computeNet(root[r,:]*(1.0-frp), age[r,:], cWo, survival=survival[r,:]))
            self.fineRootNet = np.append(self.fineRootNet, computeNet(root[r,:]*frp, age[r,:], cRo, survival=survival[r,:]))       
            self.weedAboveNet =np.append(self.weedAboveNet, computeNet(weedAbove[r,:], age[r,:], cWe, weed=True))
            self.weedBelowNet =np.append(self.weedBelowNet, computeNet(weedBelow[r,:], age[r,:], cWe, weed=True))
        
        """
        import matplotlib.pylab as plt
        fig = plt.subplot(521); plt.plot(np.ravel(age), np.ravel(foliage), 'b-')
        plt.subplot(522); plt.plot(np.ravel(age), np.ravel(self.foliageNet), 'b-')
        plt.subplot(523); plt.plot(np.ravel(age), np.ravel(branch), 'r-')
        plt.subplot(524); plt.plot(np.ravel(age), np.ravel(self.branchNet), 'r-')
        plt.subplot(525); plt.plot(np.ravel(age), np.ravel(root), 'm-')
        plt.subplot(526); plt.plot(np.ravel(age), np.ravel(self.rootNet), 'm-')
        plt.subplot(527); plt.plot(np.ravel(age), np.ravel(survival), 'c-')

        plt.subplot(528); plt.plot(np.ravel(age), np.ravel(weedAbove), 'k-')
        plt.subplot(529); plt.plot(np.ravel(age), np.ravel(self.weedAboveNet), 'k-')

        import sys; sys.exit()
        """
        #----------Gross nutrient uptake, kg ha-1 timestep-1
        """        
        self.stemGro = self.stemNet  
        self.barkGro = (self.barkNet + gy.BaLitL * (1.0-retrans) * cBa).flatten()      
        self.branchGro = (self.branchNet + gy.BrLitL * (1.0-retrans) * cBr).flatten()
        self.foliageGro = (self.foliageNet + gy.FoLitL * (1.0 - retrans) * cLe).flatten()
        self.rootGro = self.rootNet
        self.fineRootGro = (self.fineRootNet + gy.FineRLitL * (1.0 - retrans) * cRo ).flatten()   
        """
        self.weedAboveGro = (self.weedAboveNet + gy.weedALitL * (1.0 - retransWeeds) * cWe).flatten()
        self.weedBelowGro = (self.weedBelowNet + gy.weedBLitL * (1.0 - retransWeeds) * cWe).flatten()

        #----------Gross nutrient uptake, kg ha-1 timestep-1
        sh = np.shape(survival)
        BaLitL = np.reshape(gy.BaLitL,sh)
        BrLitL = np.reshape(gy.BrLitL,sh)
        FoLitL = np.reshape(gy.FoLitL,sh)
        FineRLitL = np.reshape(gy.FineRLitL,sh)

        for r in range(Nrotat):           
            frp = self.gen['FineRootPr'][int(gy.currentRowGe[r])]       
            sg = np.maximum((StemGrTree[r,:] * survival[r,:] )* cWo, np.zeros(len(survival[r,:])))           
            self.stemGro = np.append(self.stemGro, sg)
            bg = np.maximum(((BaGrTree[r,:] * survival[r,:] - BaLitL[r,:] * (retrans))* cBa ), np.zeros(len(survival[r,:])))
            self.barkGro = np.append(self.barkGro, bg)
            brg =np.maximum(((BrGrTree[r,:] * survival[r,:] - BrLitL[r,:] * (retrans))* cBr), np.zeros(len(survival[r,:])))
            self.branchGro = np.append(self.branchGro, brg)
            fg =np.maximum(((FoGrTree[r,:] * survival[r,:] - FoLitL[r,:] * (retrans))* cLe), np.zeros(len(survival[r,:])))
            self.foliageGro = np.append(self.foliageGro, fg)
            rg = np.maximum(((RootGrTree[r,:] * survival[r,:] )* cWo), np.zeros(len(survival[r,:])))
            self.rootGro = np.append(self.rootGro, rg)
            fgg = np.maximum(((FineRootGrTree[r,:] * survival[r,:] - FineRLitL[r,:] * (retrans))* cRo), np.zeros(len(survival[r,:])))
            self.fineRootGro = np.append(self.fineRootGro, fgg)       
        
    
        """
        self.stemGro = self.stemNet  
        self.barkGro = (self.barkNet + gy.BaLitL * (1.0-retrans) * cBa).flatten()      
        self.branchGro = (self.branchNet + gy.BrLitL * (1.0-retrans) * cBr).flatten()
        self.foliageGro = (self.foliageNet + gy.FoLitL * (1.0 - retrans) * cLe).flatten()
        self.rootGro = self.rootNet
        self.fineRootGro = (self.fineRootNet + gy.FineRLitL * (1.0 - retrans) * cRo ).flatten()   
        self.weedAboveGro = (self.weedAboveNet + gy.weedALitL * (1.0 - retransWeeds) * cWe).flatten()
        self.weedBelowGro = (self.weedBelowNet + gy.weedBLitL * (1.0 - retransWeeds) * cWe).flatten()
        """
        self.totalDemand = np.maximum(self.stemGro + self.branchGro + self.barkGro + \
                self.foliageGro + self.fineRootGro + self.rootGro + \
                self.weedAboveGro +self.weedBelowGro, np.zeros(len(self.stemGro)))

        self.treeDemand = np.maximum(self.stemGro + self.branchGro + self.barkGro + \
                self.foliageGro + self.fineRootGro + self.rootGro, np.zeros(len(self.stemGro)))
                
        #----------Nutrient storage in vegetation, kg ha-1
        for r in range(Nrotat):           
            frp = self.gen['FineRootPr'][int(gy.currentRowGe[r])]       
            self.stemSto = np.append(self.stemSto, stem[r,:] * cWo)
            self.barkSto = np.append(self.barkSto, bark[r,:] * cBa)
            self.branchSto = np.append(self.branchSto, branch[r,:] * cBr)
            self.foliageSto = np.append(self.foliageSto, foliage[r,:] * cLe)
            self.rootSto = np.append(self.rootSto, root[r,:]*(1.0-frp) * cWo)
            self.fineRootSto = np.append(self.fineRootSto, root[r,:]* frp * cRo)       
            self.weedAboveSto =np.append(self.weedAboveSto, weedAbove[r,:] * cWe)
            self.weedBelowSto =np.append(self.weedBelowSto, weedBelow[r,:] * cWe)


    def nutSupply(self, gy, de, pe):
        #------General parameters ------------        
        row = 0
        dt = gy.g['dt'][0]/12.0; Nrotat = int(gy.g['Nrotat'][0])
        lr = gy.g['Lrotat'][0]; dtInRot = int(lr/dt  + 1)
        micr = {'Nitrogen': 0.5, 'Phosphorus': 1.0, 'Potassium': 1.0}        
        
        #-------Nutrient parameters--------
        depo = self.para['depo'][row]; weathering = self.para['weathering'][row]    #kg ha-1 yr-1
        fixShare = self.para['fixShare'][row]; fixShareWeeds = self.para['fixShareWeeds'][row]  
        leachShare = self.para['leach'][row]        
        cWo =  self.para['concWood'][row]/100.0;	cLe= self.para['concLeaf'][row]/100.0;
        cBa = self.para['concBark'][row]/100.0; cRo = self.para['concRoot'][row]/100.0;
        cBr = self.para['concBranch'][row]/100.0; cWe = self.para['concWeeds'][row]/100.0
        retrans = self.para['retrans'][row]; retransWeeds = self.para['retransWeeds'][row]        
        cPe = self.para['peatCont'][row]

        #------Constants, kg ha-1 timestep-1--------------
        self.deposition = np.ones(Nrotat*dtInRot)*depo*dt
        self.weathering = np.ones(Nrotat*dtInRot)*weathering*dt
        self.fixing = fixShare * (self.stemGro + self.barkGro + self.branchGro 
                        + self.foliageGro + self.rootGro + self.fineRootGro)
        self.fixingWeeds = fixShareWeeds * (self.weedAboveGro + self.weedBelowGro)       
        self.retranslocation = retrans * (gy.FoLitL*cLe + gy.BaLitL*cBa 
                + gy.BrLitL*cBr + gy.FineRLitL*cRo)
        self.retranslocationWeeds = retransWeeds * (gy.weedALitL + gy.weedBLitL) * cWe
        self.releaseFert(gy)

        self.releasePeat = pe.peatCrelease * pe.CtoMass * cPe / 100.0 * micr[self.name]
        self.releaseCwd = de.cwdNutRelease[self.name]
        if self.name == 'Nitrogen': self.releaseLitter = self.ro.Nminer*10000.0
        if self.name == 'Phosphorus': self.releaseLitter = self.ro.Pminer*10000.0
        if self.name == 'Potassium': self.releaseLitter = self.ro.Kminer*10000.0
        
      
        self.totalSupply = self.deposition + self.weathering + self.fixing + \
                self.retranslocation + self.retranslocationWeeds + \
                self.releaseCwd +  self.releasePeat + self.releaseLitter +self.ferRel
        
        treeCoverage = np.minimum(gy.BiFoliage / gy.maxGreenMass, np.ones(len(gy.BiFoliage)))
        self.treeSupply = (self.deposition + self.weathering  + \
                self.releaseCwd +  self.releasePeat + \
                self.ferRel + self.releaseLitter - self.weedAboveGro - self.weedBelowGro) * treeCoverage + \
                self.fixing + self.retranslocation 
#        self.leaching = self.totalSupply * leachShare
    def releaseFert(self, gy):
        """
         fertilizer parameters: N (%), P2O5 (%), K2O (%),
             release kinetics parameters (month^-1): Kn, Kp, Kk
        """

        nStems = gy.Survival[0]                                                           # number of stems /ha
        dt = gy.g['dt'][0]; Nrotat = int(gy.g['Nrotat'][0])                     # time step (months), number of rotations
        lr = gy.g['Lrotat'][0]*12.0; dtInRot = lr/dt                            # length of rotation (months), timesteps in rotation
        length = int(Nrotat * dtInRot * dt + Nrotat)                                 # lenght of array including all timesteps   
        self.ferSto = np.zeros(length); self.ferRel = np.zeros(length)          # undissonlved fertilizer in soil kg ha-1, nutrient release  kg ha-1 timestep-1
 
        ferSto = np.zeros(length)                                               # local: fertilizer nutrient storage        
        ferRel = np.zeros(length)        
        time = np.arange(0,length) 
        
        if self.name == 'Nitrogen': name='N'; k = 'Kn'; conv = 1.0
        if self.name ==  'Phosphorus': name='P'; k = 'Kp'; conv = 0.42
        if self.name == 'Potassium': name='K'; k = 'Kk'; conv = 0.83
     
        fpara ={
             'NPK':{'N': 23, 'P':13, 'K': 7, 'Kn': 0.94, 'Kp':0.04, 'Kk': 0.3},
             'Suburin':{'N': 7, 'P':15, 'K': 8, 'Kn': 0.94, 'Kp':0.04, 'Kk': 0.3},
             'Urea':{'N': 45, 'P':0, 'K': 0, 'Kn': 0.9, 'Kp':0, 'Kk': 0},
             'RP':{'N': 0, 'P':9, 'K': 0, 'Kn': 0, 'Kp':0.04, 'Kk': 0},
             'Kcl':{'N': 0, 'P':0, 'K': 60, 'Kn': 0, 'Kp':0, 'Kk': 0.8},
             'Ash':{'N': 0, 'P':2.8/0.42, 'K': 7.6/0.83, 'Kn': 0, 'Kp':0.01, 'Kk': 0.01},
             }

        for r in range(Nrotat):                                                 # rotation
            for t in range(len(self.fer['time'])):                              # times in input data frame
                m = int(self.fer['time'][t] + (r *lr))                          # fertilization month
                for a in self.fer.keys():                                       # fertilizers in use
                    if self.fer[a][t] !=0 and a != 'time':                      # dose is > 0 
                        feNut = fpara[a][name]; kk = fpara[a][k]
                        dose = self.fer[a][t]
                        if kk > 0:                                              # nutrient is in the fertilizer
                            #print feNut, self.name, conv, nStems, kk
                            feNut = feNut * conv /100. * dose * nStems / 1000.0
                            ts = time[m:] - m
                            ferSto[m:] = ferSto[m:] + feNut * np.exp(-kk*ts*dt)
                            ferRel[m:] = ferRel[m:] + feNut*(1-np.exp(-kk*ts*dt))
        
        func = interS(time, ferRel, k=1)
        dfdt = func.derivative()
        self.ferRel = dfdt(time) * dt

        self.ferSto = ferSto
        