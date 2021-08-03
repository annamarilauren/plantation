# -*- coding: utf-8 -*-
"""
Created on Mon Feb  28 2021
Main program for plantation simulator. Uses modules:
\t Parameters
\t Growth
\t Decomposition
\t Nutrient
\t hydro_tropic
\t Figures

@author: Annamari (Ari) Lauren, 
University of Eastern Finland 
(Natural Resources Institute Finaland)
"""
import numpy as np
import datetime
import matplotlib.pylab as plt
import Parameters as para
import Growth as gr
import Decomposition as dec
import Nutrient as nut
import Figures as fgs
import Results as res
from hydro_tropic import run_striphy

class cPlantationSimulator():
    def __init__(self,si, rhoIni, meanGwl, sdGwl, Mrate, parafile = None, fold = None, peatdatafile = None):

        #******** Files and folders**************
        self.fold = fold                                                       # work folder
        self.parafile = parafile                                               # Primary parameterfile
        self.peatdatafile = peatdatafile                                       # excel datafile for hydraulic characteristis of peat
        #******** Switches*********************
        self.optStripHy = False                                                # water table caclulation using striphy simulation
        self.optNut = True                                                     # compute nutrition
        self.DrawOpt = True                                                    # draw results
        self.start = datetime.datetime(2005,5,1)                               # start time for simulation, this matches the hydrological simulation results to same time
        self.meanDwt = meanGwl                                                 # if Striphy is not used: mean water table, cm
        self.dwtAmp = sdGwl                                                    # if Striphy not used: amplitude of dwt, cm        
        self.si = si                                                           # site index in m at index age (given in parameter file)
        self.rhoIni = rhoIni                                                   # initial bulk density of peat kg m-3
        self.Mrate = Mrate                                                     # constant mortality rate stems/ month (or rather time step)
        
    def run(self, weedings=[0,3,6,9]):
        """
        Parameters
        ----------
        weedings : TYPE, optional
            DESCRIPTION. The default time for weeding is [0,3,6,9] months.

        Returns
        -------
        bal : TYPE
            DESCRIPTION Carbon balance for the run
        """
        #************Read parametrs and create parameters instance*******************
        p=para.cParameters(pfile=self.parafile') 
        nsteps =  int((p.Gen['Lrotat'][0]*12 / p.Gen['dt'][0] +1) *p.Gen['Nrotat'][0])   #  numbre of time steps in the simulation
        
        #***********Overwrite parameters if given in function call, can be disabled ***********************************************    
 
        print (p.Gen['SI'])
        p.Gen['SI'][0]=self.si                                                 # overwrite site index if given in function call
        p.Gen['Mrate'][0]=self.Mrate
        p.Dec['peatRhoInit'][0]=self.rhoIni
        p.Dec['peatRhoFinal'][0]=self.rhoIni + p.Gen['Lrotat'][0]*p.Gen['Nrotat'][0]*3.
        
        #************Growth and yield**************************
        
        gy = gr.cGrowthAndYield(p.GY, p.Gen, p.Bi, p.Thin)                     # initialize growth and yoiel instance
        gy.weedings=weedings                                                   # overwrite weedings
        gy.calibrateBA = 1.0                        
        rowGe = [0,0,0,0,0,0]                                                  # row in "General" sheet where the parameters are picked for 
                                                                               # each rotation in simulation [0,1,2] would imply first rotation from first row, secon from secon row...
        rowGr = [0,0,0,0,0,0]                                                  # row in "Growth and yiueld" sheet where the parameters are picked for                                                                            
        gy.fComputeGrowthAndYield(rowGe=rowGe, rowGr=rowGr)                    # calculation                        
        
        #************Hydrology and water tables ***********************************
        
        if self.optStripHy: 
            dfStriphy=run_striphy(gy.Hdom,                                     # domainat height array
                                  gy.LAI+gy.weedLAI,                           # total leaf area
                                  int(p.Gen['Lrotat'][0]),                     # length of roation, years
                                  optFig=self.DrawOpt)                         # option to draw hydrology figures

        #**********Decomposition of organic matter*******************************************************
        de = dec.DecomCWD(p.Gen,                                               # initialize decomposition instance with parameters and growth and yield instance
                          p.Dec, 
                          gy, 
                          p.N, 
                          p.P, 
                          p.K)                     
                                                                               # run coarse woody debris function
        de.decomCWD(gy.CWD,                                                    # CWD input
                    gy.rowGe,                                                  # parameters
                    dbh = gy.WMeanDiam.copy(),                                 # mean stem diameter as time series
                    rot = gy.rotationArr)                                      # number of rotations in the simulation
        
                                                                               # initialize Romul: decomposition of stand and weed litter
        ro = dec.Romul(p.Dec,                                                  # parameters 
                       gy,                                                     # growth and yield instance with litter input
                       de,                                                     # CWD instance
                       p.N,                                                    # nutrient contents and other nutrient parameters
                       p.P, 
                       p.K, 
                       test=False)                                             # option to test the Romul model
        ro.decomRomul()                                                        # run Romul
        
        if self.optStripHy:                                                    # peat decomposition if hydrology is calculated
            pe = dec.DecomPeat(len(gy.CWD),                                    # number of timesteps
                               1,                                              # dt, time step lenth in months
                               dwtArr=dfStriphy['gwl'].values*100.0)           # array of water tables in cm
            pe.decomposePeat(p.Dec,                                            # run peat decomposition, parameters
                             gy,                                               # litter inputs 
                             de,                                               # CWD instance
                             ro,                                               # Romul instance
                             optH=False, optCWD = True, optSine = False)       # options
        else:                                                                  # peat decomposition if water tables are exogeneous input
            dwtArr = np.random.normal(self.meanDwt, self.dwtAmp,nsteps)        # generate monthly WTs    
            pe = dec.DecomPeat(len(gy.CWD),                                    # number of timesteps 
                               1,                                              # dt, length of timestep in months
                               gwl=self.meanDwt,                               # given mean WT
                               dwtAmp=self.dwtAmp,                             # given amplitude for WT
                               dwtArr=dwtArr,                                  # WT as array  
                               peatdatafile=self.peatdatafile)                 # datafile for peat hydraulic characteristics
            pe.decomposePeat(p.Dec, gy, de, ro, optH=False, optCWD = True, optSine = False)

#**********************************************************************************************************************************

        if self.optNut == True:                                                # calclucate nutrient balance
            print ('**********************') 
            print ('NUTRIENT BALANCE')
            N = nut.cNutrient("Nitrogen",                                      # initialize nutrient for nitrogen
                              p.N,                                             # N parameters
                              p.Gen,                                           # general parameters
                              p.Fer,                                           # fertilization
                              ro)                                              # litter (Romul) instance
            N.nutDemand(gy)                                                    # calculate N demand   
            N.nutSupply(gy, de, pe)                                            # calculate N supply
            
            P = nut.cNutrient("Phosphorus",
                              p.P, 
                              p.Gen, 
                              p.Fer, 
                              ro)
            P.nutDemand(gy)    
            P.nutSupply(gy, de, pe)
            
            K = nut.cNutrient("Potassium",
                              p.K, 
                              p.Gen, 
                              p.Fer, 
                              ro)
            K.nutDemand(gy)    
            K.nutSupply(gy, de, pe)

#*****************Outputs*************************************************************************
        if self.DrawOpt== True:
            plt.close('all')
            if self.optNut == True:            
                fgs.fDrawNutrient(N, gy, de)
                fgs.fDrawNutrient(P, gy, de)
                fgs.fDrawNutrient(K, gy, de)
                fgs.fDrawTreeNutBal(N,P,K, gy,de)
                
            fgs.fDrawDecomposition(gy, de, ro, pe)            
            fgs.fDrawGrowth(gy)
            fgs.fDrawMAICAI(gy)
            fgs.drawCderivative(gy, de, ro, pe)
                    
        if self.optNut==True:
            res.save_nbal(gy,pe,N,P,K,fold=self.fold)
            res.save_results(gy,de,ro,pe,N,P,K,fold=self.fold)
 
        _, bal, _ = res.C_balance(gy, de, ro, pe)        
        
        above_stand = gy.BiBark + gy.BiBranch + gy.BiFoliage + gy.BiStem
        print ('********** Nutrients and biomass ************************')
        print ('      Biom stand above', above_stand[-1], 'kg ha-1')
        print ('      Biom weeds above', max(gy.weedAbove), 'kg ha-1')
        print ('      Stand above N', (N.barkSto+N.branchSto+N.foliageSto+N.stemSto + N.rootSto)[-1])
        print ('      Stand above P', (P.barkSto+P.branchSto+P.foliageSto+P.stemSto + P.rootSto)[-1])
        print ('      Stand above K', (K.barkSto+K.branchSto+K.foliageSto+K.stemSto + K.rootSto)[-1])

        print ('*********Litter inputs kg **************************************')
        print ('      Root Litter', sum(gy.FineRLitL + gy.RootLitD))
        print ('      AG Litter', sum(gy.BrLitD + gy.BrLitL + gy.BaLitD + gy.BaLitL + gy.FoLitD + gy.FoLitL))
        print ('      Stemlitter',  sum(gy.CWD))         
        print ('      Weed litter', sum(gy.weedALitD + gy.weedALitL + gy.weedBLitD + gy.weedBLitL))

        print ('********** Stand nutrient uptake kg ************************')
        print ('      net N', sum(N.barkNet+N.branchNet+N.foliageNet+N.stemNet + N.rootNet))
        print ('      net P', sum(P.barkNet+P.branchNet+P.foliageNet+P.stemNet + P.rootNet))
        print ('      net K', sum(K.barkNet+K.branchNet+K.foliageNet+K.stemNet + K.rootNet))

        print ('      gross N', sum(N.barkGro+N.branchGro+N.foliageGro+N.stemGro + N.rootGro))
        print ('      gross P', sum(P.barkGro+P.branchGro+P.foliageGro+P.stemGro + P.rootGro))
        print ('      gross K', sum(K.barkGro+K.branchGro+K.foliageGro+K.stemGro + K.rootGro))
        
        print ('********** Weed nutrient uptake kg ************************')
        print ('      gross N', sum(N.weedAboveGro + N.weedBelowGro))
        print ('      gross P', sum(P.weedAboveGro + P.weedBelowGro))
        print ('      gross K', sum(K.weedAboveGro + K.weedBelowGro))
        
        
        bals = res.C_balance(gy, de, ro, pe)
        sim_time_yrs = p.Gen['Lrotat'][0] * p.Gen['Nrotat'][0]
        print ('************ C balances kg ha-1 yr-1 **********************************')
        print ('      Soil', bals[0]/sim_time_yrs  )
        print ('      Stand without stems ',  bals[1]/sim_time_yrs )
        print ('      Total',        bals[2]/sim_time_yrs)
     
        

        print ('*********** Nutrient supply, N, P, K *********************')
        print ('      Deposition', sum(N.deposition), sum(P.deposition), sum(K.deposition))
        print ('      Weathering', sum(N.weathering), sum(P.weathering),sum(K.weathering)) 
        print ('      Mic fixing', sum(N.fixing), sum(P.fixing),sum(K.fixing)) 
        print ('      Retranslocation', sum(N.retranslocation), sum(P.retranslocation),sum(K.retranslocation)) 
        print ('      Retranslocation weeds', sum(N.retranslocationWeeds), sum(P.retranslocationWeeds),sum(K.retranslocationWeeds)) 
        print ('      Release CWD', sum(N.releaseCwd), sum(P.releaseCwd),sum(K.releaseCwd)) 
        print ('      Release Peat', sum(N.releasePeat), sum(P.releasePeat),sum(K.releasePeat)) 
        print ('      Release Litter', sum(N.releaseLitter), sum(P.releaseLitter),sum(K.releaseLitter)) 
        print ('      Release fertilizers', sum(N.ferRel), sum(P.ferRel),sum(K.ferRel)) 
        
        
        return bal

def runPs():
        fold = r'C:/Users/alauren/OneDrive - University of Eastern Finland/codes/plantation_simulator/'    
        parafile = r'C:/Users/alauren/OneDrive - University of Eastern Finland/codes/plantation_simulator/ps_params.xlsx'
        peatdatafile = r"C:/Users/alauren/OneDrive - University of Eastern Finland/codes/plantation_simulator/peat_data.xlsx"        
        rhoIni = 110.
        si = 21.0
        Mrate = 9.
        meanGwl = 40.
        sdGwl = meanGwl*0.2
        pl=cPlantationSimulator(si, rhoIni, meanGwl, sdGwl, Mrate, \
                                parafile=parafile, fold = fold, peatdatafile = peatdatafile)
        pl.run(weedings =[0,3,6,9,12])

runPs()
