# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 17:10:53 2016
Retrieves parameters from file and locates them into dictionary variables.
General parameters: self.Gen
        Age of the stand (years)        age	
        Index age, years	               Iage
        Hdom at index age, m	           SI
        N planted, number	          N
        Top diameter merchantability limit (cm)  Tdia
        Density of wood (kg/m3)	        WoodDens \n
        Density of bark (kg/m3)	        BarkDens \n
        Below / above ground mass proportion	BeAbRatio  \n
        Bark longevity (years)	         BarkLon  \n
        Branch longevity (years)	         BranchLon  \n
        Leaf longevity (years)	         LeafLon  \n
        Fine root longevity (years)	    FIneRLon  \n
        Fine root proportion of root mass (share)	FineRootPr \n
        Specific leaf area (m2/kg)	     SLAI  \n
        Number of weedings	             Nweedings \n
        Upper limit of diameter distribution (cm)	DiamUpperLimit  \n
        Run simulation until age (years)	RunUntil   \n
        Biomass model tree Species: A. crassicarpa April	BiomSpe \n 
        Weeds: specific leaf area       	WeedLeafA   \n
        Weeds: Below/Above ground biomass ratio	WeedRSRatio   \n
        Weeds:Leaf&root longevity (years)	  WeedsLon  \n
        Soil loss in harvesting 0…1	    SoilLoss  \n
        Length of time step (months)	       dt   \n
        Growth model tree species: Acacia crassicarpa APRIL	GrTreeS \n
        
Growth and Yield Parameters, self.GY:
        Diameter distribution (Two parameter Weibull Forss et al 1998)	
            Eq4Beta00
            Eq4Beta01
            Eq4Beta02
            Eq7Beta10
            Eq7Beta11
        Individual tree height 	
            Eq8Beta00
            Eq8Beta01
            Eq8Beta02
            Eq9Beta10
            Eq9Beta11
            Eq9Beta12
            Eq9Beta13
        Tree volume (over bark)	
            Eq10Beta0
            Eq10Beta1
            Eq10Beta2
        Merchantable volume (over bark)	
            Eq11Beta0
            Eq11Beta1
        Dominant height Eq 2	
            Eq2Beta1
            Eq2Beta2
        Eq 4 Survival	
            Eq4Beta0
            Eq4Beta1
            Eq4Beta2
            Eq4Beta3
        Eq 5 Basal area	
            Eq5Beta0
            Eq5Beta1
            Eq5Beta2
            Eq5Beta3
        Wood consumption	
            WoodBeta0
            WoodBeta1
        Weed LAI	
            WLAIBeta0
            WLAIBeta1

    Biomass functions
        Bi,H., Turner, J. & Lambert, M. 2004. Additive biomass functions for native eucalypt forest trees of temperate Australia.	
        Trees 18: 467-479.	
        Stem:  b10, b11
        Bark b21, b22 
        Branch b31, b32
        Foliage b41, b42
    
    Decomposition model (Romul) parameters:
    L refers to litter, F partly decomposed fermentation material, H humus material
    root - below ground material, above - foliage, branch, bark
    M - mass, N - nitrogen, P - phosphorus, K - potassium
    Units in input kg/ha, then converted in calculation into kg/m2
        rootLM	
        rootLN	
        rootLP	
        rootLK	
        aboveLM
        aboveLN	
        aboveLP
        aboveLK	
        rootFM	
        rootFN	
        rootFP	
        rootFK	
        aboveFM	
        aboveFN	
        aboveFP	
        aboveFK	
        rootHM	
        rootHN	
        rootHP	
        rootHK	
        abovHM	
        aboveHN	
        aboveHP	
        aboveHK
        
        aboveMass
        aboveN
        aboveP
        aboveK
        belowMass
        belowN
        belowP
        belowK
        
@author: Ari Lauren 2016
"""
import pandas as pd
import numpy as np
import datetime
import pickle
from dateutil.relativedelta import relativedelta

class cParameters():
    def __init__(self, pfile):
        ParamFile = pfile
        self.Gen = pd.read_excel(ParamFile, sheet_name='General', skiprows=1).to_dict()
        self.GY = pd.read_excel(ParamFile, sheet_name='Growth and Yield', skiprows=1).to_dict()
        self.Bi = pd.read_excel(ParamFile, sheet_name='Biomass', skiprows=1).to_dict()
        self.Dec = pd.read_excel(ParamFile, sheet_name='Decomposition', skiprows=1).to_dict()
        self.N = pd.read_excel(ParamFile, sheet_name='Nitrogen', skiprows=1).to_dict()
        self.P = pd.read_excel(ParamFile, sheet_name='Phosphorus', skiprows=1).to_dict()
        self.K = pd.read_excel(ParamFile, sheet_name='Potassium', skiprows=1).to_dict()
        self.Fer = pd.read_excel(ParamFile, sheet_name='Fertilization', skiprows=1)
        self.Time = pd.read_excel(ParamFile, sheet_name='Timeseries', skiprows=1)
        self.SurTs = pd.read_excel(ParamFile, sheet_name='SurvivalTs', skiprows=0, usecols='B')
        self.Thin = pd.read_excel(ParamFile, sheet_name='Thinning', skiprows=1, usecols='A:C', index_col='time')
    
    def composeDwtFromTs(self, time):
        from scipy.interpolate import InterpolatedUnivariateSpline as interS
        t = np.array(self.Time['time']); dwt = np.array(self.Time['wt'])*-1
        smooth = interS(t, dwt, k=1)
        return smooth(time)        
        
    def readDwt(self, f, start, nsteps):
        sd = pickle.load( open(f, "rb" ) )        
        rows=sd['rows']; cols=sd['cols']; dates=sd['dates']; ele=sd['ele']; hts=sd['hts']
        gwl = (ele[rows/2, cols/2] - hts[:,rows/2, cols/2])*-1
        gwlClose=(ele[2, 2] - hts[:,2, 2])*-1
        sdf=pd.DataFrame(dates, index=dates); sdf['gwl']=pd.Series(gwl, index=sdf.index); sdf.columns=['dates', 'gwl']
        sdf=sdf.drop('dates',1)
        sdf['gwl2']=gwlClose
        sdf=sdf.resample('M', how='mean') 
        end = start + relativedelta(months=+nsteps)
        #sdf=sdf[datetime.datetime(start):datetime.datetime(end)]
        sdf=sdf[start:end]        
        return sdf.index.values, sdf['gwl'].values, sdf['gwl2']    