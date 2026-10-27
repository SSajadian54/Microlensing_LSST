from __future__ import print_function
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarQuery
import dustmaps.bayestar
dustmaps.bayestar.fetch()
import numpy as np 
import astropy.units as units
bayestar = BayestarQuery(version='bayestar2019')
import matplotlib.pyplot as plt
import pylab as py
import math 
from pathlib import Path

import dustmaps.decaps
dustmaps.decaps.fetch()
from dustmaps.decaps import DECaPSQuery
decaps = DECaPSQuery(mean_only=True)



#####################################################################################################################


tt=np.zeros((100,2));  
for i1 in range(801): 
    i=i1+0
    for j in range(161): 
        lon=float(-100.0+0.25*i)
        lat=float(-20.0+0.25*j )
        if(lon<0.0): lon2=float(360.0+lon)
        else:        lon2=float(lon)
        #print("***********************************************************")
        #print("i, j, lon,   lat,  lon2: ",  i,   j,   lon,   lat ,   lon2  )
        
  
        #file_path=Path("./ExtBystar/Ext{0:.2f}_{1:.2f}.dat".format(float(lat),float(lon2)))
        #if(file_path.is_file()):   print("The file exists.")
        #elif(lat<=-10.0 or lat>=10.0):
        if(i==286 and j==72): 
            flag=1
            for k in range(100):
                dis=float(0.25+k*0.5)*1000.0## parsec        
                coords=SkyCoord(lon2*units.deg, lat*units.deg, distance=dis*units.pc, frame='galactic')
                #E =bayestar(coords, mode='median') ##Bayestar
                #E= decaps(coords)*(3.32/8.35816659571)## DECapS
                
                if(math.isfinite(E)): 
                    tt[k,0], tt[k,1]= float(dis/1000.0),  float(E*0.3026)## Ks_Extinction
                if(not math.isfinite(E)  and k>1):
                    dist=float(dis/1000.0);#kpc  
                    if(dist<15.0 and k<30): m=abs(tt[k-1,1]-tt[k-2,1])/(0.5)*abs(30.0-k)/30.0
                    else:                   m=0.0      
                    tt[k,0],tt[k,1]=dist , float( m*0.5+tt[k-1,1] )     
                if(not math.isfinite(E) and k<2): 
                    k=100
                    flag=-1  
            if(flag>0):   
                filew=open("./Ext{0:.2f}_{1:.2f}.dat".format(float(lat),float(lon2)),"w")           
                for k in range(100):
                    np.savetxt(filew,tt[k,:].reshape((-1,2)),fmt="%.2f  %.4f")
                filew.close()  
                
            print(flag)
#####################################################################################################################      







