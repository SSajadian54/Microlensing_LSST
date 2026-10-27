from __future__ import print_function
from astropy.coordinates import SkyCoord
import numpy as np 
import astropy.units as units
import matplotlib.pyplot as plt
import pylab as py
from pathlib import Path
import math 

import dustmaps.decaps
dustmaps.decaps.fetch()
from dustmaps.decaps import DECaPSQuery
decaps = DECaPSQuery(mean_only=True)




#####################################################################################################################
tt=np.zeros((100,2));  
for i1 in range(470): 
    i=i1+0
    for j in range(161): 
        lon=float(-100.0+0.25*i)
        lat=float(-20.0+0.25*j )
        if(lon<0.0): lon2=float(360.0+lon)
        else:        lon2=float(lon)
        print("*******************************************************")
        print("i, j, lon, lat, lon2   : ",  i,   j,   lon,   lat, lon2 )
         
        file_path = Path("./ExtDECaps/Ext{0:.2f}_{1:.2f}.dat".format(float(lat),float(lon2)))
        if(file_path.is_file()):    print("The file exists.")
        elif(lat>=10.0 or lat<=-10.0):
            print("The file does not exist.")
            flag=1
            for k in range(100):
                dis=float(0.25+k*0.5)*1000.0## parsec        
                coords=SkyCoord(lon2*units.deg, lat*units.deg, distance=dis*units.pc, frame='galactic')
                E= decaps(coords)*(3.32/8.35816659571)## extinction in ks-band 
                if(math.isfinite(E)): 
                    tt[k,0], tt[k,1]= float(dis/1000.0),  float(E*1.0)## Ks_Extinction
                if(not math.isfinite(E) and k>1):
                    dist=float(dis/1000.0);#kpc  
                    if(dist<15.0 and k<30): m=abs(tt[k-1,1] - tt[k-2,1])/(0.5)*abs(30.0-k)/30.0
                    else:                   m=0.0      
                    tt[k,0], tt[k,1]=dist ,  float( m*0.5+tt[k-1,1] )         
                if(not math.isfinite(E) and k<2):
                    k=100
                    flag=0
                    print("Error no file, ", E,  k,  lon,  lon2,  lat)
            if(flag>0):   
                filew=open("./ExtDECaps/Ext{0:.2f}_{1:.2f}.dat".format(float(lat),float(lon2)),"w")           
                for k in range(100):
                    np.savetxt(filew,tt[k,:].reshape((-1,2)),fmt="%.2f  %.4f")
                filew.close()  
        
#####################################################################################################################      
      
      
      
      
      
      
      
      
      
      
        
