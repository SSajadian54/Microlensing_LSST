from __future__ import print_function
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarQuery
import dustmaps.bayestar
dustmaps.bayestar.fetch()
import numpy as np 
import astropy.units as units
#bayestar = BayestarQuery(version='bayestar2019')##max_samples=2, 
import matplotlib.pyplot as plt
import pylab as py
from matplotlib import rcParams
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
import math
import os 
#####################################################################################################################

fil2=open("Noext.dat", "w")
fil2.close()


save=0; 
for i1 in range(801): 
    i=i1+0
    for j in range(161): 
        lon= float(-100.0+0.25*i)
        lat=  float(-20.0+0.25*j)
        print ("*****************************************************")
        print ("i, j, lon, lat: ",  i, j, lon, lat )
        
        if(lon<0.0): lon2=float(360.0+lon)## that is between [0, 360] degree
        else:        lon2=float(lon)
        tt=np.zeros((100,2))
        tt= np.loadtxt("./Ext/Ext{0:.2f}_{1:.2f}.dat".format(float(lon2),float(lat)) )
        hh=0
        for k in range(100): 
            if(math.isfinite(tt[k,1])):  hh+=1
        if(hh<100):   
            fil2=open("Noext.dat", "a+")
            arr=np.array([ i, j, lon, lat])   
            np.savetxt(fil2, arr.reshape(-1,4),fmt="%d   %d    %.2f    %.2f")
            fil2.close()
            os.remove("./Ext/Ext{0:.2f}_{1:.2f}.dat".format(float(lon2),float(lat)))
        if(hh<100 and hh>0): 
            print(lon, lat, i, j, hh, lon2)
            input("Enter a number ")


tt=np.zeros((61439,4))        
tt= np.loadtxt("Noext.dat")        

plt.cla()
plt.clf()
plt.figure(figsize=(10, 2.5))
plt.plot(tt[:,2], tt[:,3],'go', markersize=1.1, alpha=0.85)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xlim([-100.0,100.0])
plt.ylim([-20.0,20.0])
plt.title(r"$\rm{No}~\rm{Extinction}~\rm{dustmaps.bayestar()}$", fontsize=14)
plt.xlabel(r"$Galactic~longitude(deg)$", fontsize=14)
plt.ylabel(r"$Galactic~latitude(deg)$",  fontsize=14)
fig=plt.gcf()
fig.tight_layout()
fig.savefig("Noextinction.jpg",dpi=200)        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
