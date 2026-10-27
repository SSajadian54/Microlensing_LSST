import pandas as pd 
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u


'''
nr=np.zeros((801, 161))
fil1=open("./counterOpsim.txt","w")
fil1.close()
for i in range(801): 
    for j in range(161): 
        print("*************************************************************************")
        lon=float(-100.0+ 0.25*i)
        lat=float(-20.0 + 0.25*j)          
        fil1=open("./OPSIM/opm{0:.2f}_{1:.2f}.txt".format(lon, lat),"r")  
        nq=sum(1 for line in  fil1)##number of rows
        nw=nq 
        if(nq>2): 
            par =np.zeros((nq,5))
            par =np.loadtxt("./OPSIM/opm{0:.2f}_{1:.2f}.txt".format(lon, lat))
            par2=np.unique(par, axis=0)###remove duplicated rows
            if(len(par2[:,0])>1): 
                nw =int(len(par2[:,0]))
                fil1.close()
                k=np.zeros((nw))
                k=np.argsort(par2[:,2])## sorted by time
                fil2=open("./OPSIM/opm{0:.2f}_{1:.2f}.txt".format(lon, lat),"w")  
                fil2.close()
                for l in range(nw): 
                    fil2=open("./OPSIM/opm{0:.2f}_{1:.2f}.txt".format(lon, lat),"a+")  
                    np.savetxt(fil2, par2[int(k[l]),:].reshape(-1,5), fmt= "%.3f  %.3f  %.9f  %d   %.6f")
                    fil2.close()
        nr[i,j]=nw           
        print(i,     j,     lon,     lat,      nq,      nw)      
        ##########################################################################################    
        tst=np.array([ i, j, nw ])
        fil2=open("./counterOpsim.txt","a+")
        np.savetxt(fil2, tst.reshape(-1,3),fmt="%d  %d  %d")
        fil2.close()    
        
print ("maximum numer of observatios,  ", np.max(nr[:,:]) , np.min(nr[:,:]) )        
'''      

'''
This part is for latitude first and longitude then
fil1=open("./counterOpsimB.txt","w")
fil1.close()
nr=np.zeros((801,161))
for j in range(161): 
    for i in range(801): 
        print("*************************************************************************")
        lon=float(-100.0+ 0.25*i)
        lat=float(-20.0 + 0.25*j)          
        fil1=open("./OPSIM/opm{0:.2f}_{1:.2f}.txt".format(lon,lat),"r")  
        nq=sum(1 for line in  fil1)##number of rows
        tst=np.array([i, j, nq])
        fil1=open("./counterOpsimB.txt","a+")
        np.savetxt(fil1,tst.reshape(-1,3),fmt="%d  %d  %d")
        fil1.close()           
'''    

    
## that is for inclusing the cadence        
fil1=open("./counterOpsimC.txt","w")
fil1.close()

for i in range(801): 
    for j in range(161): 

        lon=float(-100.0+ 0.25*i)
        lat=float(-20.0 + 0.25*j)          
        fil1=open("./OPSIM/opm{0:.2f}_{1:.2f}.txt".format(lon,lat),"r")  
        nq=sum(1 for line in  fil1)##number of rows
        par =np.zeros((nq,5))      
        par =np.loadtxt("./OPSIM/opm{0:.2f}_{1:.2f}.txt".format(lon, lat))  
        if(nq>1):  cadence=np.mean(np.abs(np.diff(par[:,2],axis=0)))
        else:      cadence=30.0;         
        tst=np.array([i,  j,  nq,  cadence])
        fil1=open("./counterOpsimC.txt","a+")
        np.savetxt(fil1,tst.reshape(-1,4),fmt="%d  %d  %d  %.7f")
        fil1.close()            
print("*************************************************************************")          
          
          
          
                
