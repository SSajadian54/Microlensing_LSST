from __future__ import print_function
import numpy as np 
import astropy.units as units
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import math
import os 
from pathlib import Path

rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['axes.titlesize'] = 20
cmap=plt.get_cmap('viridis')

###################################################################################################



def tickfun( start, x, dd0):
    return( start+x*dd0)
    
    
###################################################################################################  
 

nx=801
ny=161
ddeg=0.25

mapp=np.zeros((ny, nx))
tt=np.zeros((100,2))
nv=9
v =np.zeros((nv))       



fil2=open("./NoExtinction.dat", "w")
fil2.close()


fil3=open("./saved_direction.txt", "w")
fil3.close()
    
no=0  
for l in range(100):         
    for i in range(nx): 
        for j in range(ny): 
            lon= float(-100.0+0.25*i)
            lat=  float(-20.0+0.25*j)
            if(lon<0.0): lon2=float(360.0+lon)
            else:        lon2=float(lon)
            #print("*****************************************************")
            #print("Ext{0:.2f}_{1:.2f}.dat".format(float(lat),float(lon2) )  )
            #print("i, j, lon, lat: ",  i,   j,   lat ,  lon2)
            file_path = Path("./ExtFinal/Ext{0:.2f}_{1:.2f}.dat".format(float(lat),float(lon2)))
            if(file_path.is_file()): 
                tt=np.loadtxt("./ExtFinal/Ext{0:.2f}_{1:.2f}.dat".format(float(lat),float(lon2)) )
                mapp[j,i]=tt[l,1] 
                
                if(mapp[j,i]>9.0 and abs(lon)>3.0): 
                    print(mapp[j,i], j, i, lon, lat, lon2, l)
                    input("Enter a number ")
                if(l==1):    
                    fil3=open("./saved_direction.txt", "a+")
                    arr=np.array([ lat, lon2 ])
                    np.savetxt(fil3, arr.reshape(-1,2), fmt="%.2f    %.2f")
                    fil3.close()
            else: 
                mapp[j,i]=0.0 
                if(l==1):    
                    fil2=open("./NoExtinction.dat", "a+")
                    arr=np.array([ i, j, lat, lon2])   
                    #fil2.write("Ext{0:.2f}_{1:.2f}.dat".format(float(lat),float(lon2)))
                    np.savetxt(fil2, arr.reshape(-1,4),fmt="%d   %d    %.2f    %.2f")
                    fil2.close() 
                    no+=1
    ###################################################################################################
    minn=float(np.min(mapp[:,:]))
    maxx=float(np.max(mapp[:,:]))
    stt =float((maxx-minn)/(nv-1.0))
    for i in range(nv):
        v[i]=round(float(minn+i*stt),2)
    print("minn, maxx, step:  ",   minn,   maxx,   stt,  v, )
    
    plt.cla()
    plt.clf()
    fig=plt.figure(figsize=(10,3.5),facecolor='w')
    ax=fig.add_subplot(111)
    plt.imshow(mapp,cmap=cmap,extent=(0.0,800.0,0.0,160.0),interpolation='nearest',aspect='auto', origin='lower')
    plt.clim()
    cbar=plt.colorbar(orientation='horizontal',shrink=0.8,pad=0.2,ticks=v)
    cbar.ax.tick_params(labelsize=16)
    plt.clim(v[0]-0.05*stt , v[nv-1]+0.05*stt)
    plt.xticks(fontsize=16, rotation=0)
    plt.yticks(fontsize=16, rotation=0)
    plt.title(r"$\rm{E(B-V)}(\rm{mag});~\rm{D}(\rm{~kpc})=$"+str(0.25+l*1*0.5),fontsize=16)
    plt.xlim(nx*0.0 , nx*1.0)
    plt.ylim(ny*0.0 , ny*1.0)
    ticc=np.array([nx*0.1, nx*0.3, nx*0.5, nx*0.7, nx*0.9 ])
    ax.set_xticks(ticc,labels=[int(j) for j in tickfun(-100.0, ticc,ddeg) ])
    ticc=np.array([ny*0.1, ny*0.3, ny*0.5, ny*0.7, ny*0.9 ])
    ax.set_yticks(ticc,labels=[int(j) for j in tickfun(-20.0,  ticc,ddeg) ])
    plt.xlabel(r"$\rm{l}(\rm{deg})$",fontsize=17,labelpad=0.05)
    plt.ylabel(r"$\rm{b}(\rm{deg})$",fontsize=17,labelpad=0.05)
    fig=plt.gcf()
    fig.tight_layout()
    fig.savefig("./figs/Ext{0:d}.jpg".format(l), dpi=200)
    
        
#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  
#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  
arr= np.zeros((no,4))
arr=np.loadtxt("./NoExtinction.dat")        
Mapp=np.zeros((ny,nx))
for i in range(no): 
    Mapp[int(arr[i,1]),int(arr[i,0])]=1.0
    
#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  
plt.cla()
plt.clf()
fig=plt.figure(figsize=(10,3.5),facecolor='w')
ax=fig.add_subplot(111)
plt.imshow(Mapp,cmap=cmap,extent=(0.0,800.0,0.0,160.0),interpolation='nearest',aspect='auto', origin='lower')
#plt.clim()
#cbar=plt.colorbar(orientation='horizontal',shrink=0.8,pad=0.2,ticks=v)
#cbar.ax.tick_params(labelsize=16)
#plt.clim(v[0]-0.05*stt , v[nv-1]+0.05*stt)
plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16, rotation=0)
plt.title(r"$No extinction$",fontsize=16)
plt.xlim(nx*0.0 , nx*1.0)
plt.ylim(ny*0.0 , ny*1.0)
ticc=np.array([nx*0.1, nx*0.3, nx*0.5, nx*0.7, nx*0.9 ])
ax.set_xticks(ticc,labels=[int(j) for j in tickfun(-100.0, ticc,ddeg) ])
ticc=np.array([ny*0.1, ny*0.3, ny*0.5, ny*0.7, ny*0.9 ])
ax.set_yticks(ticc,labels=[int(j) for j in tickfun(-20.0,  ticc,ddeg) ])
plt.xlabel(r"$\rm{l}(\rm{deg})$",fontsize=17,labelpad=0.05)
plt.ylabel(r"$\rm{b}(\rm{deg})$",fontsize=17,labelpad=0.05)
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./figs/Noextinction.jpg", dpi=200)
       
#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  
flag=0;        
for k in range(no): 
    i=int(arr[k,0])   
    j=int(arr[k,1])   
    lon0=  float(-100.0+0.25*i)
    lat0=  float(-20.0+0.25*j)
    if(lon0<0.0): lon20=float(360.0+lon0)
    else:         lon20=float(lon0)
    file_path=Path("./ExtFinal/Ext{0:.2f}_{1:.2f}.dat".format(float(lat0),float(lon20)))
    if(file_path.is_file()): flag=1 
    else:                    
        flag=0
        for l in range(12): 
            if(l==0 and flag==0):   lat, lon2= float(lat0+0.25), float(lon20+0.00)
            if(l==1 and flag==0):   lat, lon2= float(lat0-0.25), float(lon20+0.00)
            if(l==2 and flag==0):   lat, lon2= float(lat0+0.00), float(lon20+0.25)
            if(l==3 and flag==0):   lat, lon2= float(lat0+0.00), float(lon20-0.25)
            if(l==4 and flag==0):   lat, lon2= float(lat0+0.50), float(lon20+0.00)
            if(l==5 and flag==0):   lat, lon2= float(lat0-0.50), float(lon20+0.00)
            if(l==6 and flag==0):   lat, lon2= float(lat0+0.00), float(lon20+0.50)
            if(l==7 and flag==0):   lat, lon2= float(lat0+0.00), float(lon20-0.50)
            if(l==8 and flag==0):   lat, lon2= float(lat0+0.25), float(lon20+0.25)
            if(l==9 and flag==0):   lat, lon2= float(lat0-0.25), float(lon20-0.25)
            if(l==10 and flag==0):  lat, lon2= float(lat0+0.25), float(lon20-0.25)
            if(l==11 and flag==0):  lat, lon2= float(lat0-0.25), float(lon20+0.25)
            file_path2 = Path("./ExtFinal/Ext{0:.2f}_{1:.2f}.dat".format(float(lat),float(lon2)))
            if(file_path2.is_file()):    
                #print("The file exists.")
                tt=np.loadtxt("./ExtFinal/Ext{0:.2f}_{1:.2f}.dat".format(float(lat),float(lon2)))  
                filew=open("./ExtFinal/Ext{0:.2f}_{1:.2f}.dat".format(float(lat0),float(lon20)),"w")           
                for kj in range(100):
                    np.savetxt(filew,tt[kj,:].reshape((-1,2)),fmt="%.2f   %.4f")
                filew.close()  
                flag=1  
    if(flag==0): 
        print("No file:  ",   lat0,   lon20, k, flag)       
       
       
       
       
       
       
       
                
        
        
        
        
        
