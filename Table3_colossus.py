import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import pylab
import Image
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
###=============================================================================

num0=[r"$b(degree)$",r"$l(degree)$",r"$filter$",r"$\log_{10}[n_{\star,simulated}]$",r"$\epsilon_{\star,1filter}$",r"$\epsilon_{LSST}$",r"$\epsilon_{lensed}$",r"$<m_{\star}>(mag)$",r"$<m_{base,d}>(mag)$",r"$<m_{base}>(mag)$",r"$<A>(mag)$",r"$<f_{b}>$",r"$\delta_{m}(mag)$",r"$\log_{10}[N_{PSF}]$",r"$No. line$",r"$counter$",r"$T_{obs}$",r"$cadence$",r"$epoches$"]##19
num2=[r"$~,u-band$",r"$~,g-band$",r"$~,r-band$",r"$~,i-band$",r"$~,z-band$",r"$~,y-band$"]


num3=[r"$b(degree)$",r"$l(degree)$",r"$\log_{10}[n_{\star, simulated}]$",r"$\epsilon_{\star,1filter}$",r"$\log_{10}[\widetilde{\varepsilon}_\mathrm{LSST}]$",r"$\log_{10}[<t_{\mathrm{E}}>(day)]$",r"$<R_{E}>(A.U.)$",r"$<D_{s}>(kpc)$",r"$<D_{l}>(kpc)$",r"$<V_{t}>(km~s^{-1})$",r"$<v_{l}>(km~s^{-1})$",r"$<v_{s}>(km~s^{-1})$",r"$<M_{l}>(M_{\cdot})$",r"$<u_{0}>$",r"$\log_{10}[N_{\star}(deg^{-2})]$",r"$\log_{10}[M_{t,l}(M_{\odot})]$",r"$\log_{10}[\epsilon(t_{E})/t_{E}]$",r"$\log_{10}[\Gamma\rm{(star^{-1}~yr^{-1})}]$",r"$\log_{10}[\tau\times~1e6]$",r"$\log_{10}[N_{e}\rm{(deg^{-2}~yr^{-1})}]$",r"$N_{event,cum}$",r"$\log_{10}[N_{\star,\rm{LSST}}\rm{(deg^{-2})}]$","$<V-I>(mag)$","$<f_{RGB}>$",r"$\log_{10}[\tau\times~1e6]$","$<I>(mag)$",r"$No. line$",r"$No. fields$",r"$T_{obs}(year)$",r"$cadence(day)$",r"$epoches$"]##31 

Nf=int(404)#int(449)
fi=int(5)

num1=int(Nf*fi)
dat1=np.zeros((num1,19))
dat1=np.loadtxt('./files/Tab3c_colossusA.txt') 

num2=int(Nf)
dat2=np.zeros((num2,31))
dat2=np.loadtxt('./files/Tab3c_colossusB.txt') 


ave1=np.zeros((19,fi,2))
ave2=np.zeros((31,2))
number=np.zeros((2))
nevent_com=np.zeros((2))
count=int(0)

for i in range(Nf):
    if(float(dat2[i,30])<200.0):# epoch=200.0, cadence =20 days, Tobs=10 years
        nevent_com[0]+=float(dat2[i,19])
        #print "N_event: ",  dat2[i,19],  "l(degree): ",  dat2[i,1], "b(degree): ",  dat2[i,0]
        number[0] += 1.0; 
        ave2[:,0] += dat2[i,:]  
        
        for j in range(fi):
            ave1[:,j,0]+=dat1[count,:]
            count+=1   
    else:                    ### epoch =900.0,  cadence=3 days , Tobs=10 years
        nevent_com[1]+= float(dat2[i,19])
        number[1]+= 1.0; 
        ave2[:,1] += dat2[i,:]
        for j in range(fi):
            ave1[:,j,1] +=dat1[count,:]
            count+=1   

##############################################
print "count: ", count 
print "No. fields with <cadence>~ 20 days(I): ", number[0], "<cadence>: ", ave2[29,0]/number[0] 
print "No. fields with <cadence>~ 3 days(II): ", number[1], "<cadence>: ", ave2[29,1]/number[1]

print "No. comulative events (I) ",  nevent_com[0],  "No. colulative events (II)", nevent_com[1]

fil=open("Table3_colossus.dat","w")
fil.close()
sav=np.zeros((2))


for i in range(31):
    fil=open("Table3_colossus.dat","a+")
    sav[0]=float(ave2[i,0]/number[0]/1.00)
    sav[1]=float(ave2[i,1]/number[1]/1.00) 
    fil.write(str(num3[i])+ "  *********************  ")
    np.savetxt(fil,sav.reshape((1,2)),fmt="%20.10f   %20.10f") 
    fil.close()
   
   
sav1=np.zeros((fi,2))
for i in range(19):
    fil=open("Table3_colossus.dat","a+")
    sav1[:,0]=ave1[i,:,0]/number[0]
    sav1[:,1]=ave1[i,:,1]/number[1]
    fil.write(str(num0[i])+ "\n")
    np.savetxt(fil,sav1[:,0].reshape((1,fi)),fmt="%.7f   %.7f   %.7f    %.7f    %.7f") 
    np.savetxt(fil,sav1[:,1].reshape((1,fi)),fmt="%.7f   %.7f   %.7f    %.7f    %.7f") 
    fil.write("************************************************************************************\n")
    fil.close()    
    

fil=open("Table3_colossus.dat","a+")
fil.write("No.  evenets  comulative:  "+ "\n")
np.savetxt(fil,nevent_com[:].reshape((1,2)),fmt="%.7f         %.7f") 
fil.close()    
print "**************************** END ********************************"  

