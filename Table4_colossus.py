import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib as mpl
import matplotlib.cm as cm
import pylab
import Image
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import sys, csv ,operator
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
###=============================================================================

num1=[r"$b(degree)$",r"$l(degree)$",r"$\log_{10}[n_{\star, simulated}]$",r"$\epsilon_{\star,1filter}$",r"$\log_{10}[\widetilde{\varepsilon}_\mathrm{LSST}]$",r"$\log_{10}[<t_{\mathrm{E}}>(day)]$",r"$<R_{E}>(A.U.)$",r"$<D_{s}>(kpc)$",r"$<D_{l}>(kpc)$",r"$<V_{t}>(km~s^{-1})$",r"$<v_{l}>(km~s^{-1})$",r"$<v_{s}>(km~s^{-1})$",r"$<M_{l}>(M_{\cdot})$",r"$<u_{0}>$",r"$\log_{10}[N_{\star}(deg^{-2})]$",r"$\log_{10}[M_{t,l}(M_{\odot})]$",r"$\log_{10}[\epsilon(t_{E})/t_{E}]$",r"$\log_{10}[\Gamma\rm{(star^{-1}~yr^{-1})}]$",r"$\log_{10}[\tau\times~1e6]$",r"$\log_{10}[N_{e}\rm{(deg^{-2}~yr^{-1})}]$",r"$N_{event,cum}$",r"$\log_{10}[N_{\star,\rm{LSST}}\rm{(deg^{-2})}]$","$<V-I>(mag)$","$<f_{RGB}>$",r"$\log_{10}[\tau\times~1e6]$","$<I>(mag)$",r"$No. line$",r"$No. Field$",r"$T_{obs}(yr)$",r"$<cadence>(days)$",r"$epoch$"] 



fi=int(5)
fil1=open("tab4A.dat","w")
fil2=open("tab4B.dat","w")
fil1.close();  fil2.close()


NN=int(31)
num=int(405)
data1=np.zeros((num*2,NN))
data1=np.loadtxt('./files/Tab4c_colossusB_old.txt') 
data2=np.zeros((num*fi*2,19))
data2=np.loadtxt('./files/Tab4c_colossusA_old.txt')
ne_sum=np.zeros((2))

sav=np.zeros((2,num,12))
nara=int(0)
for i in range(num*2):
    j= int(i*fi+2)
    l=int(i/2)
    k=int(i%2)
################################################################################
    sav[k,l,0]=data1[i,27];###No. field
    sav[k,l,1]=data1[i,1]### longtitute
    sav[k,l,2]=data1[i,0];###latitute
    sav[k,l,3]=data1[i,5]### log(tE)
    sav[k,l,4]=data1[i,13];### u0
    sav[k,l,5]=data2[j,9]### m_bg_detected
    sav[k,l,6]=data2[j,11]### blending
    sav[k,l,7]=np.log10(data1[i,17]*10000000.0)#### log(Gama*10^7)
    sav[k,l,8]=float(data1[i,19]/(data1[i,26]*0.25*0.25))### N_e / deg^2
    sav[k,l,9]= np.log10(data1[i,4]/(data1[i,3]+0.0000047652)); ### epsilon
    sav[k,l,10]=float(data1[i,26]*0.25*0.25);### area
    sav[k,l,11]=k; ###float(data1[i,27]*0.25*0.25);### area
    if(k==0):   
        ne_sum[0] += float(sav[0,l,8])
        fil1=open("tab4A.dat","a+")
        np.savetxt(fil1,sav[0,l,:].reshape(1,12),fmt='%d   %.4f   %.4f  %.5f   %.5f  %.5f   %.5f   %.5f   %.5f  %.5f  %.5f  %d') 
        fil1.close()
        if(data1[i,30]!=160 or data1[i,29]!=21.55999 or data1[i,28]!=10.0):
            print "ERROR please check,   k", k,  "  epoch: ",  data1[i,30], "  cadence: ", data1[i,29], "  Tobs: ", data1[i,28] 
            input("Enter a number")
    elif(k==1): 
        ne_sum[1] += float(sav[1,l,8])  
        fil2=open("tab4B.dat","a+")
        np.savetxt(fil2,sav[1,l,:].reshape(1,12),fmt='%d   %.4f   %.4f  %.5f   %.5f  %.5f   %.5f   %.5f   %.5f  %.5f  %.5f  %d') 
        fil2.close() 
        if(data1[i,30]!=900 or data1[i,29]!=4.00989 or data1[i,28]!=10.0):
            print "ERROR please check,  k", k,  "  epoch: ",  data1[i,30], "  cadence: ", data1[i,29], "  Tobs: ", data1[i,28] 
            input("Enter a number")
    else:
        print "ERROR k: ", k
        input("Enter a number")
print "********************** END OF READIGN *****************************"        
print "Ne_sum[0]:   ",  ne_sum[0],  "  Ne_sum[1]:    ",  ne_sum[1]                        
            


################################################################################
ave_rate=float(0.0);   avegama=float(0.0)
nte=int(0)
nmag=int(0)
data=np.zeros((num,12));       data2=np.zeros((num,12)); 
data=np.loadtxt('tab4A.dat');  data2=np.loadtxt('tab4B.dat');  
counter=np.zeros((num))
counter=np.argsort(data[:,0])
fil1=open("tab4A2.dat","w");  fil2=open("tab4B2.dat","w")
fil1.close();    fil2.close()
fil4=open("Table6.dat","w");   fil5=open("Table6_colossusb.dat","w")
fil4.close();  fil5.close()
for i in range(num):
    k=int(counter[i])
    fil1=open("tab4A2.dat","a+") ;      fil2=open("tab4B2.dat","a+");  fil4=open("Table6.dat","a+"); fil5=open("Table6_colossusb.dat","a+")
    np.savetxt(fil1,data[k,:].reshape(1,12),fmt='%d   %.4f   %.4f    %.5f   %.5f  %.5f   %.5f   %.5f   %.5f    %.5f    %.5f  %d') 
    np.savetxt(fil2,data2[k,:].reshape(1,12),fmt='%d   %.4f   %.4f    %.5f   %.5f  %.5f   %.5f   %.5f   %.5f    %.5f    %.5f  %d')
    last=np.zeros((15))
    last[0]=int(data[k,0]);  last[1]=float(data[k,1]);  last[2]=float(data[k,2]); 
    last[3]=float(data[k,3]); last[4]=float(data2[k,3])### log10(tE)
    last[5]=float(data[k,5]); last[6]=float(data2[k,5])### mbg_detected
    last[7]=float(data[k,6]); last[8]=float(data2[k,6])
    last[9]=float(data[k,7]);last[10]=float(data2[k,7])
    last[11]=float(data[k,8]); last[12]=float(data2[k,8])### Ne(deg^-2)
    last[13]=float(data[k,9]); last[14]=float(data2[k,9])### epsilon
    avegama  +=pow(10.0,float(last[10]-last[9]))
    ave_rate +=pow(10.0,float(last[14]-last[13]))
    np.savetxt(fil4,last[1:].reshape(1,14),fmt='%6.2f,%6.2f    %4.2f,%4.2f    %5.2f,%5.2f     %4.2f,%4.2f    %4.2f,%4.2f   %6.2f,%6.2f    %5.2f,%5.2f')
    np.savetxt(fil5,last.reshape(1,15),fmt='%d  %6.2f %6.2f    %4.2f %4.2f    %5.2f %5.2f     %4.2f %4.2f    %4.2f %4.2f   %6.2f %6.2f    %5.2f %5.2f')
    fil4.close();  fil5.close()
    fil2.close();  fil1.close()


    if(float(last[4])>float(last[3])):
        nte+=1


    if(float(last[12])<float(last[11]) or float(last[14])<float(last[13]) ):
        print "ERROR row: ", i, "Ne(deg^-2)[I]: ",  last[11], "Ne(deg^-2)[II]: ", last[12], "b(deg): ",  last[2]
        print "ERROR row: ", i, "epcilon[I]: ",  last[13], "epcilon[II]: ", last[44], "b(deg): ",  last[2]
        input("Enter a number: ")

    if(float(last[6])<float(last[5])):
        print "ERROR row: ", i, "m_bg(mag)[I]: ",  last[5], "m_bg(mag)[II]: ", last[6], "b(deg): ",  last[2]
        nmag+=1


print "No.  of events tE is longer with startegy I: ", nte
print "No.  of events Mag_base is fainter for strategy I  than for strategy II: ", nmag
print "<epsilon[I]/epsilon[II]>: ", float(ave_rate/num)
print "<Gama[II]/Gama[I]>: ", float(avegama/num)
print "********************** END OF SORTING 1 *****************************"        
################################################################################



