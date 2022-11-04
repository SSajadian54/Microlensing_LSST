import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib as mpl
import matplotlib.cm as cm
import pylab
#import Image
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
###=============================================================================
#s.lat,s.lon,d.nstar,nbb,d.nevent,d.ave_te,d.ave_re,d.ave_Ds,d.ave_Dl,d.ave_vt,d.ave_vl,    
#d.ave_vs,d.ave_ml,d.ave_u0,s.Nstart,s.Rostart,d.ave_tei,d.ratio,d.ave_opt,d.nfin,Nmicro,s.nstart  #22



num1=[r"$b(degree)$",r"$l(degree)$",r"$\log_{10}[n_{\star, simulated}]$",r"$\epsilon_{\star,1filter}$",r"$\log_{10}[\widetilde{\varepsilon}_\mathrm{LSST}]$",r"$<\log_{10}[t_{\mathrm{E}}(day)]>$",r"$<R_{E}>(A.U.)$",r"$<D_{s}>(kpc)$",r"$<D_{l}>(kpc)$",r"$<V_{t}>(km~s^{-1})$",r"$<v_{l}>(km~s^{-1})$",r"$<v_{s}>(km~s^{-1})$",r"$<M_{l}>(M_{\cdot})$",r"$<u_{0}>$",r"$\log_{10}[N_{\star}(deg^{-2})]$",r"$\log_{10}[M_{t,l}(M_{\odot})]$",r"$\log_{10}[\epsilon(t_{E})/t_{E}]$",r"$\log_{10}[\Gamma_{obs}\rm{(star^{-1}~yr^{-1})}]$",r"$\log_{10}[\tau\times~1e6]$",r"$\log_{10}[N_{e}\rm{(deg^{-2}~yr^{-1})}]$",r"$N_{event,cum}$",r"$\log_{10}[N_{\star,\rm{LSST}}\rm{(deg^{-2})}]$","$<V-I>(mag)$","$<f_{RGB}>$",r"$\log_{10}[\tau_{DIA}\times~1e6]$","$<I>(mag)$",r"$\log_{10}[N_{\star,\rm{OGLE}}]$"] ##27

##print num1



NN=int(27)
num=int(801*81)
data=np.zeros((num,NN))
data=np.loadtxt('./files/res060897B.txt') 
nrows, ncols=81, 801


fil=open("Table8.dat","w");  #fih=open("Table6_paper.dat","w")
fil.close();  #fih.close()


for i in range(num):
    if (float(data[i,15])==0.0 or float(data[i,15])<0.0 or float(data[i,17])==0.0  or float(data[i,17])<0.0 or float(data[i,19])==0.0 or float(data[i,19])<0.0 or float(data[i,21])==0.0 or float(data[i,21])<0.0 or  float(data[i,24])==0.0 or float(data[i,24])<0.0): 
        print "BIG ERROR:  data: ", data[i,:]
        input("Enter a number")
    


sav=np.zeros((num,14))
sav[:,0]=data[:,0]; ## sav[:,1]=data[:,1]; 
sav[:,2]=np.log10(data[:,4]/(data[:,3]+0.0000047652)); 
sav[:,3]=data[:,5];###<log10(tE)>
sav[:,4]=data[:,6];#<RE>
sav[:,5]=data[:,7];  sav[:,6]=data[:,8];    sav[:,7]=data[:,9]; sav[:,8]=data[:,13]; 
sav[:,9]=np.log10(data[:,15]*0.25*0.25); 
sav[:,10]=np.log10(data[:,17]*10000000.0);  sav[:,11]=np.log10(data[:,19]/(10.0*0.25*0.25)); 
sav[:,12]=np.log10(data[:,21]/(0.25*0.25)); sav[:,13]=np.log10(data[:,24]); 
for i in range(num): 
    if(data[i,1]>=260.0): 
        sav[i,1]=float(data[i,1])-360.0
    else:
        sav[i,1]=float(data[i,1])
for i in range(num):
    fil=open("Table8.dat","a+")
   # fih=open("Table6_paper.dat","a+")
    np.savetxt(fil,sav[i,:].reshape(1,14),fmt=' %.2f  %.2f  %.5f  %.5f  %.5f  %.5f  %.5f  %.4f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f ') 
   # np.savetxt(fih,sav[i,:].reshape(1,14),delimiter="$ & $",fmt='%8.4f') 
        

      
    fil.close();  #fih.close()
### b(deg),  l(deg), \widetilda{\varepsilon}_{LSST}, log10(tE(day)),RE (A.U.), Ds(Kpc),Dl(kpc),V_{t}(km/s), u0, log10(M_{t,l}), 
#### log10(Gama*1.0e7), log10(Ne[deg^-2 year^-1]), log10(N_{\star,LSST}(deg^-2)), log10(\tau 10^6)
print "********************** END ********************"






data[:,19]=data[:,19]/(10.0*0.25*0.25)
data[:,21]=data[:,21]/(0.25*0.25)
#data[:,15]=data[:,15]/(0.25*0.25)

nu=int(0)
mapp=np.zeros((NN,nrows,ncols))
for i in range(nrows):
    for j in range(ncols):
        if(data[nu,1]>=260.0): 
            data[nu,1]=float(data[nu,1])-360.0
        for l in range(NN):
            if(l==15):
                if(float(data[nu,l])<0.0 or float(data[nu,l])==0.0):
                    print "ERROR data", float(data[nu,l]), "row: ", nu, "column: ", l
                    input("Enter a number ")
                mapp[l,i,j]=np.log10(float(data[nu,l])*0.25*0.25)
            elif(l==2 or l==14 or l==17 or l==18 or l==19 or l==21 or l==24 or l==16 or l==26):
                if(float(data[nu,l])<0.0 or float(data[nu,l])==0.0):
                    #print "ERROR data", float(data[nu,l]), "row: ", nu, "column: ", l, "data: ", data[nu,l]
                    data[nu,l]=1.0e-10;  
#                    input("Enter a number ")
                mapp[l,i,j]=np.log10(float(data[nu,l]))
            elif(l==3): 
                mapp[l,i,j]=float(float(data[nu,l])/float(data[nu,2]+0.0000047652))
            elif(l==4): 
                mapp[l,i,j]=float(np.log10(float(data[nu,l])/float(data[nu,3]+0.0000047652)))
            else:
                mapp[l,i,j]=float(data[nu,l])
        if(data[nu,0]!=(-10.0+ i*0.25) or data[nu,1]!=(-100.0+ j*0.25)):
            print "ERROR data", data[nu,:], "nu: ", nu
            print "first: ", -10.0+i*0.25, "second: ", -100.0+j*0.25
            input ("Enter a number")
        nu+=1
print "number of rows: ", nu
################################################################################        
tt=int(9)



for i in range(NN):
    v=np.zeros((tt))
    fig1=Figure()
    plt.figure(figsize=(10,2.5),facecolor='g')
    plt.imshow(mapp[i,:,:],cmap='viridis',extent=(-100.0,100.0,10.0,-10.0),interpolation='nearest',aspect='auto')
    plt.clim()
    plt.title(str(num1[i]),fontsize=11)
    #plt.grid(True)
    #if (i==5):
    #    minn=float(math.log10(20.0))
    #    maxx=float(math.log10(120.0))
    if (i==13):
        minn=float(0.1)
        maxx=float(0.7)
    else:
        minn=np.min(mapp[i,:,:])
        maxx=np.max(mapp[i,:,:])
    step=float((maxx-minn)/(tt-1.0));
    for m in range(tt):
        v[m]=round(float(minn+m*step),1)
    plt.colorbar(orientation='horizontal',shrink=1.0,pad=0.2,ticks=v)
    plt.clim(v[0]-0.005*step,v[tt-1]+0.005*step)
#    plt.xlim(-100.0,100.0,step=25.0)
    plt.xticks(np.arange(-100.001,100.0,25.0))
    #plt.yticks(np.arange(-10.001,+10.0,5.0))
    plt.ylim(-10.0,10.0)
    plt.xticks(fontsize=10, rotation=0)
    plt.yticks(fontsize=10, rotation=0)
    ax=plt.gca()
    ax.tick_params(direction='out',pad=5,top=False, right=False)
    #ax.tick_params(labeltop=False, labelright=False)
    plt.gca().invert_xaxis()
    fig1=plt.gcf()
    plt.xlabel('l(deg)',fontsize=10,labelpad=0.1)
    plt.ylabel('b(deg)',fontsize=10,labelpad=0.1)
    #plt.show() 
    #fig1.savefig("./figs/Map{0:d}.png".format(i))
    fig1.savefig("./figs/eps/Map{0:d}.eps".format(i), format='eps', dpi=200)
    #Image.open("./figs/Map{0:d}.png".format(i)).save("./figs/Map{0:d}.jpg".format(i),'JPEG')
    print "map is plotted  i ", i
    print "***************************************************" 


################################################################################
