import matplotlib.pyplot as plt
import numpy as np
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

# s.lat,s.lon,j,d.nstar,nbb,d.nsdet[j],d.nevent,d.ave_aps[j],d.ave_apb[j],d.ave_apbd[j],d.ave_ex[j],
# d.ave_col[j],d.ave_bl[j],d.ave_dmag[j],d.ave_npsf[j] ##15

nf=int(6)


num1=[r"$b(degree)$",r"$l(degree)$",r"$filter$",r"$\log_{10}[n_{\star,simulated}]$",r"$\epsilon_{\star,1filter}$",r"$\epsilon_{LSST}$",r"$\epsilon_{lensed}$",r"$<m_{\star}>(mag)$",r"$<m_{bg,d}>(mag)$",r"$<m_{bg}>(mag)$",r"$<A>(mag)$",r"$<f_{b}>$",r"$\delta_{m}(mag)$",r"$\log_{10}[N_{PSF}]$"]



num2=[r"$~,u-band$",r"$~,g-band$",r"$~,r-band$",r"$~,i-band$",r"$~,z-band$",r"$y-band$"]

num=int(801*nf*81)
data=np.zeros((num,14))
data=np.loadtxt('./files/res060897A.txt') 
nrows, ncols=81, 801

for i in range(num):
    if(data[i,13]<1.0 ):
        data[i,13]=1.0; 


fil=open("Table7.dat","w")
fil.close()
sav1=np.zeros((num,2));  sav2=np.chararray((num)); sav3=np.zeros((num,4))
sav1[:,0]=data[:,0];
sav3[:,0]=data[:,8]; 
sav3[:,1]=data[:,10];  
sav3[:,2]=data[:,11];  


sav3[:,3]=np.log10(data[:,13]);  

for i in range(num): 
    if(data[i,1]>=260.0): 
        sav1[i,1]=float(data[i,1])-360.0
    else:
        sav1[i,1]=float(data[i,1])
    if(int(data[i,2])==0):
        sav2[i]= 'u'
    elif(int(data[i,2])==1):
        sav2[i]= 'g'
    elif(int(data[i,2])==2):
        sav2[i]= 'r'
    elif(int(data[i,2])==3):
        sav2[i]= 'i'
    elif(int(data[i,2])==4):
        sav2[i]='z'
    elif(int(data[i,2])==5):
        sav2[i]='y'

print sav2



for i in range(num):
    fil=open("Table7.dat","a+")
    ab=np.zeros(1,dtype=[('v1', float),('v2',float),('v3', 'U6'),('v4',float),('v5',float),('v6',float),('v7',float)] )
    ab['v1']=sav1[i,0];     ab['v2']=sav1[i,1]
    ab['v3']=sav2[i]
    ab['v4']=sav3[i,0]; ab['v5']=sav3[i,1]; ab['v6']=sav3[i,2];   ab['v7']=sav3[i,3]
    np.savetxt(fil,ab,fmt="%.2f   %.2f    %1s    %.5f    %.5f    %.5f    %.5f") 
    fil.close()

print "************************* END ***************************"


nu=int(0)
mapp=np.zeros((15,nf,nrows,ncols))
for i in range(nrows):
    for j in range(ncols):
        for k in range(nf): 
            if(data[nu,1]>=260.0): 
                data[nu,1]=float(data[nu,1])-360.0
            for l in range(14):
                if(l==3 or l==13):
                    if(float(data[nu,l])<0.0 or float(data[nu,l])==0.0):
                        print "ERROR data", float(data[nu,l]), "row: ", nu, "column: ", l
                        print "filter: ", k
                        print "*****************************************"
                        data[nu,l]= float(data[nu-nf,l])   
                        ##input("Enter a number ")
                    mapp[l,k,i,j]=np.log10(float(data[nu,l]))
                elif(l==4 or l==5 or l==6): 
                    mapp[l,k,i,j]=float(float(data[nu,l])/float(data[nu,3]+0.0000047652))
                else:
                    mapp[l,k,i,j]=float(data[nu,l])
            if(data[nu,0]!=(-10.0+ i*0.25) or data[nu,1]!=(-100.0+ j*0.25)):
                print "ERROR data", data[nu,:], "nu: ", nu
                print "first: ", -10.0+i*0.25, "second: ", -100.0+j*0.25
                input ("Enter a number")
            nu+=1
print "number of rows: ", nu
################################################################################        
tt=int(8)


for i in range(14):
    print "*************************************************************" 
    print "map is plotted  i ", i, "the : ",  str(num1[i])
    if (i==0 or i==1 or i==2 or i==3 or i==4 or i==6): # no dependence on filter
        v=np.zeros((tt))
        fig1=Figure()
        plt.figure(figsize=(10,2.5),facecolor='g')
        plt.imshow(mapp[i,0,:,:],cmap='viridis',extent=(-100.0,100.0,10.0,-10.0),interpolation='nearest',aspect='auto')
        plt.clim()
        plt.title(str(num1[i]),fontsize=11)
        minn=np.min(mapp[i,0,:,:])
        maxx=np.max(mapp[i,0,:,:])
        step=float((maxx-minn)/(tt-1.0));
        for m in range(tt):
            v[m]=round(float(minn+m*step),1)
        plt.colorbar(orientation='horizontal',shrink=1.0,pad=0.2,ticks=v)
        plt.clim(v[0]-0.005*step,v[tt-1]+0.005*step)
        plt.xticks(np.arange(-100.001,100.0,25.0))
        plt.ylim(-10.0,10.0)
        plt.xticks(fontsize=10, rotation=0)
        plt.yticks(fontsize=10, rotation=0)
        ax=plt.gca()
        ax.tick_params(direction='out',pad=5,top=False, right=False)
        plt.gca().invert_xaxis()
        fig1=plt.gcf()
        plt.xlabel('l(deg)',fontsize=10,labelpad=0.1)
        plt.ylabel('b(deg)',fontsize=10,labelpad=0.1)
        #fig1.savefig("./figs/mmap{0:d}.png".format(i),dpi=200)
        #Image.open("./figs/mmap{0:d}.png".format(i)).save("./figs/mmap{0:d}.jpg".format(i),'JPEG')
        fig1.savefig("./figs/eps/mmap{0:d}.eps".format(i), format='eps', dpi=200)
    else:
        for j in range(nf):
            v=np.zeros((tt))
            fig1=Figure()
            plt.figure(figsize=(10,2.5),facecolor='g')
            plt.imshow(mapp[i,j,:,:],cmap='viridis',extent=(-100.0,100.0,10.0,-10.0),interpolation='nearest',aspect='auto')
            plt.clim()
            plt.title(str(num1[i])+ str(num2[j]),fontsize=11)
            minn=np.min(mapp[i,j,:,:])
            maxx=np.max(mapp[i,j,:,:])
            if(i==11):
                maxx=float(1.0)
            if(j==0):
                print "for the filter u-band: minn: ", minn, "Max: ", maxx 
            step=float((maxx-minn)/(tt-1.0));
            for m in range(tt):
                v[m]=round(float(minn+m*step),1)
            plt.colorbar(orientation='horizontal',shrink=1.0,pad=0.2,ticks=v)
            plt.clim(v[0]-0.005*step,v[tt-1]+0.005*step)
            plt.xticks(np.arange(-100.001,100.0,25.0))
            plt.ylim(-10.0,10.0)
            plt.xticks(fontsize=10, rotation=0)
            plt.yticks(fontsize=10, rotation=0)
            ax=plt.gca()
            ax.tick_params(direction='out',pad=5,top=False, right=False)
            plt.gca().invert_xaxis()
            fig1=plt.gcf()
            plt.xlabel('l(deg)',fontsize=10,labelpad=0.1)
            plt.ylabel('b(deg)',fontsize=10,labelpad=0.1)
            #fig1.savefig("./figs/map{0:d}_{1:d}.png".format(i,j),dpi=200)
            #Image.open("./figs/map{0:d}_{1:d}.png".format(i,j)).save("./figs/map{0:d}_{1:d}.jpg".format(i,j),'JPEG')
            fig1.savefig("./figs/eps/map{0:d}_{1:d}.eps".format(i,j), format='eps', dpi=200)



################################################################################
