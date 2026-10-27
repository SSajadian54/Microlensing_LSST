#27/06/1405 Friday 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib.figure import Figure
from numpy import conj
import pylab
from matplotlib import cm
from matplotlib import colors
cm=colors.ListedColormap(['purple', 'blue', 'darkgreen','yellowgreen', 'orange', 'red'])
from matplotlib import gridspec
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
cmap=plt.get_cmap('viridis')
mpl.rc('font',**{'family':'serif','serif':['Palatino']})
mpl.rc('text', usetex=True)
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams["font.size"] = 11.5
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Computer Modern Sans"]
mpl.rcParams["text.usetex"] = True
mpl.rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
mpl.rcParams['text.usetex'] = True


yr=365.2425

f2=open("./Test.dat","r")
nm=sum(1 for line in f2)  
mag=np.zeros((nm,3))  
mag= np.loadtxt("./Test.dat")
    
    
plt.clf()
plt.cla()
fig=plt.figure(figsize=(8,6))
plt.plot(mag[:,0],mag[:,1],"k-", lw=2.10, label=r"$\rm{Lensing}$")
plt.plot(mag[:,0],mag[:,2],"b--",lw=2.10, label=r"$\rm{Lensing}+\rm{parallax}$")
plt.xlabel(r"$\rm{time(days)}$",fontsize=18)
plt.ylabel(r"$\rm{Magnification},~r-\rm{band}$", fontsize=18)
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.xlim(0.0,10.0*yr) 
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=False,prop={"size":13.5},framealpha=0.65)
fig=plt.gcf()
fig.savefig("./Lc.jpg",dpi=200) 
    
    
    
    
    
    
    
    
    
