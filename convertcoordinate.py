from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy  as np 
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib.figure import Figure
#from IPython.display import display, Math
import math
import pylab
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
from scipy.interpolate import make_interp_spline, BSpline


#################################################################################################################
fil1=open("convertcoordinateGB3.txt","w")
fil1.close()


for j in range(161): 
    for i in range(801): 
        lon=float(-100.0+ 0.25*i)
        lat=float( -20.0+ 0.25*j) 
        
        if(lon<0.0): lon2=float(lon+360.0)   
        else:        lon2=lon; 
        
        
        c=SkyCoord(l=float(lon2), b=float(lat), unit=(u.degree ,  u.degree), frame='galactic')
        ceq=c.transform_to('icrs')
        RA= float(ceq.ra.degree)
        DEC=float(ceq.dec.degree)
        
        arr=np.array([lon, lat, RA, DEC])
        fil1=open("convertcoordinateGB3.txt","a+")
        np.savetxt(fil1,arr.reshape(-1,4),fmt= "%.2f    %.2f    %.2f    %.2f")
        fil1.close()

####################################################################################################################
####################################################################################################################












