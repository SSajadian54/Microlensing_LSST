##Elahi be omid to 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib.figure import Figure
#from IPython.display import display, Math
import math
import pylab
from matplotlib import cm
from matplotlib import colors
from matplotlib import gridspec
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import pandas as pd 
import seaborn as sns
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
import warnings
warnings.filterwarnings("ignore")

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
Mjupiter=1.89813*pow(10.0,27.0);
Msun=1.98892*pow(10.,30);
eps =float(0.00000006567576675555995)
ddeg=float(0.25)
yr  =float(365.2425);
Tobs=float(10.0*yr)
AU  =float(1.4960*np.power(10.0,11.0))
tt  =int(9)
v   =np.zeros((tt))
na  =int(14)
area=ddeg*ddeg  
nx=801
ny=161

##############################################################################################################################
#############################################################################################################################

fil1=open("Table1.txt","w")
fil2=open("Table2.txt","w")
fil3=open("Table3.txt","w")
fil4=open("Table4.txt","w")
fil5=open("Table5.txt","w")
fil1.close();   
fil2.close();
fil3.close();   
fil4.close();
fil5.close();


##############################################################################################################################
####                                        Functions       #####
##############################################################################################################################
def funcx(xarr,val):
    dd=-1
    if(val<xarr[0]  or val==xarr[0]):        dd=0
    elif(val>xarr[na-1] or val==xarr[na-1]): dd=int(na-1)
    else: 
        for i in range(na-1):  
            if(float((val-xarr[i])*(val-xarr[i+1]))<=0.0 or val==xarr[i]): 
                dd=int(i)
                break        
    if(dd<0 or dd>=na or len(xarr)!=na or math.isnan(val) ): 
        print("Error funcx,  ", dd, val, xarr, na);  
        input("Enter a number ") 
    return(int(dd))
    
#################################################
    
def tickfun(x, start, dd0, rho):
    return((start+x*dd0)*rho)
#################################################     
def plotMap(Map, pointx, pointy, i, minn, maxx , title): 
    plt.cla()
    plt.clf()
    fig=plt.figure()
    #fig, ax=fig.add_subplot(figsize=(10,3.5))
    fig, ax=plt.subplots(figsize=(10, 3.5)) 
    plt.imshow(Map, extent=[0, nx, 0, ny], cmap=cmap, interpolation='nearest', aspect='equal', origin='lower')
    plt.plot(pointx, pointy, "ro",markersize=5.0)
    plt.clim()
    step=float((maxx-minn)/(tt-1.0));
    for m in range(tt):
        v[m]=round(float(minn+m*step),1)
    cbar=plt.colorbar(orientation='horizontal',shrink=0.95,pad=0.19,ticks=v)
    cbar.ax.tick_params(labelsize=17)
    plt.clim(v[0]-0.005*step,v[tt-1]+0.005*step)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.title(title, fontsize=18)
    plt.xlim(nx*0.0 , nx*1.0)
    plt.ylim(ny*0.0 , ny*1.0)
    ticx=np.array([nx*0.1, nx*0.3, nx*0.5, nx*0.7, nx*0.9 ])
    ticy=np.array([ny*0.1, ny*0.3, ny*0.5, ny*0.7, ny*0.9 ])
    ax.set_xticks(ticx,labels=[int(j) for j in tickfun(ticx,-100.0,ddeg, 1.0 ) ])
    ax.set_yticks(ticy,labels=[int(j) for j in tickfun(ticy, -20.0,ddeg, 1.0 ) ])
    plt.xlabel(r"$l(\rm{deg})$",fontsize=20,labelpad=0.05)
    plt.ylabel(r"$b(\rm{deg})$",fontsize=20,labelpad=0.05)
    fig=plt.gcf()
    fig.tight_layout()
    fig.savefig("./files/Hist/Map"+"{0:d}.jpg".format(i),dpi=200)     
        
        
        
##############################################################################################################################s        
###########################################################################################################################################
#                     MapLMC.dat   or MapSMC.dat
###########################################################################################################################################
##############################################################################################################################


titl=[r"$nl1$", r"$nl2$", r"$\log_{10}[t_{\rm{E}}(\rm{days})]$", r"$\log_{10}[t_{\rm{E}}(\rm{days})]$",    
r"$R_{\rm{E}}(\rm{AU})$", r"$R_{\rm{E}}(\rm{AU})$", r"$\log_{10}[\pi_{\rm{E}}]$", r"$\log_{10}[\pi_{\rm{E}}]$", 
r"$\log_{10}[\theta_{\rm{E}}(\rm{mas})]$", r"$\log_{10}[\theta_{\rm{E}}(\rm{mas})]$", 
r"$V_{\rm{rel}}(\rm{km}/s)$", r"$V_{\rm{rel}}(\rm{km}/s)$", 
r"$u_{0}$", r"$u_{0}$", r"$M_{\rm{l}}(M_{\odot})$", r"$M_{\rm{l}}(M_{\odot})$",
r"$\log_{10}[\tau_{\rm{b}}(\times 10^{6})]$", r"$\log_{10}[\tau_{\rm{b}}(\times 10^{6})]$", 
r"$D_{\rm{l}}(\rm{kpc})$", r"$D_{\rm{l}}(\rm{kpc})$", r"$D_{\star}(\rm{kpc})$", r"$D_{\star}(\rm{kpc})$", 
r"$v_{\rm{l}}(\rm{km}/s)$", r"$v_{\rm{l}}(\rm{km}/s)$", r"$v_{\star}(\rm{km}/s)$",  r"$v_{\star}(\rm{km}/s)$",
r"$m_{\rm{base}}(\rm{mag}),~K-\rm{band}$", r"$m_{\rm{base}}(\rm{mag}),~K-\rm{band}$", 
r"$f_{\rm{b}}, r-\rm{band}$", r"$f_{\rm{b}}, r-\rm{band}$", 
r"$\rm{FWHM}(\rm{days})$", r"$\rm{FWHM}(\rm{days})$",
r"$\mu_{\star,~n}$", r"$\mu_{\star,~n}$", 
r"$\Delta\theta_{0}/\sigma_{\rm{A}}$",  r"$\Delta\theta_{0}/\sigma_{\rm{A}}$",  
r"$f_{\rm{Halo}}[\%]$",  r"$f_{\rm{Halo}}[\%]$",
r"$\log_{10}[\mu_{\rm{rel}}(\rm{mas}/\rm{yrs})]$", r"$\log_{10}[\mu_{\rm{rel}}(\rm{mas}/\rm{yrs})]$", 
r"$m_{\star}(\rm{mag}), r-\rm{band}$",  r"$m_{\star}(\rm{mag}), r-\rm{band}$",  
r"$N_{\rm{blend}}, r-\rm{band}$", r"$N_{\rm{blend}}, r-\rm{band}$", 
r"$\rm{Ext}(\rm{mag}), r-\rm{band}$", r"$\rm{Ext}(\rm{mag}), r-\rm{band}$",
r"$\log_{10}[\rho_{\star}]$", r"$\log_{10}[\rho_{\star}]$", 
r"$\rm{Ext}(\rm{mag}), \rm{visible}$", r"$\rm{Ext}(\rm{mag}), \rm{visible}$",  
r"$\rm{CatalogLens}$", r"$\rm{CatalogLens}$", 
r"$\log_{10}[\tau_{\rm{p}}(\times 10^{6})]$", r"$\log_{10}[\tau_{\rm{p}}(\times 10^{6})]$",
r"$\varepsilon_{\rm{b}}[\%]$",r"$\varepsilon_{\rm{p}}[\%]$",r"$\varepsilon_{\rm{D}}[\%]$",
r"$\varepsilon_{\rm{L}}[\%]$", 
r"$\log_{10}[\left<\epsilon(t_{\rm{E}})/t_{\rm{E}}\right>(\rm{yrs}^{-1})]$",
r"$\log_{10}[\Gamma(\rm{star}^{-1}.\rm{yr}^{-1})\times 10^{7}]$",
r"$\log_{10}[dN_{\rm{e}}/d\Omega(\rm{deg}^{-2})]$",
r"$\log_{10}[dM/d\Omega(M_{\odot}.\rm{deg}^{-2})]$", r"$\log_{10}[dN_{\star}/d\Omega(\rm{deg}^{-2})]$",
r"$\log_{10}[dN_{\star, \rm{Rubin}}/d\Omega(\rm{deg}^{-2})]$",
r"$\log_{10}[\delta_{u_{0}}]$", r"$\log_{10}[\delta_{t_{\rm{E}}}]$", r"$\log_{10}[\delta_{f_{\rm{b}}}]$", 
r"$\log_{10}[\delta_{\pi_{\rm{E}}}]$", r"$\log_{10}[\delta_{\theta_{\rm{E}}}]$", 
r"$\log_{10}[\delta_{M_{\rm{l}}}]$", r"$\log_{10}[\delta_{D_{\rm{l}}}]$", 
r"$\log_{10}[\delta_{\mu_{\rm{l}}}]$", r"$\log_{10}[\delta_{\mu_{\star}}]$",r"$\log_{10}[\delta_{\rho_{\star}}]$",
r"$\log_{10}[\Delta\chi^{2}_{\rm{p}}/N_{\rm{data}}]$", r"$\log_{10}[\Delta\chi^{2}_{\rm{a}}/N_{\rm{data}}]$",
r"$f_{\pi}[\%]$", r"$f_{\delta\theta_{\rm{c}}}[\%]$", r"$f_{\Delta \theta}[\%]$", 
r"$N_{\rm{data}, \Delta\theta}$", r"$m_{\rm{Image}_{1}}(\rm{mag}), r-\rm{band}$", 
r"$m_{\rm{Image}_{2}}(\rm{mag}), r-\rm{band}$", r"$\Delta\theta/\sigma_{\rm{a},2}$", r"$No. \rm{visit}$", 
r"$\rm{N}_{\rm{simulated}}$",    r"$\rm{N}_{\rm{detected}}$", r"$\rm{N}_{\rm{lensing}}$", r"$N_{\rm{Error}}$",
r"$Right Ascension$", r"$Declination$", r"$f_{\rm{BD}}[\%]$", r"$f_{\rm{BD}}[\%]$", r"$f_{\rm{MS}}[\%]$", 
r"$f_{\rm{MS}}[\%]$", r"$f_{\rm{BH}}[\%]$", r"$f_{\rm{BH}}[\%]$", 
r"$u-r(\rm{mag})$", r"$u-r(\rm{mag})$", r"$g-i(\rm{mag})$", r"$g-i(\rm{mag})$", r"$r-z(\rm{mag})$", r"$r-z(\rm{mag})$", 
r"$m_{\rm{base}}(\rm{mag}), u-\rm{band}$", r"$m_{\rm{base}}(\rm{mag}), u-\rm{band}$", 
r"$m_{\rm{base}}(\rm{mag}), g-\rm{band}$", r"$m_{\rm{base}}(\rm{mag}), g-\rm{band}$", 
r"$m_{\rm{base}}(\rm{mag}), r-\rm{band}$", r"$m_{\rm{base}}(\rm{mag}), r-\rm{band}$", 
r"$m_{\rm{base}}(\rm{mag}), i-\rm{band}$", r"$m_{\rm{base}}(\rm{mag}), i-\rm{band}$", 
r"$m_{\rm{base}}(\rm{mag}), z-\rm{band}$", r"$m_{\rm{base}}(\rm{mag}), z-\rm{band}$", 
r"$m_{\rm{base}}(\rm{mag}), y-\rm{band}$", r"$m_{\rm{base}}(\rm{mag}), y-\rm{band}$", 
r"$f_{u_{0}}[\%]$", r"$f_{t_{\rm{E}}}[\%]$", r"$f_{f_{\rm{b}}}[\%]$", r"$f_{\pi_{\rm{E}}}[\%]$", 
r"$f_{\theta_{\rm{E}}}[\%]$",  r"$f_{M_{\rm{l}}}[\%]$",  r"$f_{D_{\rm{l}}}[\%]$",  r"$f_{\mu_{\rm{l}}}[\%]$", 
r"$f_{\mu_{\star}}[\%]$",  r"$f_{\rho_{\star}}[\%]$", r"$\rm{Filter}_{1}$", r"$\rm{Filter}_{2}$", 
r"$\rm{color}_{m}$", r"$\rm{color}_{m}$", r"$flag_{\rm{color}}$" ]



nb=int(129)
f0=open("./files/3/MapA0.dat","r")
nw=sum(1 for line in f0)
fil=np.zeros((nx*ny , nb))
fil[:nw,:]=np.loadtxt("./files/3/MapA0.dat")


nb2=int(68)
f0=open("./files/3/MapB0.dat","r")
nw=sum(1 for line in f0)
gil=np.zeros((nx*ny , nb2))
gil[:nw,:]=np.loadtxt("./files/3/MapB0.dat")

Netot=0.0; 
fil[:,38:40]=fil[:,38:40]*yr # in the unit of mas per year 
fil[:,58]=fil[:,58]*yr # in the unit of per year --->  epsilon(tE)/tE
###################################################################################################################


for i in range(nw): 
    nl1,nl2, tE0,tE1, RE0,RE1, piE0,piE1, tetE0,tetE1, Vt0, Vt1 , u00,u01, Ml0,Ml1, opdB0,opdB1, Dl0,Dl1 ,Ds0,Ds1= fil[i,0:22]
    vl0,vl1, vs0,vs1, mbs0,mbs1, fb0,fb1, fwhm0,fwhm1,vsn0,vsn1, DelT0,DelT1,struc0,struc1,murel0,murel1,Map0,Map1,nbl0=fil[i,22:43]
    nbl1,Ext0,Ext1,ros0,ros1,Extv0,Extv1,catg0,catg1,opdP0,opdP1, EffiB, EffiP, EffiD, EffiL, EFF, Gamma7, Neven= fil[i,43:61]
    lRostart, lNstart, nstart, Eru0, ErtE, Erfb, ErpiE, ErtetE, Erml, Erdl, Ermul, Ermus, Eros, dchiP, dchiA=  fil[i,61:76]
    fra1,fra2,fra3, Npeak,mag1,mag2,Sepa,ndd, nsim, numd0 , numd1 , nerr ,right, decli, fraBD0,fraBD1,fraMS0,fraMS1= fil[i,76:94]
    fraBH0,fraBH1,color00,color01,color10 ,color11, color20 ,color21, mag00 ,mag01 ,mag10, mag11, mag20,mag21,mag30= fil[i,94:109]
    mag31,mag40,mag41,mag50,mag51, fra0,fra1,fra2,fra3,fra4,fra5,fra6,fra7,fra8,fra9,argm0,argm1 ,color30,color31,flag=fil[i,109:]
   
    nl1, nl2, DtE0,DtE1, DRE0,DRE1, DpiE0,DpiE1, DtetE0,DtetE1, DVt0,DVt1, Du00,Du01, DMl0,DMl1 = gil[i,0:16]
    Dopd0, Dopd1, DDl0,DDl1, DDs0,DDs1, Dvl0, Dvl1, Dvs0, Dvs1, Dmbs0, Dmbs1, Dfb0, Dfb1 =  gil[i,16:30]
    Dfwhm0, Dfwhm1, Dvsn0, Dvsn1, DDelT0, DDelT1, Dmurel0, Dmurel1, DMap0, DMap1, Dnbl0,Dnbl1=gil[i,30:42]
    DExt0,DExt1, Dros0,Dros1, DExtv0,DExtv1,Dcol00 ,Dcol01, Dcol10 ,Dcol11 ,Dcol20 ,Dcol21   =gil[i,42:54]
    Dmag00,Dmag01,Dmag10,Dmag11,Dmag20,Dmag21,Dmag30,Dmag31,Dmag40,Dmag41,Dmag50,Dmag51,Dcol30,Dcol31=gil[i,54:68]
    l= -100.0+nl1*0.25
    b= -20.0 +nl2*0.25
    Netot+= Neven*ddeg*ddeg
    
    
    ############################################# >> Table2.txt ###################################
    tes=np.array([l,b, tE1,DtE1, piE1,DpiE1, tetE1,DtetE1, np.log10(ros1),np.log10(Dros1), Ml1,DMl1, Dl1,DDl1, Ds1,DDs1, 
    u01,Du01, murel1/yr,Dmurel1/yr, vl1,Dvl1])#22
    fil2=open("Table2.txt","a+")
    np.savetxt(fil2, tes.reshape((-1,22)),fmt="$%.2f$&$%.2f$&$%.1f\pm%.1f$&$%.2f\pm%.2f$&$%.2f\pm%.2f$&$%.2f\pm%.2f$&$%.2f\pm%.2f$&$%.1f\pm%.1f$&$%.1f\pm%.1f$&$%.1f\pm%.1f$&$%.1f\pm%.1f$&$%.1f\pm%.1f$"+r"\\")
    fil2.close()
    ############################################# >> Table3.txt ###################################
    tes=np.array([l,b,mag01,Dmag01,mag11,Dmag11,mag21,Dmag21,mag31,Dmag31,mag41,Dmag41,mag51,Dmag51,Ext1,DExt1,Extv1,DExtv1])#18
    fil2=open("Table3.txt","a+")
    np.savetxt(fil2,tes.reshape((-1,18)),fmt="$%.2f$&$%.2f$&$%.2f\pm%.2f$&$%.2f\pm%.2f$&$%.2f\pm%.2f$&$%.2f\pm%.2f$&$%.2f\pm%.2f$&$%.2f\pm%.2f$&$%.2f\pm%.2f$&$%.2f\pm%.2f$"+r"\\")##18
    fil2.close()
    ############################################# >> Table4.txt ###################################
    tes=np.array([l,b,np.log10(opdP0)-6.0,np.log10(opdB0)-6.0,np.log10(Gamma7)-7.0,np.log10(Neven),lRostart,lNstart,np.log10(nstart), 
    EffiB, EffiP, EffiL,ndd])#13
    fil2=open("Table4.txt","a+")
    np.savetxt(fil2,tes.reshape((-1,13)),fmt="$%.2f$&$%.2f$&$%.2f$&$%.2f$&$%.3f$&$%.2f$&$%.3f$&$%.3f$&$%.3f$&$%.2f$&$%.2f$&$%.2f$&$%d$"+r"\\")
    fil2.close()
    ############################################# >> Table5.txt ###################################
    tes=np.array([l,b,fraBD1,fraMS1,fraBH1,fra0,fra1,fra2,fra3,fra4,fra5,fra6,fra7,fra8])#14
    fil2=open("Table5.txt","a+")
    np.savetxt(fil2,tes.reshape((-1,14)),fmt="$%.2f$&$%.2f$&$%.2f:%.2f:%.2f$&$%.2f$&$%.2f$&$%.2f$&$%.2f$&$%.2f$&$%.2f$&$%.2f$&$%.2f$&$%.2f$"+r"\\")
    fil2.close()

print("\n*********>>>> Total number of events:  ", Netot, np.sum(fil[:nw,60])*area)     
###################################################################################################################
###################################################################################################################
kil=np.zeros((nx*ny , nb));  k=0
for i in range(nw): 
    if(fil[i,3]>0.0 and fil[i,19]>0.0 and fil[i,15]>0.0):  
        kil[k,:]= fil[i,:]
        k+=1
 
print("\n\n<tE>:", np.mean(kil[:k,3]), "\pm",   np.std(kil[:k,3]))
print("<piEE>:  ", np.mean(kil[:k,7]), "\pm",   np.std(kil[:k,7]))
print("<tetE>:  ", np.mean(kil[:k,9]), "\pm",   np.std(kil[:k,9]))
print("<ros>:   ", np.mean(kil[:k,47]), "\pm",  np.std(kil[:k,47]))
print("<Ml>:   ",  np.mean(kil[:k,15]), "\pm",  np.std(kil[:k,15]))
print("<Dl>:   ",  np.mean(kil[:k,19]), "\pm",  np.std(kil[:k,19]))
print("<murel>: ", np.mean(kil[:k,39]), "\pm",  np.std(kil[:k,39]))
print("<u0>:    ", np.mean(kil[:k,13]), "\pm",  np.std(kil[:k,13]))
print("<opt_P>:" , np.mean(kil[:k,52]), "\pm",  np.std(kil[:k,52]))
print("<opt_B>: ", np.mean(kil[:k,16]), "\pm",  np.std(kil[:k,16]))
print("<Ne>:  ",   np.mean(kil[:k,60]), "\pm",  np.std(kil[:k,60]))
print("<Gamma>:  ",np.mean(kil[:k,59]), "\pm",  np.std(kil[:k,59]))
print("<Nstar>:  ",np.mean(kil[:k,63]), "\pm",  np.std(kil[:k,63]))
print("\nFraDET BD:MS:BH: ", np.mean(kil[:k,91]),  np.mean(kil[:k,93]),  np.mean(kil[:k,95]) )
#print("FraSIM BD:MS:BH: ",   np.mean(kil[:k,90]),  np.mean(kil[:k,92]),  np.mean(kil[:k,94]) )
print("\nMax_opticalDP: ", np.max(kil[:k,52]), np.argmax(kil[:k,52]), kil[np.argmax(kil[:k,52]),0:2] )
print("\nMax_opticalDB: ", np.max(kil[:k,16]), np.argmax(kil[:k,16]), kil[np.argmax(kil[:k,16]),0:2] )
print("\nMax_Nevent: ",    np.max(kil[:k,60]), np.argmax(kil[:k,60]), kil[np.argmax(kil[:k,60]),0:2] )
print("\nMax_Gamma : ",    np.max(kil[:k,59]), np.argmax(kil[:k,59]), kil[np.argmax(kil[:k,59]),0:2] )
print("\nMax_Nstar : ",    np.max(kil[:k,63]), np.argmax(kil[:k,63]), kil[np.argmax(kil[:k,63]),0:2] )
print("\nMax_murel : ",    np.max(kil[:k,39]), np.argmax(kil[:k,39]), kil[np.argmax(kil[:k,39]),0:2] )

print("\n<varepsilon_D_P>[\%]:  ",np.mean(kil[:k,55]), "\pm",  np.std(kil[:k,55]))
print("\n<varepsilon_D_B>[\%]:  ",np.mean(kil[:k,54]), "\pm",  np.std(kil[:k,54]))
print("\n<varepsilon_L>[\%]:  ",np.mean(kil[:k,57]), "\pm",  np.std(kil[:k,57]))
print("\nMax_varepsilonL: ",    np.max(kil[:k,57]), np.argmax(kil[:k,57]), kil[np.argmax(kil[:k,57]),0:2] )


print("u-mag: ", np.mean(kil[:k,103]), "\pm", np.std(kil[:k,103]), np.min(kil[:k,103]), np.max(kil[:k,103]))
print("g-mag: ", np.mean(kil[:k,105]), "\pm", np.std(kil[:k,105]), np.min(kil[:k,105]), np.max(kil[:k,105]))
print("r-mag: ", np.mean(kil[:k,107]), "\pm", np.std(kil[:k,107]), np.min(kil[:k,107]), np.max(kil[:k,107]))
print("i-mag: ", np.mean(kil[:k,109]), "\pm", np.std(kil[:k,109]), np.min(kil[:k,109]), np.max(kil[:k,109]))
print("z-mag: ", np.mean(kil[:k,111]), "\pm", np.std(kil[:k,111]), np.min(kil[:k,111]), np.max(kil[:k,111]))
print("y-mag: ", np.mean(kil[:k,113]), "\pm", np.std(kil[:k,113]), np.min(kil[:k,113]), np.max(kil[:k,113]))
##################################################################################################################
##################################################################################################################
fil[:,2]=np.log10(fil[:,2]+np.min(np.where(fil[:,2]==0,np.inf,fil[:,2])))## log(tE)
fil[:,3]=np.log10(fil[:,3]+np.min(np.where(fil[:,3]==0,np.inf,fil[:,3])))## log(tE)
fil[:,6]=np.log10(fil[:,6]+np.min(np.where(fil[:,6]==0,np.inf,fil[:,6])))## log(pi)
fil[:,7]=np.log10(fil[:,7]+np.min(np.where(fil[:,7]==0,np.inf,fil[:,7])))## log(pi)
fil[:,8]=np.log10(fil[:,8]+np.min(np.where(fil[:,8]==0,np.inf,fil[:,8])))## log(tetE)
fil[:,9]=np.log10(fil[:,9]+np.min(np.where(fil[:,9]==0,np.inf,fil[:,9])))## log(tetE)
fil[:,16]=np.log10(fil[:,16]+np.min(np.where(fil[:,16]==0,np.inf,fil[:,16])))## log(potical_depth_baseline)
fil[:,17]=np.log10(fil[:,17]+np.min(np.where(fil[:,17]==0,np.inf,fil[:,17])))## log(optical_depth_baseline)
fil[:,38]=np.log10(fil[:,38]+np.min(np.where(fil[:,38]==0,np.inf,fil[:,38])))## log(mu_rel)
fil[:,39]=np.log10(fil[:,39]+np.min(np.where(fil[:,39]==0,np.inf,fil[:,39])))## log(mu_rel)
fil[:,46]=np.log10(fil[:,46]+np.min(np.where(fil[:,46]==0,np.inf,fil[:,46])))## log(rho*)
fil[:,47]=np.log10(fil[:,47]+np.min(np.where(fil[:,47]==0,np.inf,fil[:,47])))## log(rho*)
fil[:,52]=np.log10(fil[:,52]+np.min(np.where(fil[:,52]==0,np.inf,fil[:,52])))## log(tau)
fil[:,53]=np.log10(fil[:,53]+np.min(np.where(fil[:,53]==0,np.inf,fil[:,53])))## log(tau)
fil[:,58]=np.log10(fil[:,58]+np.min(np.where(fil[:,58]==0,np.inf,fil[:,58])))## log(epsilon(tE)/tE)
fil[:,59]=np.log10(fil[:,59]+np.min(np.where(fil[:,59]==0,np.inf,fil[:,59])))## log(Gamma)
fil[:,60]=np.log10(fil[:,60]+np.min(np.where(fil[:,60]==0,np.inf,fil[:,60])))## log(Ne)
fil[:,63]=np.log10(fil[:,63]+np.min(np.where(fil[:,63]==0,np.inf,fil[:,63])))## log(NstarRubin)
fil[:,64]=np.log10(fil[:,64]+np.min(np.where(fil[:,64]==0,np.inf,fil[:,64])))## log(delta_u0)
fil[:,65]=np.log10(fil[:,65]+np.min(np.where(fil[:,65]==0,np.inf,fil[:,65])))## log(delta_t_E)
fil[:,66]=np.log10(fil[:,66]+np.min(np.where(fil[:,66]==0,np.inf,fil[:,66])))## log(delta_f_b)
fil[:,67]=np.log10(fil[:,67]+np.min(np.where(fil[:,67]==0,np.inf,fil[:,67])))## log(delta_pi_E)
fil[:,68]=np.log10(fil[:,68]+np.min(np.where(fil[:,68]==0,np.inf,fil[:,68])))## log(delta_tet_E)
fil[:,69]=np.log10(fil[:,69]+np.min(np.where(fil[:,69]==0,np.inf,fil[:,69])))## log(delta_M_l)
fil[:,70]=np.log10(fil[:,70]+np.min(np.where(fil[:,70]==0,np.inf,fil[:,70])))## log(delta_D_l)
fil[:,71]=np.log10(fil[:,71]+np.min(np.where(fil[:,71]==0,np.inf,fil[:,71])))## log(delta_mu_l)
fil[:,72]=np.log10(fil[:,72]+np.min(np.where(fil[:,72]==0,np.inf,fil[:,72])))## log(delta_mu_s)
fil[:,73]=np.log10(fil[:,73]+np.min(np.where(fil[:,73]==0,np.inf,fil[:,73])))## log(delta_ro_s)
fil[:,74]=np.log10(fil[:,74]+np.min(np.where(fil[:,74]==0,np.inf,fil[:,74])))## log(Dchi2P)
fil[:,75]=np.log10(fil[:,75]+np.min(np.where(fil[:,75]==0,np.inf,fil[:,75])))## log(Dchi2A)

pointx=[]
pointy=[]
npp=0 

param=np.zeros((ny,nx, nb))
k=0
for i in range(nx): 
    for j in range(ny):
        if(fil[k,83]>=8 and fil[k,83]<11): 
            pointx.append(i)
            pointy.append(j)  
            npp+=1
        param[j,i,:]= fil[k,:]
        k+=1

#print(pointx, pointy)        
fil2=open("point.txt","w")
for i in range(npp): 
    np.savetxt(fil2,np.array([pointx[i],pointy[i]]).reshape((-1,2)),fmt="%d    %d")
fil2.close()   

##############################################################################################################################
##############################################################################################################################        
OPT1=np.zeros((nx,2)); 
NEE1=np.zeros((nx,2));   
OPT2=np.zeros((ny,2));   
NEE2=np.zeros((ny,2));  
xax1=np.zeros((nx));
xax2=np.zeros((ny));       
for i in range(nx):              
    OPT1[i,:]=np.mean(param[:,i,52]), np.std(param[:,i,52])
    NEE1[i,:]=np.mean(param[:,i,60]), np.std(param[:,i,60])
    xax1[i]=-100.0+i*0.25
for i in range(ny):     
    OPT2[i,:]=np.mean(param[i,:,52]), np.std(param[i,:,52])
    NEE2[i,:]=np.mean(param[i,:,60]), np.std(param[i,:,60])
    xax2[i]=-20.0+i*0.25
        
        
##############################################################################################################################
nz=580
for i in range(2):  
    plt.clf()
    plt.clf()
    fig=plt.figure(figsize=(8,6))
    ax1=plt.gca()  
    ax2 = ax1.twiny() 
    ax1.set_xlim(-100.0,100.0)
    ax2.set_xlim(-20.0, 20.0)
    if(i==0):
        ax1.plot(xax1[:nz], OPT1[:nz,0], color="darkblue",lw=2.0)
        ax1.fill_between(xax1[:nz], OPT1[:nz,0]-OPT1[:nz,1], OPT1[:nz,0]+OPT1[:nz,1], color="blue", alpha=0.37)      
        ax1.set_xlabel(r"$\rm{Galactic}~\rm{Lonitude}(\rm{deg})$", color='b', fontsize=18)
        ax1.tick_params(axis='x', labelcolor='b')
        ax2.plot(xax2, OPT2[:,0], color="darkred",lw=2.0)
        ax2.fill_between(xax2, OPT2[:,0]-OPT2[:,1], OPT2[:,0]+OPT2[:,1], color="r", alpha=0.37)      
        ax2.set_xlabel(r"$\rm{Galactic}~\rm{Latitude}(\rm{deg})$", color='r', fontsize=18)
        ax2.tick_params(axis='x', labelcolor='r')
        #ax1.set_yticks(fontsize=18, rotation=0)
        #ax1.set_xticks(fontsize=18, rotation=0)
        #ax2.set_xticks(fontsize=18, rotation=0)
        ax1.set_ylabel(r"$\overline{\tau_{\rm{p}}}\times10^{6}$",fontsize=20)
    if(i==1):
        ax1.plot(xax1[:nz], NEE1[:nz,0], color="darkblue",lw=2.0)
        ax1.fill_between(xax1[:nz], NEE1[:nz,0]-NEE1[:nz,1], NEE1[:nz,0]+NEE1[:nz,1], color="blue", alpha=0.37)      
        ax1.set_xlabel(r"$\rm{Galactic}~\rm{Lonitude}(\rm{deg})$", color='b', fontsize=20)
        ax1.tick_params(axis='x', labelcolor='b')
        ax2.plot(xax2, NEE2[:,0], color="darkred",lw=2.0)
        ax2.fill_between(xax2, NEE2[:,0]-NEE2[:,1], NEE2[:,0]+NEE2[:,1], color="r", alpha=0.37)      
        ax2.set_xlabel(r"$\rm{Galactic}~\rm{Latitude}(\rm{deg})$", color='r', fontsize=20)
        ax2.tick_params(axis='x', labelcolor='r')
        ax2.set_ylim(-20.0,120.0)
        #ax1.set_yticks(fontsize=18, rotation=0)
        #ax1.set_xticks(fontsize=18, rotation=0)
        #ax2.set_xticks(fontsize=18, rotation=0)
        ax1.set_ylabel(r"$\overline{dN_{\rm{e}}/d\Omega}(\rm{deg}^{-2})$",fontsize=20)
    #plt.grid("True")
    #plt.grid(linestyle='dashed')
    fig= plt.gcf()
    fig.tight_layout()
    fig.savefig("./files/Hist/NETAU{0:d}.jpg".format(i),dpi=200)   
    
        
##############################################################################################################################
##############################################################################################################################        
for i in range(nb):   
    minn=np.min(param[:,:,i])
    maxx=np.max(param[:,:,i])
    if(i==2 or i==3): minn=1.7;  maxx=2.6
    
    if(i==102 or i==103): minn=22.9;  maxx=45.0-7## u-band
    if(i==104 or i==105): minn=20.1;  maxx=40.5-5.5##g-band
    if(i==106 or i==107): minn=19.9;  maxx=35.6-4.0##r-band
    if(i==108 or i==109): minn=20.0;  maxx=30.5##i-band
    if(i==110 or i==111): minn=21.3;  maxx=30.3##z-band 
    if(i==112 or i==113): minn=18.8;  maxx=31.4##y-band
   
    plotMap(param[:,:,i],pointx, pointy, i, minn, maxx,  titl[i])
print(">>>>>>>>>>>>>>>>>>>>Maps were plotted <<<<<<<<<<<<<<<<<<<<<<")   
print("************************************************************")   


##############################################################################################################################
####                       HISTOGRAMS 
##############################################################################################################################   
for i in range(nb):
    plt.clf()
    plt.clf()
    fig=plt.figure(figsize=(8,6))
    ax1=plt.gca()    
    plt.hist(fil[:nw,i],30,histtype='bar',ec='darkgreen',facecolor='g',alpha=0.4,rwidth=1.5)
    ax1.set_xlabel(str(titl[i]), fontsize=18, labelpad=0.1)
    yvals=ax1.get_yticks()
    ax1.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/nw)) for x in yvals]) 
    yvals = ax1.get_yticks()
    ax1.set_ylim([ np.min(yvals), np.max(yvals) ])
    plt.yticks(fontsize=18, rotation=0)
    plt.xticks(fontsize=18, rotation=0)
    ax1.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=18,labelpad=0.1)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig= plt.gcf()
    fig.tight_layout()
    fig.savefig("./files/Hist/extra/Histo{0:d}.jpg".format(i),dpi=200)    
print(">>>>>>>>>>>>>>Histograms were plotted <<<<<<<<<<<<<<<<<<<<<<")   
print("************************************************************")   




##############################################################################################################################
##############################################################################################################################
#                            EfB0.dat Efficiency versus some parameters 
##############################################################################################################################
##############################################################################################################################
#l.NdtE(0),NstE(1), NdMl(2), NsMl(3), Ndpi(4), Nspi(5), Ndu0(6), Nsu0(7), Ndmb(8), Nsmb(9),Ndfb(10), Nsfb(11), Ndmu(12),Nsmu(13), Nhalo(14), Nhalo(15), Ndisk(16), Ndisk(17), ndteBD(18) ,nsteBD(19), ndteMS(20), nsteMS(21), ndteBH(22), nsteBH(23)
 
 
parm=[r"$t_{\rm{E}}(\rm{days})$", r"$M_{\rm{l}}(M_{\odot})$", r"$\log_{10}[\pi_{\rm{rel}}(\rm{mas})]$", r"$u_{0}$", r"$m_{\rm{base}}(\rm{mag}), r-\rm{band}$", r"$f_{\rm{b}}$", r"$\mu_{\rm{rel}}(\rm{mas}/\rm{years})$", r"$t_{\rm{E}}(\rm{days})$"]##8

gg=101
efft=np.zeros((gg+1,24))
efft=np.loadtxt("./files/3/EfB0.dat")
param=np.zeros((10,gg,2))
tE_min=1.0;
tE_max=1826;
Ml_min=16.0*Mjupiter/Msun;
Ml_max=98.0;
pi_min=np.log10(1.0/8.5-1.0/9.0); ##pi_rel
pi_max=np.log10(1.0/0.05-1.0/9.0);  
u0_min=0.0; 
u0_max=1.0;
mb_min=15.1+0.5;##baseline magnitude in r-LSST band 
mb_max=24.0-0.5;  
fb_min=0.0; 
fb_max=1.0;  
mu_min=1.0;#///mu_relative (mas/years)
mu_max=100.0;  

for i in range(gg): 
    param[0,i,0:2]=float(tE_min+(tE_max-tE_min)*i/gg),  float(efft[i,0] *100.0/efft[i,1])
    param[1,i,0:2]=float(Ml_min+(Ml_max-Ml_min)*i/gg),  float(efft[i,2] *100.0/efft[i,3])
    param[2,i,0:2]=float(pi_min+(pi_max-pi_min)*i/gg),  float(efft[i,4] *100.0/efft[i,5])
    param[3,i,0:2]=float(u0_min+(u0_max-u0_min)*i/gg),  float(efft[i,6] *100.0/efft[i,7])
    param[4,i,0:2]=float(mb_min+(mb_max-mb_min)*i/gg),  float(efft[i,8] *100.0/efft[i,9])
    param[5,i,0:2]=float(fb_min+(fb_max-fb_min)*i/gg),  float(efft[i,10]*100.0/efft[i,11])
    param[6,i,0:2]=float(mu_min+(mu_max-mu_min)*i/gg),  float(efft[i,12]*100.0/efft[i,13])
    param[7,i,0:2]=float(tE_min+(tE_max-tE_min)*i/gg),  float(efft[i,18]*100.0/efft[i,19])
    param[8,i,0:2]=float(tE_min+(tE_max-tE_min)*i/gg),  float(efft[i,20]*100.0/efft[i,21])
    param[9,i,0:2]=float(tE_min+(tE_max-tE_min)*i/gg),  float(efft[i,22]*100.0/efft[i,23])
    
print("Detection efficiency for Halo-lensing: ", float(efft[0,14]*100.0/efft[0,15]))    
print("Detection efficiency for Disk-lensing: ", float(efft[0,16]*100.0/efft[0,17]))    
    
#l.NdtE(0),NstE(1), NdMl(2), NsMl(3), Ndpi(4), Nspi(5), Ndu0(6), Nsu0(7), Ndmb(8), Nsmb(9),Ndfb(10), Nsfb(11), Ndmu(12),Nsmu(13), Nhalo(14), Nhalo(15), Ndisk(16), Ndisk(17), ndteBD(18) ,nsteBD(19), ndteMS(20), nsteMS(21), ndteBH(22), nsteBH(23)
    
for i in range(8): 
    plt.cla()
    plt.clf()
    fig=plt.figure(figsize=(8,6))
    if(i<7):
        plt.plot(param[i,:,0], param[i,:,1],"k-", lw=1.6)
    else: 
        plt.plot(param[0,:,0], param[0,:,1],"k-",  lw=2.6, label=r"$\rm{Total}~\rm{lenses}$")
        plt.plot(param[0,:,0], param[7,:,1],"g--", lw=1.1, label=r"$\rm{Brown}~\rm{Dwarfs}$")
        plt.plot(param[0,:,0], param[8,:,1],"r:",  lw=1.1, label=r"$\rm{Main}~\rm{Sequences}$")
        plt.plot(param[0,:,0], param[9,:,1],"m.-", lw=1.1, label=r"$\rm{Black}~\rm{Holes}$")
    plt.xlabel(str(parm[i]),fontsize=18)
    plt.ylabel(r"$\varepsilon_{\rm{L}}[\%]$",fontsize=18)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    if(i==7): 
        plt.legend()
        plt.legend(loc='best',fancybox=True, shadow=False, prop={"size":16},framealpha=0.65)
    fig= plt.gcf()
    fig.tight_layout()
    fig.savefig("./files/Hist/Effi{0:d}.jpg".format(i),dpi=200) 






##############################################################################################################################
#                        LpLMC.dat or LpSMC.dat    Parameters of detected lensing events 
##############################################################################################################################

#nl1(0),nl2(1),ndd(2), Dl(3), Ml(4), mu_lens(5), u0(6), catgLens(7), Ds(8),mbaseG(9), mbaseR(10),mbaseI(11),mbaseZ(12), fbR(13),tE(14),piE(15), tetE(16), ros(17), murel(18), opt(19),ErtE(20),Erfb(21),ErpiE(22),ErtetE(23), ErMl(24),ErDl(25),ErMuL(26),ErRos(27) ##28
   
nam4=[r"$M_{\rm{l}}$",r"$D_{\rm{l}}$",r"$t_{\rm{E}}$",r"$u_{0}$", r"$m_{\rm{base}}$", r"$N_{\rm{data}}$", 
r"$\mu_{\rm{rel}}$",r"$\delta_{t_{\rm{E}}}$",r"$\delta_{\pi_{\rm{E}}}$",r"$\delta_{\theta_{\rm{E}}}$",r"$\delta_{\rho_{\star}}$"]##11

pd.options.display.precision = 3
df= pd.read_csv("./files/3/Lpar0.dat",sep=" ",skipinitialspace=True,header=None,usecols=[4,3,14,6,10,2,18,20,22,23,27],names=nam4)
corrM = df.corr()
fig, ax = plt.subplots(figsize=(10, 8))
corrM.style.background_gradient(cmap='coolwarm').format(precision=3)
df.style.format("{:.3f}")
ax=sns.heatmap(corrM, annot=True, xticklabels=nam4, yticklabels=nam4, annot_kws={"size":15.5},square=True,linewidth=1.0, cbar_kws={"shrink":.99}, linecolor="k",fmt=".1f", cbar=True, vmax=1, vmin=-1, center=0.0, robust=True)#ax=None,
cbar=ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=17.5, pad=0.0)
plt.xticks(rotation=35,horizontalalignment='right',fontweight='light', fontsize=15)
plt.yticks(rotation=0, horizontalalignment='right',fontweight='light', fontsize=15)
plt.title(r"$\rm{Correlation}~\rm{Matrix}$", fontsize=17)
pd.set_option('display.precision', 3)
fig.tight_layout()
plt.savefig("./files/Hist/Correlation.jpg",dpi=200)
print("**** Correlation matrix was calculated ******** ")



##############################################################################################################################
##############################################################################################################################  
#     Table1.txt saving the fractions  
##############################################################################################################################
##############################################################################################################################
#nl1(0),nl2(1),ndd(2), Dl(3), Ml(4), mu_lens(5), u0(6), catgLens(7), Ds(8),mbaseG(9), mbaseR(10),mbaseI(11),mbaseZ(12), fbR(13),tE(14),piE(15), tetE(16), ros(17), murel(18), opt(19),ErtE(20),Erfb(21),ErpiE(22),ErtetE(23), ErMl(24),ErDl(25),ErMuL(26),ErRos(27) ##28

f0=open("./files/3/Lpar0.dat","r")
ns=sum(1 for line in f0)
par=np.zeros((ns,28))
par=np.loadtxt("./files/3/Lpar0.dat")    

thre=np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.50])##10
frace=np.zeros((9,10))

for i in range(ns):
    for k in range(10): #threshold
        for j in range(8):#Relative_Errors  
            if(par[i,20+j]<thre[k]):    frace[j,k]+=1.0
        if(par[i,24]<thre[k] and par[i,25]<thre[k] and par[i,26]<thre[k]): frace[8,k]+=1.0
        
for k in range(10): 
    for j in range(9):  
        frace[j,k]=float(frace[j,k]*100.0/ns)
fil5=open("Table1.txt","w")             
for i in range(10):
    sav=np.array([thre[i]*100,frace[0,i],frace[1,i],frace[2,i],frace[3,i],frace[4,i],frace[5,i],frace[6,i],
    frace[7,i],frace[8,i]])        
    np.savetxt(fil5,sav.reshape((-1,10)),fmt="$%d$&$%.2f$&$%.2f$&$%.2f$&$%.2f$&$%.2f$&$%.2f$&$%.2f$&$%.2f$&$%.2f$"+r"\\")  
fil5.close()


##############################################################################################################################
##    population  
##############################################################################################################################   
#nl1(0),nl2(1),ndd(2), Dl(3), Ml(4), mu_lens(5), u0(6), catgLens(7), Ds(8),mbaseG(9), mbaseR(10),mbaseI(11),mbaseZ(12), fbR(13),tE(14),piE(15), tetE(16), ros(17), murel(18), opt(19),ErtE(20),Erfb(21),ErpiE(22),ErtetE(23), ErMl(24),ErDl(25),ErMuL(26),ErRos(27) ##28
nam5=[r"$\log_{10}[t_{\rm{E}}(\rm{days})]$",r"$\log_{10}[\pi_{\rm{E}}]$", r"$\log_{10}[\theta_{\rm{E}}(\rm{mas})]$", r"$\log_{10}[\rho_{\star}]$", r"$D_{\rm{l}}(\rm{kpc})$", r"$u_{0}$", r"$D_{\star}(\rm{kpc})$", r"$m_{\rm{base}}(\rm{mag}), r-\rm{band}$"]##8

pop=np.zeros((ns,8,4))
n1=int(0); n2=int(0); n3=int(0); n4=int(0);  

for i in range(ns):
    if(par[i,14]>0.0 and par[i,15]>0.0 and par[i,16]>0.0 and par[i,17]): 
        if(par[i,7]==1):pop[n1,:,0]=np.array([par[i,14],par[i,15],par[i,16],par[i,17],par[i,3],par[i,6],par[i,8],par[i,10]]); n1+=1
        if(par[i,7]==2):pop[n2,:,1]=np.array([par[i,14],par[i,15],par[i,16],par[i,17],par[i,3],par[i,6],par[i,8],par[i,10]]); n2+=1
        if(par[i,7]==3):pop[n3,:,2]=np.array([par[i,14],par[i,15],par[i,16],par[i,17],par[i,3],par[i,6],par[i,8],par[i,10]]); n3+=1
        pop[n4,:,3]=np.array([par[i,14],par[i,15],par[i,16],par[i,17],par[i,3],par[i,6],par[i,8],par[i,10]]); n4+=1

n1=int(n1);  n2=int(n2);  
n3=int(n3);  n4=int(n4)    

pop[:n1,0:4,0]=np.log10(pop[:n1,0:4,0])
pop[:n2,0:4,1]=np.log10(pop[:n2,0:4,1])
pop[:n3,0:4,2]=np.log10(pop[:n3,0:4,2])
pop[:n4,0:4,3]=np.log10(pop[:n4,0:4,3]) 


ave=np.zeros((8,4))
std=np.zeros((8,4))
for i in range(8): 
    ave[i,0]=np.mean(pop[:n1,i,0]);  std[i,0]=np.std(pop[:n1,i,0]);  
    ave[i,1]=np.mean(pop[:n2,i,1]);  std[i,1]=np.std(pop[:n2,i,1])
    ave[i,2]=np.mean(pop[:n3,i,2]);  std[i,2]=np.std(pop[:n3,i,2])
    ave[i,3]=np.mean(pop[:n4,i,3]);  std[i,3]=np.std(pop[:n4,i,3])
    ave[i,:]=np.round(ave[i,:],2);   std[i,:]=np.round(std[i,:],2);

for i in range(8): 
    plt.clf()
    plt.clf()
    fig=plt.figure(figsize=(8,6))
    ax1=plt.gca()    
    plt.hist(pop[:n1,i,0],30,histtype='step',lw=1.7,ls=':',color='g', label=r"$\rm{BDs},~$"+str(ave[i,0])+ r"$\pm$"+str(std[i,0]))
    plt.hist(pop[:n2,i,1],30,histtype='step',lw=1.7,ls='--',color='r',label=r"$\rm{MSs},~$"+str(ave[i,1])+ r"$\pm$"+str(std[i,1]))
    plt.hist(pop[:n3,i,2],30,histtype='step',lw=1.7,ls='-.',color='m',label=r"$\rm{BHs},~$"+str(ave[i,2])+ r"$\pm$"+str(std[i,2]))
    plt.hist(pop[:n4,i,3],30,histtype='step',lw=2.0,ls='-',color='k',label=r"$\rm{Total},~$"+str(ave[i,3])+r"$\pm$"+str(std[i,3]))    
    ax1.set_xlabel(str(nam5[i]), fontsize=18)
    yvals=ax1.get_yticks()
    ax1.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/ns)) for x in yvals]) 
    yvals = ax1.get_yticks()
    ax1.set_ylim([1.0, np.max(yvals) ])
    plt.yscale('log')
    plt.yticks(fontsize=18,rotation=0)
    plt.xticks(fontsize=18,rotation=0)
    ax1.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=18)
    #plt.grid("True")
    #plt.grid(linestyle='dashed')
    plt.legend()
    plt.legend(loc='best',fancybox=True, shadow=False,prop={"size":15}, framealpha=0.6)
    fig= plt.gcf()
    fig.tight_layout()
    fig.savefig("./files/Hist/Population{0:d}.jpg".format(i),dpi=200)   
print(">>>>>>>> population_plots were ploted <<<<<<<<<<<<<<<<<<<")

       
   
##############################################################################################################################
##    dependence of errors to parameters 
##############################################################################################################################   

laby=[r"$\log_{10}[\delta_{t_{\rm{E}}}]$",r"$\log_{10}[\delta_{f_{\rm{b}}}]$",r"$\log_{10}[\delta_{\pi_{\rm{E}}}]$",r"$\log_{10}[\delta_{\theta_{\rm{E}}}]$", r"$\log_{10}[\delta_{M_{\rm{l}}}]$", r"$\log_{10}[\delta_{D_{\rm{\rm{l}}}}]$", r"$\log_{10}[\delta_{\mu_{\rm{l}}}]$"]#7

labx=[r"$\log_{10}[M_{\rm{l}}(M_{\odot})]$",r"$D_{\rm{l}}(\rm{kpc})$",r"$\log_{10}[\rho_{\star}]$",r"$f_{\rm{b}}$",
r"$m_{\rm{base}}(\rm{mag})$", r"$\mu_{\rm{rel}}(\rm{mas}/\rm{years})$", 
r"$u_{0}$",r"$\rm{Galactic}~\rm{Longitude}(\rm{deg})$",r"$\rm{Galactic}~\rm{Latitude}(\rm{deg})$",
 r"$\log_{10}[\pi_{\rm{E}}]$",r"$\log_{10}[\theta_{\rm{E}}(\rm{mas})]$", r"$t_{\rm{E}}(\rm{days})$"]##12
#logMl(0), Dl(1), LogRho*(2),fb(3), mbase(4), Murel(6), u0(7), l(8), b(9) LogpiE(10), logtetE(11), tE(12)  
#r"$\log_{10}[t_{\rm{E}}(\rm{days})]$",
 


par[:ns,0]  =-100.0+par[:ns,0]*0.25 ## longitude 
par[:ns,1]  =-20.00+par[:ns,1]*0.25 ## latitude 
par[:ns,4]  =np.log10(par[:ns,4] +  np.min(np.where(par[:ns,4]==0,np.inf,par[:ns,4])  )) # log10(Ml)
par[:ns,17] =np.log10(par[:ns,17] + np.min(np.where(par[:ns,17]==0,np.inf,par[:ns,17])))#log10(rho*)
par[:ns,15] =np.log10(par[:ns,15] + np.min(np.where(par[:ns,15]==0,np.inf,par[:ns,15])))# log(piE)
par[:ns,16] =np.log10(par[:ns,16] + np.min(np.where(par[:ns,16]==0,np.inf,par[:ns,16]))) #log(tetE)
par[:ns,20:]=np.log10(par[:ns,20:])#piE,tetE,delta_all


xar= np.zeros((12,na,19))##X-axis(12),  Ngrid(14),  Y-axis(7)
xar[0,:,0]= np.array([-2.0, -1.6,  -1.4,  -1.1,  -0.8, -0.5, -0.2,  0.1,  0.4, 0.7, 1.0,  1.3, 1.6, 1.9])#14 log10(Ml)
xar[1,:,0]= np.array([0.5,  0.9 ,  1.3,  1.7, 2.1, 2.5, 2.9, 3.3,  3.7, 4.0, 4.4, 4.7, 5.0, 5.5])#DL
xar[2,:,0]= np.array([-3.5, -3.35, -3.2,  -3.05, -2.90, -2.75, -2.60, -2.45,  -2.30, -2.2, -2.1, -2.0,-1.9, -1.8])#log10(Rho_*)  
xar[3,:,0]= np.array([0.0, 0.08, 0.16,  0.24, 0.32, 0.4, 0.48, 0.56,  0.64, 0.72, 0.8, 0.88, 0.95,  1.0])#fb     
xar[4,:,0]= np.array([17.0,  17.5,  18.0,  18.5, 19.0, 19.5, 20.0, 20.5,  21.0, 21.5, 22.0, 22.5, 23.0, 24.0])#mbase
#xar[5,:,0]= np.array([-0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7,  0.9,  1.1, 1.3, 1.5,  1.7, 1.9,  2.1])#log10(tE)
xar[5,:,0]= np.array([0.0,  2.5, 5.0, 7.5,10.0, 12.5, 15.0, 17.5, 20.0,22.5, 25.0,  27.5, 30.0,32.5])#mu_rel
xar[6,:,0]= np.array([0.0, 0.08, 0.16,  0.24, 0.32, 0.4, 0.48, 0.56,  0.64, 0.72, 0.8, 0.88, 0.96,  1.0])#u0
xar[7,:,0]= np.array([-75.0,-64, -53, -42, -31, -20 ,-9, 2, 13, 24, 35, 46, 57, 68])#longitude
xar[8,:,0]= np.array([-20.0,-17.0, -14.0, -11.0,-8.0, -5.0, -2.0, 1.0, 4.0, 7.0, 10.0, 13.0, 16.0, 19.0])#latitude
xar[9,:,0]=np.array([-0.35, -0.27, -0.2,  -0.15, -0.1, -0.05, 0.0, 0.05,  0.1, 0.15, 0.2,  0.25,  0.3, 0.35])#log10(pi_E)
xar[10,:,0]=np.array([-0.2,-0.15, -0.1, -0.05, 0.0, 0.05,  0.1,  0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.46])#log10(tetE)
xar[11,:,0]=np.array([1.0, 15.0, 30.0,  45.0, 60.0, 75.0, 90.0,  105.0, 120.0, 135.0,150.0, 165.0, 180.0, 195.0])#tE

#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
#nl1(0),nl2(1),ndd(2), Dl(3), Ml(4), mu_lens(5), u0(6), catgLens(7), Ds(8),mbaseG(9), mbaseR(10),mbaseI(11),mbaseZ(12), fbR(13),tE(14),piE(15), tetE(16), ros(17), murel(18), opt(19),ErtE(20),Erfb(21),ErpiE(22),ErtetE(23), ErMl(24),ErDl(25),ErMuL(26),ErRos(27) ##28
t=[4,3, 17, 13, 10, 18, 6, 0, 1, 15, 16, 14]
for i in range(ns):#Npoint 
    for k in range(12):#X-axis  
        gg= funcx(xar[k,:,0], float(par[i,int(t[k])]) );  
        xar[k,gg,1:9] +=par[i,20:]
        xar[k,gg,9]   +=1.0
       
for k in range(12): 
    for i in range(na): 
        for j in range(8): 
            xar[k,i,j+1]=float(xar[k,i,j+1]/(xar[k,i,9]+eps))
#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
'''
for i in range(ns): 
    for k in range(12):  
        gg=funcx(xar[k,:,0], float(par[i,int(t[k])])  );  
        xar[k,gg,10:18]+=np.power(par[i,20:]-xar[k,gg,1:9],2.0);#STD 
        xar[k,gg,18]   +=1.0
        
for k in range(12): 
    for i in range(na): 
        for j in range(8): 
            xar[k,i,j+10]=np.sqrt(xar[k,i,j+10]/( xar[k,i,18]+eps))
'''            
###################################################################################################
CC=['darkred','darkblue','k','darkgreen','purple','teal','darkorange']
DD=['red','b','gray','g', 'pink','cyan','orange']

for j in range(12):#X-axis
    dx=abs(xar[j,1,0]-xar[j,0,0])/2.0
    plt.clf()
    plt.cla()
    fig=plt.figure(figsize=(8,6))
    
    for i in range(7):
        if(j==0 or j==3 or j==6 or j==9):   
            plt.plot(xar[j,1:13,0]+dx,xar[j,1:13,i+1],marker='o',mfc=DD[i],mec=CC[i],ms=7.3,ls='-',c=CC[i],lw=2.0,label=str(laby[i]))   
        else: 
            plt.plot(xar[j,:,0]+dx,xar[j,:,i+1],marker='o',mfc=DD[i],mec=CC[i],ms=7.3,ls='-',c=CC[i],lw=2.0,label=str(laby[i]))


    plt.xlabel(str(labx[j]), fontsize=22.5)
    plt.ylabel(r"$\rm{Logarithmic}~\rm{Relative}~\rm{Errors}$", fontsize=20)##,labelpad=0.1)
    #plt.xlim([ xa, xb])
    #plt.ylim([ ya, yb])
    plt.xticks(fontsize=20, rotation=0)
    plt.yticks(fontsize=20, rotation=0)  
    plt.legend(bbox_to_anchor = (1, 1))        
    plt.legend()
    plt.legend(loc='best',fancybox=True, shadow=False,prop={"size":16}, framealpha=0.25)
    fig=plt.gcf()
    fig.tight_layout()
    fig.savefig("./files/Hist/ErrorPH{0:d}.jpg".format(j),dpi=200)
###################################################################################################
    
      
