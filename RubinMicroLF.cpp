#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
#include <cstdlib>  
#include <ctime>
#include <cmath> 
#include <string>
#include <vector>
#include "VBBinaryLensingLibrary.h"

using std::cout;
using std::endl;
using std::cin;
using std::vector;



#define Nx 6 ///u0, tE, fb, piE, ksi, ros
#define Ny 4 ///tetE, mus1, mus2, piE
#define MIN(a,b) ((a) < (b) ? (a) : (b))


const double u0m=1.0;  
const double RAa=double(180.0/M_PI);
const double pi =M_PI; 
const double BF =double(2.0/3.0);//binary fraction
const double velocity=299792458.0;//velosity of light
const double Msun=1.98892*pow(10.,30); //in [kg].
const double Rsun=6.957*pow(10.0,8.0); ///solar radius [meter]
const double KP  =3.08568025*pow(10.,19); // in meter.
const double GN  =6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double AU  =1.4960*pow(10.0,11.0);
const double vro_sun=226.0;
const double VSunR=11.1;
const double VSunT=vro_sun*(1.00762+0.00712)+ 12.24;
const double VSunZ=7.25;
const double year =365.2425;//days 
const double eps  =double(0.000000000000005463263454624313654);
const double Mjupiter=1.89813*pow(10.0,27.0);///kilograms
const double Mmoon= 7.34767309*pow(10.0,22.0);//kilograms
///============================ Besancon constant ==========================///

const double Dsun   =8.0;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]  ={0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419,6.2/3.52112,4.0/2.27237,5.8/3.29525,4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const double Rv[4]  ={3.1, 2.5 ,3.1 ,3.1 };//Disk, Bulge, Thick, Halo 


///=============================LSST constant===============================///
//                      https://rubinobservatory.org/for-scientists/rubin-101/key-numbers
////                    https://arxiv.org/pdf/2603.23786v1  for FWHM values its table (3)
///   Saturation limit  https://arxiv.org/html/2507.01343v1

const int    M=6;    ///No. of filter  ugrizy  of LSST
const double Tobs   =10.0*year;///LSST observational time 10 years
const double delta2 =0.005;///systematic errors
const double gama[M]   ={0.038,0.039,0.039,0.039,0.039,0.039};///afunction of sky brightness and airmass in wavelengths.
const double sigma[M+1]={0.022,0.02,0.017,0.017,0.027,0.027,0.04};
const double thre[M+1] ={23.8, 24.5, 24.0, 23.4, 22.7, 22.0, 23.0};///depth of single visit in ugrizy+  K-bad(
const double satu[M+1] ={16.0 ,15.7, 15.1, 15.2, 14.5, 14.7 , 12.0};///saturation limit of single visit in ugrizy+ K-band
const double FWHM[M+1] ={1.48 ,1.17 ,1.12, 1.03, 1.11 , 1.01,  5.0*0.004};//LSST [arcsec] ugrizy + K-band
const double a0[M+1]   ={0.9429,1.0138,0.94027, 0.8139, 0.6641, 0.5703,      0.1615};//for calculating the extinction +K-band
const double b0[M+1]   ={1.9788, 0.5575,-0.2197, -0.4982, -0.6097, -0.5236, -0.1483};


const int coun=3000;
const int YZ  =int(3578);   //No.yzma.txt rows
const int met =70;   //No rows metal.txt  
const int Nel =int(7);//rows in "sigma_ELT.txt" (https://academic.oup.com/mnras/article/494/3/4413/5813442)
const int Na  =int(6);     //rows in "sigmaA_LSST.txt"
const int N1  =27004, N2=36558, N3=2612, N4=491;///CMD_BESANCON, thinD, bulge, thickD, halo
const int nq  =int(16);     //resu 
const int Nl  =int(3000); //GBbaselines.dat


const double tetp= double(M_PI/3.0);     //parallax
const double omegae=double(2.0*M_PI/year);//radian per day
const double vearth=omegae;                //radian per day 
const float  sig[2] ={+1.0, -1.0};  
const float  sig2[2]={+0.5, +1.0};  

const int    GG=int(101); 
const double tE_min=1.0;///days 
const double tE_max=1826;//days  

const double Ml_min=16.0*Mjupiter/Msun;///based on the Func_lens 
const double Ml_max=98.0;//in sun mass  

const double pi_min=log10(1.0/8.5-1.0/9.0); 
const double pi_max=log10(1.0/0.05-1.0/9.0);  

const double u0_min=0.0; 
const double u0_max=u0m;

const double mb_min=satu[2]+0.5;//baseline magnitude in r-LSST band 
const double mb_max=thre[2]-0.5;  

const double fb_min=0.0; 
const double fb_max=1.0;  

const double mu_min=1.0;///mu_relative (mas/years)
const double mu_max=100.0;  

////=====================#== Galactic bulge ====================================

const int    Num =int(800);
const double MaxD=double(12.5);///kpc
const double step=double(MaxD/Num/1.0);///step in kpc
const double ddeg=float(0.25);  
const double FoV =double(3.5/2.0);  //the radius of teh Rubin Field of View 

///============================================================================
struct source{
   int    nums,struc, cl;
   double nssim[Num], nsdet[Num];  
   double Ds,TET, FI, lat, lon, vs;
   double right, decli;
   double od_disk, od_ThD, od_bulge, od_halo[2], opt;
   double rho_disk[Num], rho_ThD[Num], rho_halo[Num][2], rho_stars[Num], rho_bulge[Num];
   double Rostar0[Num][2], Rostari[Num][2],Nstari[Num];
   double Nstart, Rostart[2], Romaxs[2],  Romins[2], nstart;
   double nsbl[M+1], blend[M+1], Fluxb[M+1], magb[M+1], Ai[M+1], Mab[M+1], Map[M+1], Mpeak[M];
   double mass, logT, typ, age, ros, Rstar;
   double mus1, mus2, mus;
   double xv, yv, zv; 
   double Av; 
   double SV_n1, LV_n1, VSun_n1;
   double SV_n2, LV_n2, VSun_n2;
   double pos1b, pos2b, pos1c, pos2c, def1a,  def2a,  def1c,  def2c;  
   double errM, errA, FWHM;  
   double ut, ut0, Astar, xi, ux, uy; 
   double fb[2], mbs[2];
};

struct lens{
   int    numl,struc;
   double Ml,Dl,vl,Vt, xls, u0, A0, mi1, mi2;
   double rhomaxl, tE, RE, t0, murel, DeltaT, murel1, murel2;
   double piE, pirel, tetE, pos1, pos2, mul1, mul2, mul; 
   double betal, betas, deltal, deltas, deltao; 
   double nstE[GG+1],    ndtE[GG+1];
   double nste[3][GG+1], ndte[3][GG+1];//BD, MS, BH
   double NstE[GG+1], NdtE[GG+1], tEs[GG+1], EFTE[GG+1];  
   double NsMl[GG+1], NdMl[GG+1], Mls[GG+1];  
   double Nspi[GG+1], Ndpi[GG+1], pis[GG+1];  
   double Nsu0[GG+1], Ndu0[GG+1], u0s[GG+1];  
   double Nsmb[GG+1], Ndmb[GG+1], mbs[GG+1];  
   double Nsfb[GG+1], Ndfb[GG+1], fbs[GG+1];  
   double Nsmu[GG+1], Ndmu[GG+1], mus[GG+1];  
   double Nhalo[2], Ndisk[2];  
   double timn[coun];
   double magn[coun];
   double soux[coun];
   double souy[coun];
   double errm[coun];
   double erra[coun];
   int    tele[coun]; 
   int catg; 
};

struct astromet{
   double Ve_n1, Ve_n2;
   double ue_n1, ue_n2;
};

struct CMD{  
   double logT_d[N1], mass_d[N1], Mk_d[N1],typ_d[N1], radi_d[N1]; int cl_d[N1],age_d[N1];///thin disk
   double logT_b[N2], mass_b[N2], Mk_b[N2],typ_b[N2], radi_b[N2]; int cl_b[N2],age_b[N2];/// bulge
   double logT_t[N3], mass_t[N3], Mk_t[N3],typ_t[N3], radi_t[N3]; int cl_t[N3],age_t[N3];///thick disk
   double logT_h[N4], mass_h[N4], Mk_h[N4],typ_h[N4], radi_h[N4]; int cl_h[N4],age_h[N4];/// halo
   double **Mab_d = new double*[N1];
   double **Mab_b = new double*[N2];
   double **Mab_t = new double*[N3];
   double **Mab_h = new double*[N4];
};

struct lsst{
   double maga[Na],  erra[Na];//LSST_astrometrc_error
   double mage[Nel], ere1[Nel], ere2[Nel];//ELT   
   int    filter[Nl];
   int    Nfilt[6],Nfilt2[6],argm[6]; 
   double RA, DEC, tim[Nl], sig5[Nl]; 
   int flag; 
};

struct covarian{
   double Era[Nx], Erb[Ny], sigmul1, sigmul2,f1, f2;   
   double magw, derm1[2], derm2[2], dera1[2], derb1[2], dera2[2], derb2[2]; 
   double derm1f, derm2f, dera1f, derb1f, dera2f, derb2f, diff; 
   double Delta1[Nx], Delta2[Ny], bb[2];
   double resu[nq]; 
   double inputA[Nx][Nx];  
   double inverA[Nx][Nx]; 
   double adjunA[Nx][Nx]; 
   double tempoA[Nx][Nx]; 
   double inputB[Ny][Ny];  
   double inverB[Ny][Ny]; 
   double adjunB[Ny][Ny]; 
   double tempoB[Ny][Ny]; 
   double A[5][5];  
   double summ[Nx][Nx]; 
   double deter;  
   int    sign;
   int    flagi;  
};

struct yfilter{
    double   Age[YZ]; 
    double     B[YZ];  
    double     M[YZ];   
    double    mm[YZ]; 
    int   number[met];
    int    count[met];   
    double Metal[met];
}; 

struct extinc{
   double dis[100];///distance 
   double Extks[100];///ks-band extinction
   double Aks;
};

///===================== FUNCTION ===========================================//
int    Funcu0(lens & l); 
int    FunctE(lens & l); 
int    FuncMl(lens & l); 
int    FuncPi(lens & l); 
int    FuncMu(lens & l); 
int    FuncMb(lens & l, double);
int    FuncFb(lens & l, double);
void   read_cmd(CMD & cm,yfilter & yf);
void   optical_depth(source & s); 
void   func_source(source & s,CMD & cm, extinc & ex);
void   func_lens(lens & l, source & s);
void   vrel(source & s, lens & l);
void   Disk_model(source & s, int);
void   ErrorCal(covarian & co, lens &l, source &s);
void   lightcurve(source & s, lens & l, astromet & as, double); 
void   FisherM(source & s, lens & l, astromet & as,  covarian & co, int,  VBBinaryLensing & vbb);
double ylsst(yfilter & yf, double , double , double);
double errlsstM(double,int, double); 
double errlsstA(lsst & ls,  double);
double errlsstA2(double, double, double); 
double Interpol(double ds, extinc & ex);
int    Extinction(extinc & ex,source & s);
void   getCofactorA(double input[Nx][Nx], double temp[Nx][Nx], int , int , int );
void   getCofactorB(double input[Ny][Ny], double temp[Ny][Ny], int , int , int ); 
double determinantA(double input[Nx][Nx], int); 
double determinantB(double input[Ny][Ny], int); 
double deter55(double AA[5][5]); 
void   inverse(covarian & co, int );  
double RandN(double , double);
double RandR(double , double);
void   sortarray(lsst & ls);
///=========================================================================================

    time_t _timeNow;
    unsigned int _randVal;
    unsigned int _dummyVal;
    FILE * _randStream;
///==============================================================//
///                                                              //
///                  Main program                                //
///                                                              //
///==============================================================//
int main(){

   srand(static_cast<unsigned int>(time(0)));
   
   time(&_timeNow);
   _randStream=fopen("/dev/urandom", "r");
   _dummyVal  =fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   
   time(&_timeNow);
   _randStream=fopen("/dev/urandom", "r");
   _dummyVal  =fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);


   
   VBBinaryLensing vbb;
   vbb.Tol=1.e-4;
   vbb.a1 =0.0;  
   vbb.LoadESPLTable("./ESPL.tbl");


   source   s;
   lens     l;
   CMD      cm;  
   extinc   ex; 
   lsst     ls; 
   yfilter  yf;  
   astromet as;
   covarian co; 

   
   
   for(int i=0; i<N1; i++){cm.Mab_d[i] = new double[M];}
   for(int i=0; i<N2; i++){cm.Mab_b[i] = new double[M];}
   for(int i=0; i<N3; i++){cm.Mab_t[i] = new double[M];}
   for(int i=0; i<N4; i++){cm.Mab_h[i] = new double[M];}
  
  
   FILE* fil0;   FILE* fil1;    FILE* fil2;  FILE* fil2b;  FILE *fpd;    FILE* fildd; 
   FILE* fil3;   FILE* fil3B;  FILE* fil3C;  FILE* fil4;   FILE* fil5;   FILE* filg;   
   FILE* fildf;  FILE* filg2;
   
   char fnLDt[40];  
   char fnEff[40];
   char fnEffB[40];    
   char fnGam[40];
   char fnGamB[40];
   char fnGamC[40];
   char testfil[40];    
   char testfil2[40];    
   char opsim[40];
   char filnam1[40];
   char filnam2[40];
  
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
  
   //    https://scixplorer.org/abs/2019ApJ...873..111I/abstract table (3 Refree wanted)  
   FILE* fil; 
   fil=fopen("./files/sigmaA_LSSTB.txt","r");
   if(!fil){cout<<"cannot read sigmaA_LSSTB.txt: "<<endl; int yye;  cin>>yye;}
   for(int i=0; i<Na;++i){
   fscanf(fil,"%lf   %lf\n",&ls.maga[i], &ls.erra[i]); ///unit of erra is [mas] 
   if(ls.maga[i]<19.0 or ls.maga[i]>24.0 or ls.erra[i]>74.0 or ls.erra[i]<10.0){ 
   cout<<"ERROR SigmaA_LSSTB.txt: "<<ls.maga[i]<<"\t erra: "<<ls.erra[i]<<endl;  int yye; cin>>yye;}}
   cout<<"**** File sigmaA_LSSTB.txt was read ****"<<endl;
   fclose(fil);
    
    
   /*
   fil=fopen("./files/sigma_ELT.txt","r");
   if(!fil){cout<<"cannot read sigma_ELT.txt: "<<endl; int yye;  cin>>yye;}
   for(int i=0; i<Nel;++i){
   fscanf(fil,"%lf   %lf   %lf\n",&ls.mage[i], &ls.ere1[i], &ls.ere2[i] );
   ls.ere2[i]=ls.ere2[i]*0.001;//mili arcsec
   cout<<"mage:  "<<ls.mage[i]<<"\t err_magni(mag): "<<ls.ere1[i]<<"\t err_astro(mili arcsec):  "<<ls.ere2[i]<<endl;}
   cout<<"**** File sigma_ELT.txt was read ****"<<endl;
   fclose(fil);*/
   
   
   
   fil=fopen("./files/metal.txt","r"); 
   if(!fil){cout<<"cannot read metal.txt: "<<endl; int yye;  cin>>yye;}
   for(int i=0; i<met; ++i){
   fscanf(fil,"%lf   %d    %d\n",&yf.Metal[i], &yf.count[i], &yf.number[i]);}
   //cout<<"metal:  "<<yf.Metal[i]<<"\t count:  "<<yf.count[i]<<"\t number "<<yf.number[i]<<endl;
   cout<<"**** File metal.txt was read ****"<<endl;
   fclose(fil); 
    
    
   fil=fopen("./files/yzma.txt", "r"); 
   if(!fil){cout<<"cannot read yzma.txt: "<<endl; int yye;  cin>>yye;}
   for(int i=0; i<YZ; ++i){
   fscanf(fil,"%lf   %lf   %lf  %lf\n",&yf.Age[i],&yf.mm[i],&yf.B[i],&yf.M[i]);}
   cout<<"**** File yzma.txt was read ****"<<endl;
   fclose(fil); 

  
   read_cmd(cm,yf); 
   cout<<"************ read_cmd was done ************"<<endl; 
  
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  
  
   int    save=2000, flagm=0, gg, ss, qq, ww, vv, zz, pp, counter;    
   int    icon=0, nlens=0, hh=0;  
   int    gi, FFG, ndw=0, ndd=0;
   int    flag_det=0, nml=0;
   int    flagf, fi, datf1, datf2;
   int    detectF[6], za[2];
   int    n11, n22, flage, conn, fdetB, fdetP;  
   double errs, errs2, errg, fel, mind, minc, cade;
   double magnio, test, dist, deltaA;
   double magni[M+1], magni0[M+1], Astar0, As1, As0;
   double trajm, trajp, ddf; 
   double chi1, chi2, chi3, chi1a, chi2a, chi3a, sil, sil2; 
   double dchiL, dchiP, dchiA; 
   double flag0, flag1, flag2;
   double vs1, vs2, def1p, def2p, vsave, dt, Mpeak; 
   double ErtE,  ErpiE, ErtetE, Erml, Erdl, Ermul, Ermus, Eru0, Erfb, nsim, Eros;  
   double tE[2], RE[2], piE[2], tetE[2], Vt[2], u0[2], Ml[2], opdP[2], Dl[2], Ds[2], vl[2], vs[2], ros[2], Extv[2], opdB[2];
   double numd[2], Struc[2],  murel[2], vsn[2], DelT[2], mbs[2], fb[2], fwhm[2], Map[2], nbl[2], Ext[2],  catg[2];
   double DtE[2], DRE[2], DpiE[2], DtetE[2], DVt[2], Dros[2], Du0[2], DMl[2], Dopd[2], DDl[2], DDs[2], DExtv[2]; 
   double Dvl[2], Dvs[2], Dmbs[2], Dfb[2], Dfwhm[2], Dvsn[2], DDelT[2], Dmurel[2], DMap[2], Dnbl[2], DExt[2], numB[2], numP[2]; 
   double fraBD[2], fraMS[2], fraBH[2]; 
   double color[4][2], Dcol[4][2], mag[6][2],Dmag[6][2]; 
   double fra[10];
   double mbase, fblend, Gamma,  Neven, EFF, EffiD, EffiL, EffiB, EffiP, nerr=0.0;  
   double dchi1, dchi2, fra1, fra2, fra3, npeak, mi1, mi2, sepa; 
   double Npeak, mag1,  mag2, Sepa, Mi1, Mi2, numdts;  
   double longi, lati, lonn,  Efi; 
   double sig, delt, Lon, Lat; 
   double timeobs, ntot, gtot, nsto, ndto; 
   int nlongi=0, nlongis=0, nlongid=0;
   double shib, effi;  


   
   
   for(int i=0; i<=GG;++i){
   l.NstE[i]=0.0;     l.NsMl[i]=0.0;     l.Nspi[i]=0.0; l.Nsu0[i]=0.0;  l.Nsmb[i]=0.0;  l.Nsfb[i]=0.0;  l.Nsmu[i]=0.0;
   l.NdtE[i]=0.0;     l.NdMl[i]=0.0;     l.Ndpi[i]=0.0; l.Ndu0[i]=0.0;  l.Ndmb[i]=0.0;  l.Ndfb[i]=0.0;  l.Ndmu[i]=0.0;
   l.nste[0][i]=0.0;  l.nste[1][i]=0.0;  l.nste[2][i]=0.0;  
   l.ndte[0][i]=0.0;  l.ndte[1][i]=0.0;  l.ndte[2][i]=0.0; 
   l.tEs[i]=double(tE_min+(tE_max-tE_min)*i/GG);
   l.Mls[i]=double(Ml_min+(Ml_max-Ml_min)*i/GG);
   l.pis[i]=double(pi_min+(pi_max-pi_min)*i/GG);
   l.u0s[i]=double(u0_min+(u0_max-u0_min)*i/GG);
   l.mbs[i]=double(mb_min+(mb_max-mb_min)*i/GG);
   l.fbs[i]=double(fb_min+(fb_max-fb_min)*i/GG);
   l.mus[i]=double(mu_min+(mu_max-mu_min)*i/GG);}
   for(int i=0; i<2; ++i){l.Nhalo[i]=0.0; l.Ndisk[i]=0.0;}
   
  
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  
   
   fil=fopen("./files/EftERubinC.txt","r");
   if(!fil){cout<<"cannot read EftERubinC.txt: "<<endl;  int rre;   cin>>rre;}
   else{
   for(int i=0; i<=GG; ++i){ 
   fscanf(fil,"%lf %lf\n",&l.tEs[i],&l.EFTE[i]);}
   fclose(fil);
   cout<<"File EftERubinC.txt is read "<<endl;  int yye; cin>>yye;}
   
   
   
  
   int IMnum=6;  
   int FlagF=1;  
   
   fil=fopen("./files/MONTGB/files/EfB6.dat","r");
   if(!fil){cout<<"cannot read EfB6.dat: "<<endl; }
   else{
   for(int i=0; i<=GG; ++i){
   fscanf(fil,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %lf %lf %lf %lf\n",
   &l.NdtE[i],&l.NstE[i], &l.NdMl[i],&l.NsMl[i], &l.Ndpi[i],&l.Nspi[i],   &l.Ndu0[i], &l.Nsu0[i], &l.Ndmb[i],&l.Nsmb[i], 
   &l.Ndfb[i],&l.Nsfb[i], &l.Ndmu[i],&l.Nsmu[i], &l.Nhalo[1],&l.Nhalo[0], &l.Ndisk[1],&l.Ndisk[0],
   &l.ndte[0][i],&l.nste[0][i],&l.ndte[1][i],&l.nste[1][i],&l.ndte[2][i],&l.nste[2][i]);}//24
   fclose(fil); 
   cout<<"File EfB6.dat is read "<<endl;  int yye; cin>>yye;}
   
   
   
   if(FlagF==1){
   fil1=fopen("./files/MONTGB/files/MONTsave.dat","a+");//2000 example saved events
   fclose(fil1);}
   
   sprintf(fnEff,"./files/MONTGB/files/%c%c%c%d.dat",'E','f','A',int(IMnum));//efficiency A
   fil2=fopen(fnEff,"w");
   fclose(fil2); 
  
   sprintf(fnEffB,"./files/MONTGB/files/%c%c%c%d.dat",'E','f','B',int(IMnum));//efficiency B
   fil2b=fopen(fnEffB,"w");
   fclose(fil2b); 
  

   sprintf(fnLDt,"./files/MONTGB/files/%c%c%c%c%d.dat",'L','p','a','r', int(IMnum));//Lensing paramaters + errors 
   fil0=fopen(fnLDt,"a+");
   fclose(fil0);
   
   sprintf(fnGam,"./files/MONTGB/files/%c%c%c%c%d.dat",'M','a','p','A', int(IMnum));//Maps
   fil3=fopen(fnGam,"a+");//Maps
   fclose(fil3);

   sprintf(fnGamB,"./files/MONTGB/files/%c%c%c%c%d.dat",'M','a','p','B', int(IMnum));//Maps
   fil3B=fopen(fnGamB,"a+");//Maps
   fclose(fil3B);
   
   sprintf(fnGamC,"./files/MONTGB/files/%c%c%c%c%d.dat",'M','a','p','C', int(IMnum));//Maps
   fil3C=fopen(fnGamC,"a+");//Maps
   fclose(fil3C);
   
   sprintf(testfil,"./%c%c%c%c%d.dat",'t','r','s','t',int(IMnum));//test.dat
   filg=fopen(testfil,"w");
   fclose(filg);  
   
   sprintf(testfil2,"./%c%c%c%c%d.dat",'T','r','s','t',int(IMnum));//test.dat
   filg2=fopen(testfil2,"w");
   fclose(filg2);
   
 
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH   
///                        Monte Carlo Simulation                                                                          ///
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
 
   for(int nl1=550; nl1<801;  nl1++){///801
   
 
   filg2=fopen(testfil2,"w");
   fclose(filg2);  
   nlongi=0; 
   nlongid=0; 
   nlongis=0;
   ntot=0.0; 
   gtot=0.0;   
   
   
   for(int nl2=0; nl2<161;  nl2++){



   hh=int(nl2+nl1*161+1); 
   fildf=fopen("./files/convertcoordinateGB2.txt", "r"); 
   if(!fildf){cout<<"cannot read convertCoordinateGB.txt: "<<endl; int yye;  cin>>yye;}
   fildd=fopen("./files/counterOpsim.txt", "r"); 
   if(!fildd){cout<<"cannot read counterOpsim.txt: "<<endl;   int yye;  cin>>yye;}
   for(int h=1; h<=hh; ++h){
   fscanf(fildf,"%lf   %lf   %lf  %lf\n",&longi,  &lati,  &s.right,  &s.decli);
   fscanf(fildd,"%d    %d    %d\n",&n11,  &n22,  &ndd);}
   fclose(fildf);
   fclose(fildd);  
   
   
   

   lonn= float(-100.0+nl1*ddeg);  
   s.lat=float( -20.0+nl2*ddeg);  
   if(lonn<=0.0)  s.lon=float(360.0+lonn);
   else           s.lon=lonn;
   s.TET=double((360.0-s.lon)*M_PI/180.0);///radian
   s.FI= double(s.lat*M_PI/180.0);///radian
   Disk_model(s,1000); 
   flage=Extinction(ex,s);  
   
   if(n11!= nl1 or n22!=nl2 or fabs(longi-lonn)>ddeg or fabs(lati-s.lat)>ddeg or ndd<0 or flage<1 or s.lon<0.0){
   cout<<"Error n11: "<<n11<<"\t n22:  "<<n22<<"\t nl1: "<<nl1<<"\t nl2: "<<nl2<<"\t flage:  "<<flage<<endl; 
   cout<<"Error Lonn:  "<<lonn<<"\t longi:  "<<longi<<"\t s.lat:  "<<s.lat<<"\t lati:  "<<lati<<endl; 
   int yye; cin>>yye;}
   
   
   
   cout<<">>>>>>>>>> NEW STEP <<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
   cout<<"Nl1:      "<<nl1<<"  \t Nl2:    "<<nl2<<endl;
   cout<<"latitude: "<<s.lat<<"\t longtitude: "<<lonn<<endl;
   cout<<"longi:    "<<longi<<"\t Lati:  "<<lati<<endl;
   cout<<"ndata:    "<<ndd<<"  \t flage:  "<<flage<<endl; 
   cout<<"RA:      "<<s.right<<"\t    DEC:   "<<s.decli<<endl;
   cout<<">>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;

   
   minc=10000.0;  
   if(ndd>=3){
   if(lonn<0.0) sig=-1.0;
   else sig=1.0;
   double delt=fabs(lonn)-floor(fabs(lonn));
   if(delt>1.0 or delt<0.0){cout<<"ERROR longtitude: delt: "<<delt<<"\t lonn: "<<lonn<<endl; int uue; cin>>uue; }
   else if(delt<0.25) Lon=(floor(fabs(lonn))+0.00)*sig;
   else if(delt<0.50) Lon=(floor(fabs(lonn))+0.25)*sig;
   else if(delt<0.75) Lon=(floor(fabs(lonn))+0.50)*sig;
   else               Lon=(floor(fabs(lonn))+0.75)*sig;

   if(s.lat<0.0) sig=-1.0;
   else sig=1.0;
   delt=fabs(s.lat)-floor(fabs(s.lat));
   if(delt>1.0 or delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lat: "<<s.lat<<endl; int uue;  cin>>uue;}
   else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
   else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
   else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
   else                Lat=(floor(fabs(s.lat))+0.75)*sig;
   if(Lon>100.00 or Lon<-100.00 or fabs(Lat)>20.0){
   cout<<"BIG error (stopped program) Lon: "<<Lon<<"\t Lat: "<<Lat<<endl;       int uue;  cin>>uue;}

   sprintf(opsim,"./files/OPSIM/%c%c%c%.2lf%c%.2lf.txt",'o','p','m',Lon,'_',Lat);
   fpd=fopen(opsim,"r");
   if(!fpd){
   cout<<"OPsim_file name: longitude  "<<Lon<<"\t latitude    "<<Lat<<endl;
   cout<<"cannot read opsim.txt: "<<endl;             int yye;  cin>>yye;}
   
   
   
   for(int i=0;  i<Nl;  i++){ ls.tim[i]=-1.0;   ls.filter[i]=-1;   ls.sig5[i]=0.0;}
   for(int i=0;  i<6; ++i)    ls.Nfilt[i]=0;   
   int i1=0;  
   while(!feof(fpd)){
   fscanf(fpd,"%lf    %lf   %lf    %d    %lf\n",&ls.RA,   &ls.DEC,   &ls.tim[i1],   &ls.filter[i1],   &ls.sig5[i1]);
   if(i1>0 and minc>fabs(ls.tim[i1]-ls.tim[i1-1]))  minc=fabs(ls.tim[i1]-ls.tim[i1-1]);///cadence_minimum  
   if(ls.filter[i1]>=6 or ls.filter[i1]<0 or fabs(ls.RA-s.right)>float(FoV*10.1) or fabs(ls.DEC-s.decli)>float(FoV*10.1) or 
   ls.sig5[i1]<=0.0  or ls.tim[i1]<0.0  or ls.tim[i1]>Tobs){ 
   cout<<"ERROR tim[i1]: "<<ls.tim[i1]<<"\t filter: "<<ls.filter[i1]<<"\t i:  "<<i1<<endl;  
   cout<<"ls.RA:   "<<ls.RA<<"\t ls.DEC:  "<<ls.DEC<<"\t sig5:  "<<ls.sig5[i1]<<endl;  
   cout<<"RA:  "<<s.right<<"\t DEC:   "<<s.decli<<endl; int yye; cin>>yye;}
   ls.Nfilt[int(ls.filter[i1])]+=1; 
   i1+=1;}
   if(ndd!=i1 or ls.tim[0]<0.0 or ls.filter[0]<0 or 
   ndd!=int(ls.Nfilt[0]+ls.Nfilt[1]+ls.Nfilt[2]+ls.Nfilt[3]+ls.Nfilt[4]+ls.Nfilt[5])){
   cout<<"Error OPsim_file name: longitude  "<<Lon<<"\t latitude    "<<Lat<<endl;
   cout<<"ndd:  "<<ndd<<"\t i1:  "<<i1<<"\t filter "<<ls.filter[0]<<"\t tim:  "<<ls.tim[0]<<endl;
   cout<<"Nfilt[0]:  "<<ls.Nfilt[0]<<"\t Nfilt[1]:  "<<ls.Nfilt[1]<<"filt[2]:  "<<ls.Nfilt[2]<<endl; 
   cout<<"Nfilt[3]:  "<<ls.Nfilt[3]<<"\t Nfilt[4]:  "<<ls.Nfilt[4]<<"filt[5]:  "<<ls.Nfilt[5]<<endl; 
   cout<<"  File OpSimBaselin.dat was read ****"<<endl;    int yye; cin>>yye;} 
   fclose(fpd);
   }  
  
  
   sortarray(ls);
   if(ls.Nfilt2[int(ls.argm[0])]==0 or ls.Nfilt2[int(ls.argm[1])]==0) ls.flag=0; 
   else   ls.flag=1;
   filg=fopen(testfil,"w");
   fclose(filg);  
   icon=0; 
   nlens=0; 
   nsim=0.0;   
   nerr=0.0;  
   for(int i=0; i<Num; ++i){s.nssim[i]=0.0;  s.nsdet[i]=0.0;}
   for(int i=0; i<=GG;++i) { l.nstE[i]=0.0;   l.ndtE[i]=0.0;}  
 
 
 
///############################################################################################################################
///            Loop of lensing for a direction  
///############################################################################################################################
         


do{

   nsim+=1.0;
   nlongis+=1;
   //do{
   func_source(s,cm, ex);
   func_lens(l,s);
   //}while(l.tE>900.0);
   
   optical_depth(s);
   s.nssim[s.nums]+=1.0;  
   l.A0=vbb.ESPLMag2(l.u0,s.ros);

  
   fdetB=0;  
   fdetP=0;  
   for(int i=0; i<M; ++i){  
   detectF[i]=0;  
   s.Mpeak[i]=double(s.magb[i]-2.5*log10(s.blend[i]*l.A0+1.0-s.blend[i]));  
   if(RandR(0.0,1.0)<=s.blend[i]){
   if(s.magb[i]<=thre[i]){fdetB+=1;  detectF[i]=1;} 
   if(s.Mpeak[i]<=thre[i])fdetP+=1;}}
   
     
   if(l.A0<1.0 or s.Mpeak[0]>s.magb[0] or s.Mpeak[1]>s.magb[1] or s.Mpeak[2]>s.magb[2] or fdetP<fdetB){
   cout<<"Error_0 magb[0]:  "<<s.magb[0]<<"\t magb[1]:  "<<s.magb[1]<<"\t magb[2]:  "<<s.magb[2]<<endl;
   cout<<"thre[0]:  "<<thre[0]<<"\t  thre[1]:  "<<thre[1]<<"\t thre[2]:  "<<thre[2]<<endl;
   cout<<"l.A0:  "<<l.A0<<"\t blend: "<<s.blend[2]<<"\t fdetB: "<<fdetB<<"\t fdetP"<<fdetP<<endl;
   int uue; cin>>uue;}
   
   



   s.errM=0.0;     s.errA=0.0;   vsave=0.0;  npeak=0.0;
   chi1=0.0;       chi2=0.0;     chi3=0.0;   ndw=0;
   chi1a=0.0;      chi2a=0.0;    chi3a=0.0;  Mi2=0.0; 
   gi=0;           sepa=0.0;     Mi1=0.0;    FFG=0;
   dchiL=0.0;      dchiP=0.0;    dchiA=0.0;  flagf=0;   
   flag0=0.0;      flag1=0.0;    flag2=0.0;  flagm=0; 
   def1p=0.0;      s.def1c=0.0;  dt=3.0;///days 
   def2p=0.0;      s.def2c=0.0;  flag_det=0;
   for(int i=0; i<coun; ++i){
   l.timn[i]=0.0;    l.magn[i]=0.0;   l.soux[i]=0.0;   l.souy[i]=0.0;  
   l.errm[i]=0.0;    l.erra[i]=0.0;   l.tele[i]=-1;}
   for(int i=0;  i<nq;  ++i)      co.resu[i]=-1.0;   
   timeobs=double(fabs(ls.tim[ndd-1]-ls.tim[0])); 
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   if((fdetP>=2 or fdetB>=2) and ndd>=2){///Visible source stars
   s.nsdet[s.nums]+=1.0;
   flagf=1; 
   icon+=flagf;  
   nlongid+=1;
   flagm=0; 
   if(RandR(0.0,100.0)<1.0 and save<2000 and FlagF>0 and ndd>=10){
   save+=1;  
   sprintf(filnam1,"./files/MONTGB/files/1/%c%c%c%c%d.dat",'m','a','g','C', int(save));
   sprintf(filnam2,"./files/MONTGB/files/1/%c%c%c%c%d.dat",'d','a','t','C', int(save));
   fil4=fopen(filnam1,"w");    
   fil5=fopen(filnam2,"w");
   flagm=1;}    
    
    
   for(double tim=float(0.0*year-20.0);  tim<float(10.0*year+20.0); tim=tim+dt){
   lightcurve(s,l,as,tim);
   if(s.ut0<double(30.0*s.ros))    Astar0= vbb.ESPLMag2(s.ut0, s.ros); 
   else                            Astar0= double(s.ut0*s.ut0+2.0)/sqrt(s.ut0*s.ut0*(s.ut0*s.ut0+4.0));
   if(s.ut<double(30.0*s.ros))     s.Astar=vbb.ESPLMag2(s.ut, s.ros);
   else                            s.Astar=double(s.ut*s.ut+2.0)/sqrt(s.ut*s.ut*(s.ut*s.ut+4.0));
   As0    =double( Astar0*s.blend[2] + 1.0 - s.blend[2]);
   As1    =double(s.Astar*s.blend[2] + 1.0 - s.blend[2]);//LSST r-band 
   vs1    =double(s.mus1+(s.def1c-def1p)/dt);//[mas/days]  
   vs2    =double(s.mus2+(s.def2c-def2p)/dt);//[mas/days] 
   def1p  =s.def1c;//pervious 
   def2p  =s.def2c;
   
   if(flagm>0){
   fprintf(fil4,"%.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n",
   tim,s.def1c,s.def2c,As0,As1,s.pos1b,s.pos2b,s.pos1c,s.pos2c,s.def1a,s.def2a,l.pos1,l.pos2);}
   for(int i=0; i<int(M+1); ++i){
   magni0[i] =s.magb[i]-2.5*log10( Astar0*s.blend[i] + 1.0 - s.blend[i]);
   magni[i]  =s.magb[i]-2.5*log10(s.Astar*s.blend[i] + 1.0 - s.blend[i]);}
   

   if(tim>=0.0 and tim<=Tobs and tim>=ls.tim[0] and tim<=ls.tim[int(ndd-1)] and gi<ndd and gi>=0 and tim>=ls.tim[gi]){
   fi=int(ls.filter[gi]);
   
   
   if(magni[fi]>=satu[fi] and magni[fi]<=thre[fi]){  
   errg=errlsstM(magni[fi], int(fi), double(ls.sig5[gi]));//[mag]
   errs=errlsstA(ls, magni[2]); ///[mas]
   deltaA=fabs(pow(10.0,-0.4*errg)-1.0)*(s.blend[fi]*s.Astar+1.0-s.blend[fi]);
   magnio=magni[fi]+RandN(errg,3.0);
   
   chi1 +=fabs( (magnio-  magni[fi])*(magnio-  magni[fi])/(errg*errg)); ///real 
   chi2 +=fabs( (magnio- magni0[fi])*(magnio- magni0[fi])/(errg*errg)); ///real no parallax 
   chi3 +=fabs( (magnio- s.magb[fi])*(magnio- s.magb[fi])/(errg*errg)); ///baseline
   
   sil= RandN(errs*sqrt(2.0), 3.0);
   sil2=RandN(errs,3.0);
   
   trajm=sqrt(s.pos1b*s.pos1b + s.pos2b*s.pos2b);///stright+parallax   
   trajp=sqrt(s.pos1c*s.pos1c + s.pos2c*s.pos2c);///stright+parallax+lensing  
   
   chi1a+= fabs( (trajp+sil-trajp)*(trajp+sil-trajp)/(errs*errs*2.0) ); ///real trajectory
   chi2a+= fabs( (trajp+sil-trajm)*(trajp+sil-trajm)/(errs*errs*2.0) ); ///real without lensing 
   //chi3a+= fabs( sil2*sil2/(errs*errs) ); 
     
   l.timn[ndw]=tim;
   l.magn[ndw]=magni[2];///scale to r-LSST band, but the errors are different 
   l.errm[ndw]=errg;
   l.soux[ndw]=s.pos1c;
   l.souy[ndw]=s.pos2c;
   l.erra[ndw]=errs;
   l.tele[ndw]=0; 
   
   
   
   
   mi1 =double(s.magb[fi]-2.5*log10(fabs((s.Astar+1.0)*0.5*s.blend[fi] + 1.0 - s.blend[fi])));//image_1
   mi2 =double(s.magb[fi]-2.5*log10(fabs((s.Astar-1.0)*0.5*s.blend[fi] + 1.0 - s.blend[fi])));//image_2   
   errs2= fabs(errlsstA2(fabs(mi1-mi2), errs, errg)); 
   if(mi1>=satu[fi] and mi1<=thre[fi] and mi2>=satu[fi] and mi2<=thre[fi] and double(sqrt(s.ut*s.ut+4.0)*l.tetE/errs2)>=1.0){
   npeak+=1.0; 
   Mi1 +=double(s.magb[2]-2.5*log10(fabs((s.Astar+1.0)*0.5*s.blend[2] + 1.0 - s.blend[2])));
   Mi2 +=double(s.magb[2]-2.5*log10(fabs((s.Astar-1.0)*0.5*s.blend[2] + 1.0 - s.blend[2])));
   sepa+=double(sqrt(s.ut*s.ut+4.0)*l.tetE/errs2);}
   
   flag2=0.0;
   if(fabs(magnio-s.magb[fi])>fabs(3.0*errg))   flag2=1.0;
   if(ndw>2 and float(flag0+flag1+flag2)>2.0)   flag_det=1;
   flag0=flag1;
   flag1=flag2;
   
   if(flagm>0){
   fprintf(fil5,"%.4lf %.6lf %.6lf %.4lf %.4lf %.4lf %.4lf %.6lf %d\n",tim,As1+RandN(deltaA,3.0),deltaA,s.pos1c+RandN(errs,3.0),
   s.pos2c+RandN(errs,3.0),s.def1c+RandN(errs,3.0),s.def2c+RandN(errs,3.0),errs,fi);}
   
   if(gi<0 or gi>int(Nl-1) or fi<0 or fi>5 or ndw>ndd or errs<0.0 or errs2<0.0 or errg<0.0 or deltaA<0.0 or gi>ndd){
   cout<<"Error(LSST):  filter:  "<<fi<<"\t counter_filter:  "<<gi<<"\t errs2:  "<<errs2<<endl;
   cout<<"Err_Astrometry:  "<<errs<<"\t Err_magnitude:  "<<errg<<"\t Err_Magnification "<<deltaA<<endl;
   cout<<"ndd:  "<<ndd<<"\t gi:  "<<gi<<"\t ndw:  "<<ndw<<endl;  int yys; cin>>yys;}
   
   s.errM+=deltaA; 
   s.errA+=errs; 
   vsave+=sqrt(vs1*vs1+vs2*vs2); 
   ndw+=1;
   }//magnitude limit
   gi+=1;}
   if(tim<-50.0 or tim>(10.0*year+50.0))             dt=10.0; //days
   else if(tim<ls.tim[0] or tim>ls.tim[int(ndd-1)] ) dt=3.0;///days
   else{
   if(gi>0 and gi<ndd) cade=float(ls.tim[gi]-ls.tim[gi-1]);///days
   else                cade=minc; 
   if(flagm>0)         dt=double(cade*0.1); 
   else                dt=double(cade);}}
   

   Mi1   =double(Mi1/(npeak+eps));
   Mi2   =double(Mi2/(npeak+eps));
   sepa  =double(sepa/(npeak+eps));
   vsave =double( vsave/(ndw+eps));
   s.errM=double(s.errM/(ndw+eps)); 
   s.errA=double(s.errA/(ndw+eps)); 
   dchiL =fabs(chi3- chi1);//lensing+parallax_effect
   dchiP =fabs(chi2- chi1);//parallax_effect 
   dchiA =fabs(chi2a-chi1a);//deflection_effect
   
   
   if(s.FWHM<=timeobs and  dchiL>=float(400.0) and flag_det>0 and ndd>=10){//lensing
   FFG=1;//Lensing
   nlens+=1;
   
   
   FisherM(s, l, as, co, ndw, vbb); 
   if(co.flagi>0){
   nerr+=1.0;  
   ErrorCal(co,l,s);
   if(nerr<20){
   fil0=fopen(fnLDt,"a+");//>>>>>>>>>>>> Lpar1.dat   <<<<<<<<<<<<<<<   
   fprintf(fil0,"%d %d %d %.4lf %.5lf %.4lf %.4lf %d %.4lf %.4lf %.4lf %.4lf %.4lf %.5lf %.5lf %.6lf %.6lf %.6lf %.5lf %.7lf %.4lf %.4lf %.4lf %.4lf %.4lf  %.4lf  %.4lf  %.4lf\n",
   nl1,nl2,ndd,l.Dl,l.Ml,l.mul,l.u0,l.catg,s.Ds,s.magb[1],s.magb[2],s.magb[3],s.magb[4],s.fb[0],l.tE,l.piE,l.tetE,s.ros,l.murel,
   s.opt*1.0e6,co.resu[1],co.resu[2],co.resu[3],co.resu[5],co.resu[9],co.resu[10],co.resu[13],co.resu[15]);//28
   fclose(fil0);}}
   
   
   
   filg2=fopen(testfil2,"a+");////TEST.txt a temporary file for saving for each lonitude value 
   fprintf(filg2,
   "%d %d %.8lf %.6lf %.8lf %.7lf %.6lf %.6lf %.6lf %.7lf %.6lf %.6lf %.6lf %.6lf %.6lf %.7lf %d  %d %.8lf  "//19
   "%.6lf %.6lf  %.6lf  %.8lf  %.8lf %.8lf %.8lf %.8lf  "//27
   "%.8lf  %.8lf  %.8lf  %.8lf  %.7lf  %.2lf %d  %.6lf %.4lf %.4lf "//37
   "%.2lf %.8lf %.8lf %.8lf %.7lf   %.6lf  %.8lf  %d  %d  %d %.5lf  %.5lf  %.5lf %.5lf  %.5lf  %.5lf\n",//53
   icon,FFG,l.tE,l.RE/AU,l.piE,l.tetE,l.Vt,l.u0,l.Ml,s.opt*1.0e6,l.Dl,s.Ds,l.vl,s.vs,s.mbs[1],s.fb[0],gg,l.struc,s.FWHM,//19
   vsave/(s.mus+eps),l.DeltaT/(s.errA+eps),l.murel*year,co.resu[0],co.resu[1],co.resu[2],co.resu[3],co.resu[5],//27
   co.resu[9],co.resu[10],co.resu[13],co.resu[14],s.Map[2],s.nsbl[2],co.flagi,s.Ai[2],dchiP/(ndw+eps),dchiA/(ndw+eps),//37
   npeak,Mi1,Mi2,sepa,s.ros,s.Av, co.resu[15], l.catg, fdetB, fdetP,s.magb[0],s.magb[1],s.magb[2],s.magb[3],s.magb[4],s.magb[5]);//53
   nlongi +=1; 
   fclose(filg2);}
   
   
   gg=int(FunctE(l));
   ss=int(FuncMl(l));
   qq=int(FuncPi(l));
   ww=int(Funcu0(l));
   vv=int(FuncMu(l));
   zz=int(FuncMb(l, s.magb[2]));
   pp=int(FuncFb(l,s.blend[2])); 
   if(l.struc==4 or l.struc==3) l.Nhalo[0]+=1.0;///stellar and dark-matter Halos
   else                         l.Ndisk[0]+=1.0;//Bulge and thin and thick disk
   l.nstE[gg]+=1.0;  l.NsMl[ss]+=1.0;  l.Nspi[qq]+=1.0;  l.Nsu0[ww]+=1.0; l.Nsmb[zz]+=1.0; l.Nsfb[pp]+=1.0; l.Nsmu[vv]+=1.0;
   if(FFG>0){
   if(l.struc==4 or l.struc==3) l.Nhalo[1]+=1.0;  
   else                         l.Ndisk[1]+=1.0;
   l.ndtE[gg]+=1.0;  l.NdMl[ss]+=1.0;  l.Ndpi[qq]+=1.0;  l.Ndu0[ww]+=1.0; l.Ndmb[zz]+=1.0; l.Ndfb[pp]+=1.0; l.Ndmu[vv]+=1.0;}
   
   
   filg=fopen(testfil,"a+");
   fprintf(filg,
   "%d %d %.8lf %.6lf %.8lf %.7lf %.6lf %.6lf %.6lf %.7lf %.6lf %.6lf %.6lf %.6lf %.6lf %.7lf %d  %d %.8lf  "//19
   "%.6lf %.6lf  %.6lf  %.8lf  %.8lf %.8lf %.8lf %.8lf  "//27
   "%.8lf  %.8lf  %.8lf  %.8lf  %.7lf  %.2lf %d  %.6lf %.4lf %.4lf "//37
   "%.2lf %.8lf %.8lf %.8lf %.7lf   %.6lf  %.8lf  %d  %d  %d %.5lf  %.5lf  %.5lf %.5lf  %.5lf  %.5lf\n",//53
   icon,FFG,l.tE,l.RE/AU,l.piE,l.tetE,l.Vt,l.u0,l.Ml,s.opt*1.0e6,l.Dl,s.Ds,l.vl,s.vs,s.mbs[1],s.fb[0],gg,l.struc,s.FWHM,//19
   vsave/(s.mus+eps),l.DeltaT/(s.errA+eps),l.murel*year,co.resu[0],co.resu[1],co.resu[2],co.resu[3],co.resu[5],//27
   co.resu[9],co.resu[10],co.resu[13],co.resu[14],s.Map[2],s.nsbl[2],co.flagi,s.Ai[2],dchiP/(ndw+eps),dchiA/(ndw+eps),//37
   npeak,Mi1,Mi2,sepa,s.ros,s.Av, co.resu[15], l.catg, fdetB, fdetP,s.magb[0],s.magb[1],s.magb[2],s.magb[3],s.magb[4],s.magb[5]);//53
   fclose(filg);
   
   
   
   if(flagm>0){
   fil1=fopen("./files/MONTGB/files/MONTsave.dat","a+");
   fprintf(fil1,
   "%d %d %d  %d  %.5lf   %.5lf   %.5lf  %.5lf  "//7
   "%d   %.5lf   %.4lf   %.4lf  %d  %.5lf   %.5lf   %.4lf  %.4lf   %d  %.4lf   %.1lf  " //19
   "%.6lf  %.6lf  %.5lf  %.5lf   %.5lf  %.6lf   %.7lf  %.7lf    %.7lf   %.7lf   %.6lf  "  //30
   "%.4lf  %.4lf  %.4lf  %.4lf   %.4lf  %.4lf   %.4lf  %.4lf  %.4lf  %.4lf   %.4lf  %.4lf  " //42
   "%.4lf  %.4lf  %.4lf  %.4lf   %.4lf  %.4lf   %.5lf  %.5lf  %.5lf  %.5lf   %.5lf  %.5lf  " //54
   "%.5lf  %.5lf  %.5lf  %.5lf   %.5lf  %.5lf   %.3lf  %.3lf  %.3lf  %.3lf   %.3lf  %.3lf  " //66
   "%d  %.3lf  %.7lf   %.7lf  %.7lf  %.7lf   %.7lf  %.6lf %d %d %d %d %d %d %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf "
   "%d  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.3lf  %.3lf  %.1lf %.5lf  %.5lf  %.5lf  %.5lf\n",
   save, FFG, nl1, nl2, s.lat, s.lon, s.right, s.decli,//7
   l.struc, l.Ml, l.Dl, l.vl, s.struc, s.mass,s.Ds, s.vs, s.logT, s.cl, s.Rstar, s.age,//19
   l.tE, l.RE/AU, l.t0, l.murel, l.Vt, l.u0, s.opt*1.0e6,  l.tetE, l.piE, s.ros, s.FWHM,//30
   s.Mab[0],s.Mab[1],s.Mab[2],s.Mab[3],s.Mab[4], s.Mab[5],s.Map[0],s.Map[1],s.Map[2],s.Map[3],s.Map[4], s.Map[5],//42
   s.magb[0],s.magb[1],s.magb[2],s.magb[3],s.magb[4],s.magb[5],s.blend[0],s.blend[1],s.blend[2],s.blend[3],s.blend[4],s.blend[5],//54  
   s.Ai[0], s.Ai[1], s.Ai[2], s.Ai[3],s.Ai[4],s.Ai[5],s.nsbl[0],s.nsbl[1],s.nsbl[2],s.nsbl[3],s.nsbl[4], s.nsbl[5], //66
   ndw,dchiL/(ndw+eps),s.Rostart[0],s.Nstart,l.deltao*180.0/M_PI,l.murel1,l.murel2, 
   s.xi*180.0/M_PI,detectF[0],detectF[1],detectF[2],detectF[3],
   detectF[4],detectF[5],co.resu[1],co.resu[2],co.resu[3],co.resu[5],co.resu[9],co.resu[10],co.resu[13],co.resu[14],co.resu[15],//89
   l.catg, s.errM, s.errA, Mi1, Mi2, sepa, dchiP/(ndw+eps), dchiA/(ndw+eps), npeak, s.mus1, s.mus2, l.mul1, l.mul2);//98+4
   fclose(fil1);
   fclose(fil4);  
   fclose(fil5);}}//if star is visible
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH   
   conn=-1;  
   if(ndd>=20 and nlens>=120 and nerr>35)     conn=1; //lensing may happen  
   if(ndd>=10 and ndd<20   and nlens>=45)     conn=1; //lensing may happen  
   if(ndd<10 and ndd>=3    and icon>=120)     conn=1; //detectin may happen
   if(ndd<3)                                  conn=1; //No lensing no detection
}while(conn<0); //end of loop icon
   
   





///############################################################################################################################
///            Indicating & Saving the Efficiency of Lensing 
///############################################################################################################################
   fil2 =fopen( fnEff , "w"); 
   fil2b=fopen(fnEffB , "w");       
   for(int i=0; i<=GG; ++i){
   l.NstE[i]+= l.nstE[i];  
   l.NdtE[i]+= l.ndtE[i]; 
   if(l.catg==1){l.nste[0][i]+=l.nstE[i];   l.ndte[0][i]+=l.ndtE[i]; }//BD
   if(l.catg==2){l.nste[1][i]+=l.nstE[i];   l.ndte[1][i]+=l.ndtE[i]; }//MS
   if(l.catg==3){l.nste[2][i]+=l.nstE[i];   l.ndte[2][i]+=l.ndtE[i]; }//BH
   l.ndtE[i]=double(l.ndtE[i]/(l.nstE[i]+eps));//[0.0,1.0]
   if(l.ndtE[i]==0.0 and l.nstE[i]==0.0 and i>0){
   nsto=0.0;  ndto=0.0; 
   for(int j=0; j<i; ++j){
   nsto+=l.nstE[j];
   ndto+=l.ndtE[j];}
   l.ndtE[i]=double(ndto/(nsto+eps));}
   
   fprintf(fil2,
   "%.4lf %.6lf  %.4lf %.6lf %.4lf %.6lf  %.4lf %.6lf  %.4lf %.6lf  %.5lf %.5lf  %.4lf %.6lf %.6lf %.6lf  %.6lf  %.3lf %.3lf\n",
   l.Mls[i],double(l.NdMl[i]*100.0/(l.NsMl[i]+eps)), 
   l.pis[i],double(l.Ndpi[i]*100.0/(l.Nspi[i]+eps)), 
   l.u0s[i],double(l.Ndu0[i]*100.0/(l.Nsu0[i]+eps)), 
   l.mbs[i],double(l.Ndmb[i]*100.0/(l.Nsmb[i]+eps)), 
   l.fbs[i],double(l.Ndfb[i]*100.0/(l.Nsfb[i]+eps)),
   l.mus[i],double(l.Ndmu[i]*100.0/(l.Nsmu[i]+eps)), 
   l.tEs[i],double(l.NdtE[i]*100.0/(l.NstE[i]+eps)), double(l.ndte[0][i]*100.0/(l.nste[0][i]+eps)),
   double(l.ndte[1][i]*100.0/(l.nste[1][i]+eps)),    double(l.ndte[2][i]*100.0/(l.nste[2][i]+eps)), 
   double(l.Nhalo[1]*100.0/l.Nhalo[0]), double(l.Ndisk[1]*100.0/l.Ndisk[0]) ); //19
   
   fprintf(fil2b,"%.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf\n",//24
   l.NdtE[i],l.NstE[i], l.NdMl[i],l.NsMl[i],  l.Ndpi[i], l.Nspi[i],   l.Ndu0[i], l.Nsu0[i], l.Ndmb[i],l.Nsmb[i], 
   l.Ndfb[i],l.Nsfb[i], l.Ndmu[i],l.Nsmu[i], l.Nhalo[1], l.Nhalo[0],  l.Ndisk[1],l.Ndisk[0],
   l.ndte[0][i],l.nste[0][i], l.ndte[1][i],l.nste[1][i], l.ndte[2][i],l.nste[2][i]);}//24
   fclose(fil2);  
   fclose(fil2b);  
  

  
   s.nstart=0.0;  
   for(int i=0; i<Num; ++i){
   test=double(s.nsdet[i]/(s.nssim[i]+eps)); 
   if(s.nsdet[i]==0.0 and s.nssim[i]==0.0 and i>0){
   nsto=0.0;  ndto=0.0; 
   for(int j=0; j<i; ++j){
   nsto+=s.nssim[j];
   ndto+=s.nsdet[j];}
   test=double(ndto/(nsto+eps));}
   s.nstart+=s.Rostari[i][0]*(s.Nstart/s.Rostart[0])*test;}//Number/deg^{2}
   
   
   

///############################################################################################################################
///                Indicating the average of Lensing parameters 
///############################################################################################################################

   for(int i=0;i<2; ++i){
   fraBD[i]=0.0; fraMS[i]=0.0;  fraBH[i]=0.0; 
   tE[i]=0.0,    RE[i]=0.0,     piE[i]=0.0,  tetE[i]=0.0;  Vt[i]=0.0;     ros[i]=0.0;  Extv[i]=0.0; 
   u0[i]=0.0,    Ml[i]=0.0,     opdB[i]=0.0; opdP[i]=0.0;  Dl[i]=0.0,     Ds[i]=0.0;  catg[i]=0.0;  
   vl[i]=0.0,    vs[i]=0.0,     mbs[i]=0.0,  fb[i]=0.0;    numd[i]=0.0;   numB[i]=0.0;  numP[i]=0.0;  
   fwhm[i]=0.0;  vsn[i]=0.0;    DelT[i]=0.0; Struc[i]=0.0; murel[i]=0.0; 
   Map[i]=0.0;   nbl[i]=0.0;    Ext[i]=0.0;
   for(int j=0;j<4;++j){color[j][i]=0.0;     Dcol[j][i]=0.0;  }
   for(int j=0;j<6;++j){  mag[j][i]=0.0;     Dmag[j][i]=0.0;  }
   DtE[i]=0.0,   DRE[i]=0.0,    DpiE[i]=0.0, DtetE[i]=0.0, DVt[i]=0.0,   Dros[i]=0.0; 
   Du0[i]=0.0,   DMl[i]=0.0,    Dopd[i]=0.0, DDl[i]=0.0,   DDs[i]=0.0,   DExtv[i]=0.0; 
   Dvl[i]=0.0,   Dvs[i]=0.0,    Dmbs[i]=0.0, Dfb[i]=0.0,   Dfwhm[i]=0.0, Dvsn[i]=0.0; 
   DDelT[i]=0.0, Dmurel[i]=0.0, DMap[i]=0.0, Dnbl[i]=0.0,  DExt[i]=0.0;}
   for(int i=0;i<10; ++i) fra[i]=0.0; 
   EffiL=0.0;    EFF=0.0;    Gamma=0.0;   Neven=0.0; EffiD=0.0; EffiB=0.0;  EffiP=0.0;     
   ErtE=0.0;     ErpiE=0.0;  ErtetE=0.0;  Erml=0.0; Erfb=0.0;  nerr=0.0;
   Erdl=0.0;     Ermul=0.0;  Ermus=0.0;   Eru0=0.0; Eros=0.0;   
   dchi1=0.0;    dchi2=0.0;  fra1=0.0;    fra2=0.0;  fra3=0.0;   
   Npeak=0.0;    mag1=0.0;   mag2=0.0;    Sepa=0.0;  numdts=0.0;  
   
   
if(icon>=1 and nlens>=2){
   filg=fopen(testfil,"r");
   for(int i=0; i<icon; ++i){
   fscanf(filg,
   "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
   "%lf %lf %lf %lf %lf %d %lf %lf %lf  %lf %lf %lf %lf %lf %lf  %lf %d  %d  %d %lf %lf %lf %lf %lf %lf\n",//53
   &counter, &FFG, &l.tE, &l.RE, &l.piE, &l.tetE, &l.Vt, &l.u0, &l.Ml,    
   &s.opt, &l.Dl, &s.Ds, &l.vl, &s.vs, &mbase, &fblend, &gg, &l.struc, &s.FWHM, &vsave, &l.DeltaT, &l.murel,
   &co.resu[0],  &co.resu[1], &co.resu[2], &co.resu[3], &co.resu[5], &co.resu[9], &co.resu[10], &co.resu[13],
   &co.resu[14], &s.Map[2], &s.nsbl[2], &co.flagi, &s.Ai[2], &dchiP, &dchiA, &npeak, &Mi1, &Mi2, &sepa, &s.ros, &s.Av, &co.resu[15],
   &l.catg, &fdetB, &fdetP, &s.magb[0],&s.magb[1],&s.magb[2],&s.magb[3],&s.magb[4],&s.magb[5]);//53
   
   
   
   if(l.struc==4 or l.struc==3)          Struc[0]+=1.0;//Halo-lensing fraction 
   tE[0]+=l.tE;       RE[0]+=l.RE;       piE[0]+=l.piE;     tetE[0]+=l.tetE;  Vt[0]+=l.Vt;    ros[0]+=s.ros; 
   u0[0]+=l.u0;       Ml[0]+=l.Ml;       Dl[0]+= l.Dl;      Ds[0]+=s.Ds;    Extv[0]+=s.Av;
   vl[0]+=l.vl;       vs[0]+=s.vs;       mbs[0]+=mbase;     fb[0]+=fblend;  numd[0]+=1.0;
   fwhm[0]+=s.FWHM;   vsn[0]+=vsave;     DelT[0]+=l.DeltaT; murel[0]+=l.murel; 
   Map[0]+=s.Map[2];  nbl[0]+=s.nsbl[2]; Ext[0]+=s.Ai[2];   catg[0]+=l.catg*1.0;
   if(fdetB>=2)       {opdB[0]+=s.opt;  numB[0]+=1.0; }
   if(fdetP>=2)       {opdP[0]+=s.opt;  numP[0]+=1.0; }
   if(l.catg==1)      fraBD[0]+=1.0; 
   if(l.catg==2)      fraMS[0]+=1.0;  
   if(l.catg==3)      fraBH[0]+=1.0; 
   color[0][0]+=float(s.magb[0]-s.magb[2]);///u-r
   color[1][0]+=float(s.magb[1]-s.magb[3]);///g-i
   color[2][0]+=float(s.magb[2]-s.magb[4]);///r-z   
   color[3][0]+=float(s.magb[ls.argm[0]]-s.magb[ls.argm[1]]);///r-z   
   
   mag[0][0]  +=float(s.magb[0]);//g-band
   mag[1][0]  +=float(s.magb[1]);//i-band
   mag[2][0]  +=float(s.magb[2]);//z-band
   mag[3][0]  +=float(s.magb[3]);//z-band
   mag[4][0]  +=float(s.magb[4]);//z-band
   mag[5][0]  +=float(s.magb[5]);//z-band
   
   if(FFG>0){
   if(l.struc==4 or l.struc==3)          Struc[1]+=1.0;//Halo-lensing fraction 
   tE[1]+=l.tE;       RE[1]+=l.RE;       piE[1]+=l.piE;     tetE[1]+=l.tetE;   Vt[1]+=l.Vt;  ros[1]+=s.ros; 
   u0[1]+=l.u0;       Ml[1]+=l.Ml;       Dl[1]+= l.Dl;      Ds[1]+=s.Ds;     Extv[1]+=s.Av;
   vl[1]+=l.vl;       vs[1]+=s.vs;       mbs[1]+=mbase;     fb[1]+=fblend;   numd[1]+=1.0;
   fwhm[1]+=s.FWHM;   vsn[1]+=vsave;     DelT[1]+=l.DeltaT; murel[1]+=l.murel;
   Map[1]+=s.Map[2];  nbl[1]+=s.nsbl[2]; Ext[1]+=s.Ai[2];   catg[1]+=l.catg*1.0;
   if(fdetB>=2)       {opdB[1]+=s.opt;   numB[1]+=1.0; }
   if(fdetP>=2)       {opdP[1]+=s.opt;   numP[1]+=1.0; }
   if(l.catg==1)      fraBD[1]+=1.0; 
   if(l.catg==2)      fraMS[1]+=1.0;  
   if(l.catg==3)      fraBH[1]+=1.0; 
   color[0][1]+=float(s.magb[0]-s.magb[2]);///u-r 
   color[1][1]+=float(s.magb[1]-s.magb[3]);///g-i
   color[2][1]+=float(s.magb[2]-s.magb[4]);///r-z  
   color[3][1]+=float(s.magb[ls.argm[0]]-s.magb[ls.argm[1]]);///r-z   
   mag[0][1]  +=float(s.magb[0]);//g-band
   mag[1][1]  +=float(s.magb[1]);//i-band
   mag[2][1]  +=float(s.magb[2]);//z-band
   mag[3][1]  +=float(s.magb[3]);//g-band
   mag[4][1]  +=float(s.magb[4]);//i-band
   mag[5][1]  +=float(s.magb[5]);//z-band
   
   dchi1+=dchiP;      
   dchi2+=dchiA; 
   Npeak+=npeak;
   if(Mi1>0.0 and Mi2>0.0 and sepa>0.0 and npeak>0.0){ 
   numdts+=1.0;  
   mag1 +=Mi1;   
   mag2 +=Mi2;    
   Sepa +=sepa;}
   if(dchiP>=2.0) fra1+=1.0;///parallax detected >2Ndata
   if(dchiA>=2.0) fra2+=1.0;///astrometric trajectory detected >3Ndata
   if(npeak>=2.0) fra3+=1.0;///two images are resolved in at least 3 data points 
   
   ///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   effi=-1.0; shib=0.0; 
   if(l.tE<=l.tEs[0]){
   shib=double(l.EFTE[1]-l.EFTE[0])/(l.tEs[1]-l.tEs[0]);
   effi=double(l.EFTE[0]+shib*(l.tE-l.tEs[0]))*0.01;
   if(effi<0.0) effi=0.0;}
   
   else if(l.tE>=l.tEs[GG]){
   shib=double(l.EFTE[GG]-l.EFTE[GG-1])/(l.tEs[GG]-l.tEs[GG-1]);
   effi=double(l.EFTE[GG]+shib*(l.tE-l.tEs[GG]))*0.01;
   if(effi<0.0) effi=0.0;}
   
   else{
   for(int ij=1; ij<=GG; ++ij){
   if(double((l.tE-l.tEs[ij])*(l.tE-l.tEs[ij-1]))<0.0 or l.tE==l.tEs[ij-1]){
   shib=double(l.EFTE[ij]-l.EFTE[ij-1])/(l.tEs[ij]-l.tEs[ij-1]);
   effi=double(l.EFTE[ij-1]+shib*(l.tE-l.tEs[ij-1]))*0.01;
   if(effi<0.0) effi=0.0; ij=GG+5;}}}
   
   if(effi<0.0 or effi>1.0){
   cout<<"There is a problem:  "<<effi<<"\t tE:  "<<l.tE<<endl; 
   int uue;    cin>>uue;}
   ///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   
   EFF+=double(effi/l.tE);//1/days
   if(co.flagi>0){
   if(co.resu[0]>0.0 and co.resu[0]<0.5)   fra[0]+=1.0;//u0   
   if(co.resu[1]>0.0 and co.resu[1]<0.5)   fra[1]+=1.0;//tE
   if(co.resu[2]>0.0 and co.resu[2]<0.5)   fra[2]+=1.0;//fb
   if(co.resu[3]>0.0 and co.resu[3]<0.5)   fra[3]+=1.0;//piE
   if(co.resu[5]>0.0 and co.resu[5]<0.5)   fra[4]+=1.0;//tetE
   if(co.resu[9]>0.0 and co.resu[9]<0.5)   fra[5]+=1.0;//ML
   if(co.resu[10]>0.0 and co.resu[10]<0.5) fra[6]+=1.0;//Dl
   if(co.resu[13]>0.0 and co.resu[13]<0.5) fra[7]+=1.0;//muL
   if(co.resu[14]>0.0 and co.resu[14]<0.5) fra[8]+=1.0;//Mus
   if(co.resu[15]>0.0 and co.resu[15]<0.5) fra[9]+=1.0;//ros_star
   nerr+=1.0;  
   Eru0+=co.resu[0];  ErtE+=co.resu[1];  Erfb+=co.resu[2];    ErpiE+=co.resu[3];  ErtetE+=co.resu[5];  
   Erml+=co.resu[9];  Erdl+=co.resu[10], Ermul+=co.resu[13];  Ermus+=co.resu[14]; Eros  +=co.resu[15];}}  }
   fclose(filg);
   
   
   Eru0=    double(Eru0/ (nerr+eps));  
   Erfb=    double(Erfb/ (nerr+eps));  
   ErtE=    double(ErtE/ (nerr+eps));  
   ErpiE=  double(ErpiE/ (nerr+eps));  
   ErtetE=double(ErtetE/ (nerr+eps));  
   Erml=    double(Erml/ (nerr+eps));  
   Erdl=    double(Erdl/ (nerr+eps));  
   Ermul=   double(Ermul/(nerr+eps));  
   Ermus=   double(Ermus/(nerr+eps));  
   Eros=    double(Eros/ (nerr+eps)); 
   dchi1=double(dchi1/(numd[1]+eps));///<dchiP/Ndata>parallax
   dchi2=double(dchi2/(numd[1]+eps));///<dchiA/Ndata> Astrometric deflection
   Npeak=double(Npeak/(numd[1]+eps));///<Npeak/Ndata>
   mag1 =double(mag1/(numdts+eps));///<mag1/Ndata>
   mag2 =double(mag2/(numdts+eps));///<mag2/Ndata>
   Sepa =double(Sepa/(numdts+eps));///<Separation/Ndata>
   fra1 =double(fra1*100.0/(numd[1]+eps));///fraction of parallax detected
   fra2 =double(fra2*100.0/(numd[1]+eps));///fraction of astrometric detected
   fra3 =double(fra3*100.0/(numd[1]+eps));///fraction of images resolved
   for(int j=0;j<10; ++j){ 
   fra[j]=double(fra[j]*100.0/(numd[1]+eps));///fraction of measuring erroso
   }
   //colo1=double(colo1/(numd[1]+eps));  
   //colo2=double(colo2/(numd[1]+eps));  
  
    for(int i =0; i<2; ++i){
   tE[i]= double(  tE[i]/(numd[i]+eps));//[days]  
   RE[i]= double(  RE[i]/(numd[i]+eps));//[AU] 
  ros[i]= double( ros[i]/(numd[i]+eps));//[]   
 Extv[i]= double(Extv[i]/(numd[i]+eps));//[]   
  piE[i]= double( piE[i]/(numd[i]+eps));//[]     
 tetE[i]= double(tetE[i]/(numd[i]+eps));//[mas]  
   Vt[i]= double(  Vt[i]/(numd[i]+eps));//[km/s]  
   u0[i]= double(  u0[i]/(numd[i]+eps));//[]  
   Ml[i]= double(  Ml[i]/(numd[i]+eps));//[Msun]  
  opdB[i]= double( opdB[i]/(numB[i]+eps))*u0m*u0m;//[] x10^{6}
  opdP[i]= double( opdP[i]/(numP[i]+eps))*u0m*u0m;//[] x10^{6}
   Dl[i]= double(  Dl[i]/(numd[i]+eps));//[kpc]  
   Ds[i]= double(  Ds[i]/(numd[i]+eps));//[kpc]  
   vl[i]= double(  vl[i]/(numd[i]+eps));//[km/s]  
   vs[i]= double(  vs[i]/(numd[i]+eps));//[km/s]  
  mbs[i]= double( mbs[i]/(numd[i]+eps));//[mag]  
   fb[i]= double(  fb[i]/(numd[i]+eps));//[]
 fwhm[i]= double(fwhm[i]/(numd[i]+eps));//[days]
  vsn[i]= double( vsn[i]/(numd[i]+eps));//[]
 DelT[i]= double(DelT[i]/(numd[i]+eps));//[]
murel[i]=double(murel[i]/(numd[i]+eps));//[mas/year]
  Map[i]=double(  Map[i]/(numd[i]+eps));//[mag]
  nbl[i]=double(  nbl[i]/(numd[i]+eps));//[]
  Ext[i]=double(  Ext[i]/(numd[i]+eps));//[mag]
 catg[i]=double( catg[i]/(numd[i]+eps));   
Struc[i]=double(Struc[i]*100.0/(numd[i]+eps));
fraBD[i]=double(fraBD[i]*100.0/(numd[i]+eps));///fraction of parallax detected
fraMS[i]=double(fraMS[i]*100.0/(numd[i]+eps));///fraction of astrometric detected
fraBH[i]=double(fraBH[i]*100.0/(numd[i]+eps));///fraction of images resolve
color[0][i]=double(color[0][i]/(numd[i]+eps));///u-r 
color[1][i]=double(color[1][i]/(numd[i]+eps));///g-i
color[2][i]=double(color[2][i]/(numd[i]+eps));///r-z 
color[3][i]=double(color[3][i]/(numd[i]+eps));/// mag[argm0]-mag[argm1]
  mag[0][i]=double(mag[0][i]/(numd[i]+eps));//u-band
  mag[1][i]=double(mag[1][i]/(numd[i]+eps));//g-band
  mag[2][i]=double(mag[2][i]/(numd[i]+eps));//r-band
  mag[3][i]=double(mag[3][i]/(numd[i]+eps));//i-band
  mag[4][i]=double(mag[4][i]/(numd[i]+eps));//z-band
  mag[5][i]=double(mag[5][i]/(numd[i]+eps));}//y-band  
  
  

///###########################################################################################################################3
///            Indicating the ranges of Lensing parameters 
///###########################################################################################################################
   filg=fopen(testfil,"r");
   for(int j=0; j<icon; ++j){
   fscanf(filg,
   "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
   "%lf %lf %lf %lf %lf %d %lf %lf %lf  %lf %lf %lf %lf %lf %lf  %lf %d  %d  %d %lf %lf %lf %lf %lf %lf\n",//53
   &counter, &FFG, &l.tE, &l.RE, &l.piE, &l.tetE, &l.Vt, &l.u0, &l.Ml,    
   &s.opt, &l.Dl, &s.Ds, &l.vl, &s.vs, &mbase, &fblend, &gg, &l.struc, &s.FWHM, &vsave, &l.DeltaT, &l.murel,
   &co.resu[0],  &co.resu[1], &co.resu[2], &co.resu[3], &co.resu[5], &co.resu[9], &co.resu[10], &co.resu[13],
   &co.resu[14], &s.Map[2], &s.nsbl[2], &co.flagi, &s.Ai[2], &dchiP, &dchiA, &npeak, &Mi1, &Mi2, &sepa, &s.ros, &s.Av, &co.resu[15],
   &l.catg, &fdetB, &fdetP, &s.magb[0],&s.magb[1],&s.magb[2],&s.magb[3],&s.magb[4],&s.magb[5]);//53
   
   
   za[0]=1;
   za[1]=FFG;  

  
   for(int i=0; i<2; ++i){
   DtE[i]   +=pow(l.tE-tE[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);       
   DRE[i]   +=pow(l.RE-RE[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);       
   DpiE[i]  +=pow(l.piE-piE[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);     
   DtetE[i] +=pow(l.tetE-tetE[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);  
   DVt[i]   +=pow(l.Vt-Vt[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);    
   Dros[i]  +=pow(s.ros-ros[i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Du0[i]   +=pow(l.u0-u0[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);       
   DMl[i]   +=pow(l.Ml-Ml[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);       
   Dopd[i]  +=pow(s.opt-opdP[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);     
   DDl[i]   +=pow(l.Dl-Dl[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);     
   DDs[i]   +=pow(s.Ds-Ds[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);    
   DExtv[i] +=pow(s.Av-Extv[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);
   Dvl[i]   +=pow(l.vl-vl[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);       
   Dvs[i]   +=pow(s.vs-vs[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);       
   Dmbs[i]  +=pow(mbase-mbs[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);     
   Dfb[i]   +=pow(fblend-fb[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);    
   Dfwhm[i] +=pow(s.FWHM-fwhm[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);   
   Dvsn[i]  +=pow(vsave-vsn[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);     
   DDelT[i] +=pow(l.DeltaT-DelT[i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Dmurel[i]+=pow(l.murel-murel[i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   DMap[i]  +=pow(s.Map[2]-Map[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);  
   Dnbl[i]  +=pow(s.nsbl[2]-nbl[i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   DExt[i]  +=pow(s.Ai[2]-Ext[i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Dcol[0][i]+=pow(s.magb[0]-s.magb[2]-color[0][i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Dcol[1][i]+=pow(s.magb[1]-s.magb[3]-color[1][i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Dcol[2][i]+=pow(s.magb[2]-s.magb[4]-color[2][i],2.0)*za[i]/fabs(numd[i]+eps-1.0);
   Dcol[3][i]+=pow(s.magb[ls.argm[0]]-s.magb[ls.argm[1]]-color[3][i],2.0)*za[i]/fabs(numd[i]+eps-1.0);
   
   Dmag[0][i]+=pow(s.magb[0]-mag[0][i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Dmag[1][i]+=pow(s.magb[1]-mag[1][i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Dmag[2][i]+=pow(s.magb[2]-mag[2][i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Dmag[3][i]+=pow(s.magb[3]-mag[3][i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Dmag[4][i]+=pow(s.magb[4]-mag[4][i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Dmag[5][i]+=pow(s.magb[5]-mag[5][i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
    
   } }
   fclose(filg);
   
   for(int i=0; i<2; ++i){
   DtE[i]=sqrt(fabs(DtE[i]));      DRE[i]=sqrt(fabs(DRE[i]));      DpiE[i]=sqrt(fabs(DpiE[i]));  
   DtetE[i]=sqrt(fabs(DtetE[i]));  DVt[i]=sqrt(fabs(DVt[i]));      Dros[i]=sqrt(fabs(Dros[i])); 
   Du0[i]=sqrt(fabs(Du0[i]));      DMl[i]=sqrt(fabs(DMl[i]));      Dopd[i]=sqrt(fabs(Dopd[i])); 
   DDl[i]=sqrt(fabs(DDl[i]));      DDs[i]=sqrt(fabs(DDs[i]));      DExtv[i]=sqrt(fabs(DExtv[i])); 
   Dvl[i]=sqrt(fabs(Dvl[i]));      Dvs[i]=sqrt(fabs(Dvs[i]));      Dmbs[i]=sqrt(fabs(Dmbs[i])); 
   Dfb[i]=sqrt(fabs(Dfb[i]));      Dfwhm[i]=sqrt(fabs(Dfwhm[i]));  Dvsn[i]=sqrt(fabs(Dvsn[i])); 
   DDelT[i]=sqrt(fabs(DDelT[i]));  Dmurel[i]=sqrt(fabs(Dmurel[i])); DMap[i]=sqrt(fabs(DMap[i])); 
   Dnbl[i]=sqrt(fabs(Dnbl[i]));    DExt[i]=sqrt(fabs(DExt[i]));
   Dcol[0][i]=sqrt(fabs(Dcol[0][i]));   Dcol[1][i]=sqrt(fabs(Dcol[1][i]));   
   Dcol[2][i]=sqrt(fabs(Dcol[2][i]));   Dcol[3][i]=sqrt(fabs(Dcol[3][i]));
   Dmag[0][i]=sqrt(fabs(Dmag[0][i]));   Dmag[1][i]=sqrt(fabs(Dmag[1][i]));   Dmag[2][i]=sqrt(fabs(Dmag[2][i]));
   Dmag[3][i]=sqrt(fabs(Dmag[3][i]));   Dmag[4][i]=sqrt(fabs(Dmag[4][i]));   Dmag[5][i]=sqrt(fabs(Dmag[5][i]));} 

}
   
   
   
   

///###########################################################################################################################3
///            Saving them in files  MapA1.dat   &   MapB1.dat 
///###########################################################################################################################
   EFF=  double(EFF/(numd[1]+eps));//1/[days]
   Gamma=double(2.0/M_PI * opdP[0] *1.0e-6 *EFF)/u0m * year;//1/[star*year]  
   EffiL=double(numd[1]*100.0/(numd[0]+eps)); //probability of detecting lensing  
   EffiD=double(numd[0]*100.0/(nsim   +eps)); //probability of detecting stars  
   EffiB=double(numB[0]*100.0/(nsim   +eps)); //probability of detecting stars  
   EffiP=double(numP[0]*100.0/(nsim   +eps)); //probability of detecting stars  
   Neven=double(s.nstart * Gamma * 10.0);//deg^{-2}  
   gtot+=Gamma; 
   ntot+=Neven;  
   fil3=fopen(fnGam,"a+");//MapA1.dat
   fprintf(fil3,
   "%d  %d  %.6lf %.6lf  %.6lf %.6lf  %.7lf %.7lf   %.7lf %.7lf  %.6lf %.6lf  %.6lf %.6lf  %.5lf %.5lf  %.8lf %.8lf  "//18
   "%.5lf %.5lf  %.5lf %.5lf  %.6lf %.6lf  %.6lf %.6lf  %.5lf %.5lf  %.7lf %.7lf  %.6lf %.6lf  %.6lf %.6lf  %.6lf %.6lf  "//36
   "%.6lf %.6lf  %.7lf %.7lf  %.5lf %.5lf  %.4lf  %.4lf  %.4lf  %.4lf  %.8lf  %.8lf %.4lf %.4lf  %.2lf %.2lf  %.7lf  %.7lf  "//54
   "%.9lf  %.9lf  %.9lf %.9lf  %.9lf %.9lf  %.9lf %.9lf  %.9lf  %.9lf  " //64
   "%.5lf  %.5lf  %.5lf %.5lf  %.5lf %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  " //74
   "%.4lf  %.4lf  %.5lf  %.5lf  %.5lf  %.4lf  %.5lf  %.5lf %.5lf %d  %d  %d  %d  %d   %.4lf  %.4lf  "//90
   "%.3lf  %.3lf  %.3lf  %.3lf  %.3lf  %.3lf  %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf  "//102
   "%.5lf %.5lf   %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf %.5lf %.5lf %.5lf "//114
   "%.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %d  %d  %.5lf  %.5lf %d\n",//129
   nl1,nl2, tE[0],tE[1], RE[0],RE[1], piE[0],piE[1], tetE[0],tetE[1], Vt[0],Vt[1], u0[0],u0[1], Ml[0],Ml[1], opdB[0],opdB[1],//18
   Dl[0],Dl[1], Ds[0],Ds[1], vl[0],vl[1], vs[0],vs[1], mbs[0],mbs[1], fb[0],fb[1], fwhm[0],fwhm[1],vsn[0],vsn[1], DelT[0],DelT[1],//36
   Struc[0],Struc[1],murel[0],murel[1],Map[0],Map[1],nbl[0],nbl[1],Ext[0],Ext[1],ros[0],ros[1],Extv[0],Extv[1],catg[0],catg[1],
   opdP[0],opdP[1],//54
   EffiB, EffiP, EffiD, EffiL, EFF, Gamma*1.0e7, Neven, log10(s.Rostart[0]+eps), log10(s.Nstart+eps),s.nstart,//64
   Eru0, ErtE, Erfb, ErpiE, ErtetE, Erml, Erdl, Ermul, Ermus, Eros,//74
   dchi1, dchi2,fra1,fra2,fra3,Npeak,mag1,mag2,Sepa,int(ndd),int(nsim),int(numd[0]),int(numd[1]),int(nerr),s.right,s.decli,//90
   fraBD[0],fraBD[1],fraMS[0],fraMS[1],fraBH[0],fraBH[1],color[0][0],color[0][1],color[1][0],color[1][1],color[2][0],color[2][1],//102 
   mag[0][0],mag[0][1],mag[1][0],mag[1][1], mag[2][0],mag[2][1],mag[3][0],mag[3][1],mag[4][0],mag[4][1],mag[5][0],mag[5][1],//114
   fra[0],fra[1],fra[2],fra[3],fra[4],fra[5],fra[6],fra[7],fra[8],fra[9],ls.argm[0],ls.argm[1],color[3][0],color[3][1],ls.flag);//129
   fclose(fil3);  
   
   fil3B=fopen(fnGamB,"a+");//MapB1.dat
   fprintf(fil3B,"%d  %d  %.6lf  %.6lf %.6lf  %.6lf %.7lf  %.7lf  %.7lf  %.7lf  %.6lf %.6lf  %.6lf %.6lf  %.5lf %.5lf  " //16
   "%.8lf   %.8lf  %.5lf %.5lf  %.5lf %.5lf  %.6lf  %.6lf  %.6lf %.6lf  %.5lf %.5lf  %.5lf %.5lf  "//30
   "%.7lf   %.7lf  %.6lf %.6lf  %.6lf %.6lf  %.6lf  %.6lf  %.6lf %.6lf  %.3lf %.3lf "//42
   "%.5lf   %.5lf  %.8lf %.8lf  %.5lf  %.5lf %.5lf  %.5lf  %.5lf %.5lf  %.5lf  %.5lf  "//54
   "%.5lf  %.5lf %.5lf  %.5lf %.5lf  %.5lf %.5lf  %.5lf %.5lf  %.5lf %.5lf  %.5lf  %.5lf  %.5lf\n",//68
    nl1, nl2, DtE[0],DtE[1], DRE[0],DRE[1], DpiE[0],DpiE[1], DtetE[0],DtetE[1], DVt[0],DVt[1], Du0[0],Du0[1], DMl[0],DMl[1],//16
    Dopd[0], Dopd[1], DDl[0],DDl[1], DDs[0],DDs[1], Dvl[0], Dvl[1], Dvs[0], Dvs[1], Dmbs[0], Dmbs[1], Dfb[0], Dfb[1], //30
    Dfwhm[0], Dfwhm[1], Dvsn[0], Dvsn[1], DDelT[0], DDelT[1], Dmurel[0], Dmurel[1], DMap[0], DMap[1], Dnbl[0],Dnbl[1],//42
    DExt[0],DExt[1], Dros[0],Dros[1], DExtv[0],DExtv[1],Dcol[0][0],Dcol[0][1],Dcol[1][0],Dcol[1][1],Dcol[2][0],Dcol[2][1],//54
    Dmag[0][0],Dmag[0][1],Dmag[1][0],Dmag[1][1],Dmag[2][0],Dmag[2][1],Dmag[3][0],Dmag[3][1],Dmag[4][0],Dmag[4][1],
    Dmag[5][0],Dmag[5][1],Dcol[3][0],Dcol[3][1]);//68
    fclose(fil3B);  
   
   cout<<"nsim:  "<<nsim<<"\t Ndetected:  "<<icon<<"\t Nlensing:  "<<nlens<<"\t NError:  "<<nerr<<endl;
   cout<<"numd0:  "<<numd[0]<<"\t numd1:  "<<numd[1]<<endl;
   cout<<"EFF:  "<<EFF<<"\t Gamma:  "<<Gamma<<"\t Neven:  "<<Neven<<endl;
   cout<<"l.tE:  "<<l.tE<<"\t l.Ml:  "<<l.Ml<<"\t catgL:  "<<l.catg<<endl;
   cout<<"========================================================================="<<endl;  
   if(icon>=1 and nlens>=2 and
   (EFF<0.0 or Gamma<0.0 or EffiL<0.0 or EffiD<0.0 or Neven<0.0 or DelT[0]<0.0 or DelT[1]<0.0 or vsn[0]<0.0 or vsn[1]<0.0 or 
   fwhm[0]<0.0 or fwhm[1]<0.0 or fb[0]<0.0 or fb[0]>1.0 or fb[1]<0.0 or fb[1]>1.0 or mbs[0]<0.0 or mbs[1]<0.0  or 
   vs[0]<0.0 or vs[1]<0.0 or vl[0]<0.0 or vl[1]<0.0 or Ds[0]<0.0 or Ds[1]<0.0 or Dl[0]<0.0 or Dl[1]<0.0 or piE[1]<0.0 
   or tetE[1]<0.0 or murel[0]<0.0 or murel[1]<0.0 or ErtE<0.0 or ErpiE<0.0 or ErtetE<0.0 or Erml<0.0 or Erdl<0.0 or 
   Ermul<0.0 or Ermus<0.0 or Eru0<0.0 or icon!=numd[0] or numd[1]!=nlens or Extv[0]<0.0 or Extv[1]<0.0 or DtE[0]<0.0 or DtE[1]<0.0)){
   cout<<"Extv[0]:  "<<Extv[0]<<"\t Extv[1]:  "<<Extv[1]<<endl;
   cout<<"TT:  "<<EFF<<"\t Gamma:  "<<Gamma<<"\t Effic: "<<EffiD<<EffiL<<"\t Neven "<<Neven<<endl;  
   cout<<"ErtE: "<<ErtE<<"\t ErpiE:  "<<ErpiE<<"\t ErtetE: "<<ErtetE<<"\t nlens:  "<<nlens<<endl; 
   cout<<"icon:  "<<icon<<"\t numd[0]:  "<<numd[0]<<"\t numd[1]:  "<<numd[1]<<"\t nsim:  "<<nsim<<endl;  
   cout<<"Erml:  "<<Erml<<"\t Erdl:  "<<Erdl<<"\t Ermul:  "<<Ermul<<"\t Ermus:  "<<Ermus<<endl;
   cout<<"fra3:  "<<fra3<<"\t numd[1]:  "<<numd[1]<<endl;
   cout<<"Npeak:  "<<Npeak<<"\t numdts:   "<<numdts<<endl;
   cout<<"mag1:  "<<mag1<<"\t mag2: "<<mag2<<"\t Sepa:  "<<Sepa<<endl; int rdr; cin>>rdr;}
   }//loop of galactic latitude
   
   
   
   
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///            Calculating the average and STD values for a fixed LONGITUDE values
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   for(int i=0;i<2; ++i){
   fraBD[i]=0.0; fraMS[i]=0.0;  fraBH[i]=0.0; 
   tE[i]=0.0,    RE[i]=0.0,     piE[i]=0.0,  tetE[i]=0.0;  Vt[i]=0.0;  ros[i]=0.0;  Extv[i]=0.0; 
   u0[i]=0.0,    Ml[i]=0.0,     opdB[i]=0.0; opdP[i]=0.0;  Dl[i]=0.0,    Ds[i]=0.0;  catg[i]=0.0;  
   vl[i]=0.0,    vs[i]=0.0,     mbs[i]=0.0,  fb[i]=0.0;    numd[i]=0.0;   numB[i]=0.0;  numP[i]=0.0;  
   fwhm[i]=0.0;  vsn[i]=0.0;    DelT[i]=0.0; Struc[i]=0.0; murel[i]=0.0; 
   Map[i]=0.0;   nbl[i]=0.0;    Ext[i]=0.0;
   for(int j=0;j<4;++j){color[j][i]=0.0;      Dcol[j][i]=0.0;  }
   for(int j=0;j<6;++j){  mag[j][i]=0.0;      Dmag[j][i]=0.0;  }
   DtE[i]=0.0,   DRE[i]=0.0,    DpiE[i]=0.0, DtetE[i]=0.0, DVt[i]=0.0,   Dros[i]=0.0; 
   Du0[i]=0.0,   DMl[i]=0.0,    Dopd[i]=0.0, DDl[i]=0.0,   DDs[i]=0.0,   DExtv[i]=0.0; 
   Dvl[i]=0.0,   Dvs[i]=0.0,    Dmbs[i]=0.0, Dfb[i]=0.0,   Dfwhm[i]=0.0, Dvsn[i]=0.0; 
   DDelT[i]=0.0, Dmurel[i]=0.0, DMap[i]=0.0, Dnbl[i]=0.0,  DExt[i]=0.0;}      
   for(int i=0;i<10; ++i)    fra[i]=0.0;   
   ErtE=0.0;     ErpiE=0.0;  ErtetE=0.0;  Erml=0.0; Erfb=0.0;   nerr=0.0; 
   Erdl=0.0;     Ermul=0.0;  Ermus=0.0;   Eru0=0.0; Eros=0.0;   EffiL=0.0;
   dchi1=0.0;    dchi2=0.0;  fra1=0.0;    fra2=0.0;  fra3=0.0;   
   Npeak=0.0;    mag1=0.0;   mag2=0.0;    Sepa=0.0;  numdts=0.0; 
   filg2=fopen(testfil2,"r");
   for(int i=0; i<nlongi; ++i){
   fscanf(filg2,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
   "%lf %lf %lf %lf %lf %d %lf %lf %lf  %lf %lf %lf %lf %lf %lf  %lf %d  %d  %d %lf %lf %lf %lf %lf %lf\n",//53
   &counter, &FFG, &l.tE, &l.RE, &l.piE, &l.tetE, &l.Vt, &l.u0, &l.Ml,    
   &s.opt, &l.Dl, &s.Ds, &l.vl, &s.vs, &mbase, &fblend, &gg, &l.struc, &s.FWHM, &vsave, &l.DeltaT, &l.murel,
   &co.resu[0],  &co.resu[1], &co.resu[2], &co.resu[3], &co.resu[5], &co.resu[9], &co.resu[10], &co.resu[13],
   &co.resu[14], &s.Map[2], &s.nsbl[2], &co.flagi, &s.Ai[2], &dchiP, &dchiA, &npeak, &Mi1, &Mi2, &sepa, &s.ros, &s.Av, &co.resu[15],
   &l.catg, &fdetB, &fdetP, &s.magb[0],&s.magb[1],&s.magb[2],&s.magb[3],&s.magb[4],&s.magb[5]);//53
   
   if(FFG>0){
   if(l.struc==4 or l.struc==3)          Struc[1]+=1.0;//Halo-lensing fraction 
   tE[1]+=l.tE;       RE[1]+=l.RE;       piE[1]+=l.piE;     tetE[1]+=l.tetE;   Vt[1]+=l.Vt;  ros[1]+=s.ros; 
   u0[1]+=l.u0;       Ml[1]+=l.Ml;       Dl[1]+= l.Dl;      Ds[1]+=s.Ds;  Extv[1]+=s.Av;
   vl[1]+=l.vl;       vs[1]+=s.vs;       mbs[1]+=mbase;     fb[1]+=fblend;     numd[1]+=1.0;
   fwhm[1]+=s.FWHM;   vsn[1]+=vsave;     DelT[1]+=l.DeltaT; murel[1]+=l.murel;
   Map[1]+=s.Map[2];  nbl[1]+=s.nsbl[2]; Ext[1]+=s.Ai[2];   catg[1]+=l.catg*1.0;
   if(fdetB>=2)       {opdB[1]+=s.opt;   numB[1]+=1.0; }
   if(fdetP>=2)       {opdP[1]+=s.opt;   numP[1]+=1.0; }
   if(l.catg==1)      fraBD[1]+=1.0; 
   if(l.catg==2)      fraMS[1]+=1.0;  
   if(l.catg==3)      fraBH[1]+=1.0; 
   
   color[0][1] +=float(s.magb[0]-s.magb[2]);///u-r 
   color[1][1] +=float(s.magb[1]-s.magb[3]);///g-i
   color[2][1] +=float(s.magb[2]-s.magb[4]);///r-z
      
   mag[0][1]   +=float(s.magb[0]);//u-band
   mag[1][1]   +=float(s.magb[1]);//g-band
   mag[2][1]   +=float(s.magb[2]);//r-band
   mag[3][1]   +=float(s.magb[3]);//i-band
   mag[4][1]   +=float(s.magb[4]);//z-band
   mag[5][1]   +=float(s.magb[5]);//y-band
   
   dchi1+=dchiP;      
   dchi2+=dchiA; 
   Npeak+=npeak;
   if(Mi1>0.0 and Mi2>0.0 and sepa>0.0 and npeak>0.0){ 
   numdts+=1.0;  
   mag1 +=Mi1;   
   mag2 +=Mi2;    
   Sepa +=sepa;}
   if(dchiP>=2.0) fra1+=1.0;///parallax detected >2Ndata
   if(dchiA>=2.0) fra2+=1.0;///astrometric trajectory detected >3Ndata
   if(npeak>=2.0) fra3+=1.0;///two images are resolved in at least 3 data points 

   if(co.flagi>0){
   if(co.resu[0]>0.0 and co.resu[0]<0.5)   fra[0]+=1.0;//u0
   if(co.resu[1]>0.0 and co.resu[1]<0.5)   fra[1]+=1.0;//tE
   if(co.resu[2]>0.0 and co.resu[2]<0.5)   fra[2]+=1.0;//fb
   if(co.resu[3]>0.0 and co.resu[3]<0.5)   fra[3]+=1.0;//piE
   if(co.resu[5]>0.0 and co.resu[5]<0.5)   fra[4]+=1.0;//tetE
   if(co.resu[9]>0.0 and co.resu[9]<0.5)   fra[5]+=1.0;//ML
   if(co.resu[10]>0.0 and co.resu[10]<0.5) fra[6]+=1.0;//Dl
   if(co.resu[13]>0.0 and co.resu[13]<0.5) fra[7]+=1.0;//muL
   if(co.resu[14]>0.0 and co.resu[14]<0.5) fra[8]+=1.0;//Mus
   if(co.resu[15]>0.0 and co.resu[15]<0.5) fra[9]+=1.0;//ros_star
   nerr +=1.0;  
   Eru0+=co.resu[0];  ErtE+=co.resu[1];  Erfb+=co.resu[2];    ErpiE+=co.resu[3];  ErtetE+=co.resu[5];  
   Erml+=co.resu[9];  Erdl+=co.resu[10], Ermul+=co.resu[13];  Ermus+=co.resu[14]; Eros  +=co.resu[15];}}}
   fclose(filg2); 
   
   Eru0=    double(Eru0/ (nerr+eps));  
   Erfb=    double(Erfb/ (nerr+eps));  
   ErtE=    double(ErtE/ (nerr+eps));  
   ErpiE=  double(ErpiE/ (nerr+eps));  
   ErtetE=double(ErtetE/ (nerr+eps));  
   Erml=    double(Erml/ (nerr+eps));  
   Erdl=    double(Erdl/ (nerr+eps));  
   Ermul=   double(Ermul/(nerr+eps));  
   Ermus=   double(Ermus/(nerr+eps));  
   Eros=    double(Eros/ (nerr+eps)); 
   
   dchi1=double(dchi1/(numd[1]+eps));///<dchiP/Ndata>parallax
   dchi2=double(dchi2/(numd[1]+eps));///<dchiA/Ndata> Astrometric deflection
   Npeak=double(Npeak/(numd[1]+eps));///<Npeak/Ndata>
   
   mag1 =double(mag1/(numdts+eps));///<mag1/Ndata>
   mag2 =double(mag2/(numdts+eps));///<mag2/Ndata>
   Sepa =double(Sepa/(numdts+eps));///<Separation/Ndata>
   fra1 =double(fra1*100.0/(numd[1]+eps));///fraction of parallax detected
   fra2 =double(fra2*100.0/(numd[1]+eps));///fraction of astrometric detected
   fra3 =double(fra3*100.0/(numd[1]+eps));///fraction of images resolved
   for(int j=0;j<10; ++j){ fra[j]=double(fra[j]*100.0/(numd[1]+eps));}
   
   for(int i=1; i<2; ++i){
   tE[i]= double(  tE[i]/(numd[i]+eps));//[days]  
   RE[i]= double(  RE[i]/(numd[i]+eps));//[AU] 
  ros[i]= double( ros[i]/(numd[i]+eps));//[]   
 Extv[i]= double(Extv[i]/(numd[i]+eps));//[]   
  piE[i]= double( piE[i]/(numd[i]+eps));//[]     
 tetE[i]= double(tetE[i]/(numd[i]+eps));//[mas]  
   Vt[i]= double(  Vt[i]/(numd[i]+eps));//[km/s]  
   u0[i]= double(  u0[i]/(numd[i]+eps));//[]  
   Ml[i]= double(  Ml[i]/(numd[i]+eps));//[Msun]  
  opdB[i]= double( opdB[i]/(numB[i]+eps))*u0m*u0m;//[] x10^{6}
  opdP[i]= double( opdP[i]/(numP[i]+eps))*u0m*u0m;//[] x10^{6}
   Dl[i]= double(  Dl[i]/(numd[i]+eps));//[kpc]  
   Ds[i]= double(  Ds[i]/(numd[i]+eps));//[kpc]  
   vl[i]= double(  vl[i]/(numd[i]+eps));//[km/s]  
   vs[i]= double(  vs[i]/(numd[i]+eps));//[km/s]  
  mbs[i]= double( mbs[i]/(numd[i]+eps));//[mag]  
   fb[i]= double(  fb[i]/(numd[i]+eps));//[]
 fwhm[i]= double(fwhm[i]/(numd[i]+eps));//[days]
  vsn[i]= double( vsn[i]/(numd[i]+eps));//[]
 DelT[i]= double(DelT[i]/(numd[i]+eps));//[]
murel[i]=double(murel[i]/(numd[i]+eps));//[mas/year]
  Map[i]=double(  Map[i]/(numd[i]+eps));//[mag]
  nbl[i]=double(  nbl[i]/(numd[i]+eps));//[]
  Ext[i]=double(  Ext[i]/(numd[i]+eps));//[mag]
 catg[i]=double( catg[i]/(numd[i]+eps));   
Struc[i]=double(Struc[i]*100.0/(numd[i]+eps));
fraBD[i]=double(fraBD[i]*100.0/(numd[i]+eps));///fraction of parallax detected
fraMS[i]=double(fraMS[i]*100.0/(numd[i]+eps));///fraction of astrometric detected
fraBH[i]=double(fraBH[i]*100.0/(numd[i]+eps));///fraction of images resolve
color[0][i]=double(color[0][i]/(numd[i]+eps));///u-r 
color[1][i]=double(color[1][i]/(numd[i]+eps));///g-i
color[2][i]=double(color[2][i]/(numd[i]+eps));///r-z 

  mag[0][i]=double(mag[0][i]/(numd[i]+eps));//u-band
  mag[1][i]=double(mag[1][i]/(numd[i]+eps));//g-band
  mag[2][i]=double(mag[2][i]/(numd[i]+eps));//r-band
  mag[3][i]=double(mag[3][i]/(numd[i]+eps));//i-band
  mag[4][i]=double(mag[4][i]/(numd[i]+eps));//z-band
  mag[5][i]=double(mag[5][i]/(numd[i]+eps));}//y-band
  
      EffiL=double(nlongi*100.0/(nlongid+eps)); //probability of detecting lensing  
      EffiD=double(nlongid*100.0/(nlongis+eps)); //probability of detecting lensing 
  ///################################################################################################
   filg2=fopen(testfil2,"r");
   for(int j=0; j<nlongi; ++j){
   fscanf(filg2,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
   "%lf %lf %lf %lf %lf %d %lf %lf %lf  %lf %lf %lf %lf %lf %lf  %lf %d  %d  %d %lf %lf %lf %lf %lf %lf\n",//53
   &counter, &FFG, &l.tE, &l.RE, &l.piE, &l.tetE, &l.Vt, &l.u0, &l.Ml,    
   &s.opt, &l.Dl, &s.Ds, &l.vl, &s.vs, &mbase, &fblend, &gg, &l.struc, &s.FWHM, &vsave, &l.DeltaT, &l.murel,
   &co.resu[0],  &co.resu[1], &co.resu[2], &co.resu[3], &co.resu[5], &co.resu[9], &co.resu[10], &co.resu[13],
   &co.resu[14], &s.Map[2], &s.nsbl[2], &co.flagi, &s.Ai[2], &dchiP, &dchiA, &npeak, &Mi1, &Mi2, &sepa, &s.ros, &s.Av, &co.resu[15],
   &l.catg, &fdetB, &fdetP, &s.magb[0],&s.magb[1],&s.magb[2],&s.magb[3],&s.magb[4],&s.magb[5]);//53
   
   za[0]=1;
   za[1]=FFG;
   for(int i=1; i<2; ++i){
   DtE[i]   +=pow(l.tE-tE[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);       
   DRE[i]   +=pow(l.RE-RE[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);       
   DpiE[i]  +=pow(l.piE-piE[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);     
   DtetE[i] +=pow(l.tetE-tetE[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);  
   DVt[i]   +=pow(l.Vt-Vt[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);    
   Dros[i]  +=pow(s.ros-ros[i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Du0[i]   +=pow(l.u0-u0[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);       
   DMl[i]   +=pow(l.Ml-Ml[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);       
   Dopd[i]  +=pow(s.opt-opdP[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);     
   DDl[i]   +=pow(l.Dl-Dl[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);     
   DDs[i]   +=pow(s.Ds-Ds[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);    
   DExtv[i] +=pow(s.Av-Extv[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);
   Dvl[i]   +=pow(l.vl-vl[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);       
   Dvs[i]   +=pow(s.vs-vs[i]  ,2.0)*za[i]/fabs(numd[i]+eps-1.0);       
   Dmbs[i]  +=pow(mbase-mbs[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);     
   Dfb[i]   +=pow(fblend-fb[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);    
   Dfwhm[i] +=pow(s.FWHM-fwhm[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);   
   Dvsn[i]  +=pow(vsave-vsn[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);     
   DDelT[i] +=pow(l.DeltaT-DelT[i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Dmurel[i]+=pow(l.murel-murel[i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   DMap[i]  +=pow(s.Map[2]-Map[i],2.0)*za[i]/fabs(numd[i]+eps-1.0);  
   Dnbl[i]  +=pow(s.nsbl[2]-nbl[i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   DExt[i]  +=pow(s.Ai[2]-Ext[i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   
   Dcol[0][i]+=pow(s.magb[0]-s.magb[2]-color[0][i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Dcol[1][i]+=pow(s.magb[1]-s.magb[3]-color[1][i],2.0)*za[i]/fabs(numd[i]+eps-1.0); 
   Dcol[2][i]+=pow(s.magb[2]-s.magb[4]-color[2][i],2.0)*za[i]/fabs(numd[i]+eps-1.0);
   
   Dmag[0][i]+=pow(s.magb[0]-mag[0][i],2.0)*za[i]/fabs(numd[i]+eps-1.0);
   Dmag[1][i]+=pow(s.magb[1]-mag[1][i],2.0)*za[i]/fabs(numd[i]+eps-1.0);
   Dmag[2][i]+=pow(s.magb[2]-mag[2][i],2.0)*za[i]/fabs(numd[i]+eps-1.0);
   Dmag[3][i]+=pow(s.magb[3]-mag[3][i],2.0)*za[i]/fabs(numd[i]+eps-1.0);
   Dmag[4][i]+=pow(s.magb[4]-mag[4][i],2.0)*za[i]/fabs(numd[i]+eps-1.0);
   Dmag[5][i]+=pow(s.magb[5]-mag[5][i],2.0)*za[i]/fabs(numd[i]+eps-1.0);
                  
   }}
   fclose(filg2);
   for(int i=1; i<2; ++i){
   DtE[i]=sqrt(fabs(DtE[i]));      DRE[i]=sqrt(fabs(DRE[i]));      DpiE[i]=sqrt(fabs(DpiE[i]));  
   DtetE[i]=sqrt(fabs(DtetE[i]));  DVt[i]=sqrt(fabs(DVt[i]));      Dros[i]=sqrt(fabs(Dros[i])); 
   Du0[i]=sqrt(fabs(Du0[i]));      DMl[i]=sqrt(fabs(DMl[i]));      Dopd[i]=sqrt(fabs(Dopd[i])); 
   DDl[i]=sqrt(fabs(DDl[i]));      DDs[i]=sqrt(fabs(DDs[i]));      DExtv[i]=sqrt(fabs(DExtv[i])); 
   Dvl[i]=sqrt(fabs(Dvl[i]));      Dvs[i]=sqrt(fabs(Dvs[i]));      Dmbs[i]=sqrt(fabs(Dmbs[i])); 
   Dfb[i]=sqrt(fabs(Dfb[i]));      Dfwhm[i]=sqrt(fabs(Dfwhm[i]));  Dvsn[i]=sqrt(fabs(Dvsn[i])); 
   DDelT[i]=sqrt(fabs(DDelT[i]));  Dmurel[i]=sqrt(fabs(Dmurel[i])); DMap[i]=sqrt(fabs(DMap[i])); 
   Dnbl[i]=sqrt(fabs(Dnbl[i]));    DExt[i]=sqrt(fabs(DExt[i]));
   Dcol[0][i]=sqrt(fabs(Dcol[0][i])); Dcol[1][i]=sqrt(fabs(Dcol[1][i]));  Dcol[2][i]=sqrt(fabs(Dcol[2][i]));
   
   Dmag[0][i]=sqrt(fabs(Dmag[0][i])); Dmag[1][i]=sqrt(fabs(Dmag[1][i]));  Dmag[2][i]=sqrt(fabs(Dmag[2][i]));
   Dmag[3][i]=sqrt(fabs(Dmag[3][i])); Dmag[4][i]=sqrt(fabs(Dmag[4][i]));  Dmag[5][i]=sqrt(fabs(Dmag[5][i])); } 
   fil3C=fopen(fnGamC,"a+");//Map1C.dat
    
    
   fprintf(fil3C,
   "%d %.5lf %.5lf  %.5lf %.5lf  %.7lf %.7lf   %.7lf %.7lf  %.5lf %.5lf  %.5lf %.5lf  %.5lf %.5lf  %.5lf %.5lf  "//17
   "%.5lf %.5lf  %.5lf %.5lf  %.5lf %.5lf  %.5lf %.5lf  %.5lf %.5lf  %.5lf %.5lf  %.6lf %.6lf  %.6lf  %.6lf "//33
   "%.6lf %.6lf  %.5lf %.5lf  %.5lf %.5lf  %.4lf  %.4lf  %.7lf  %.7lf  %.5lf  %.5lf %.8lf %.8lf  %.8lf  "//48
   "%.5lf  %.5lf %.6lf  %.6lf  %.6lf %.6lf  %.6lf %.6lf  %.6lf %.6lf  %.6lf  %.6lf %.6lf  %.7lf  %.7lf  " //63
   "%.4lf  %.4lf  %.5lf  %.5lf  %.5lf  %.4lf  %.5lf  %.5lf %.5lf %d   %d  "//74
   "%.3lf %.3lf %.3lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf " //85
   "%.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf  " //95
   "%.5lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %d  %d  %d\n",//108
   nl1, tE[1],DtE[1], RE[1],DRE[1], piE[1],DpiE[1], tetE[1],DtetE[1], Vt[1],DVt[1], u0[1],Du0[1], Ml[1],DMl[1],Dl[1],DDl[1],//17 
   Ds[1],DDs[1], vl[1],Dvl[1], vs[1],Dvs[1], mbs[1],Dmbs[1], fb[1],Dfb[1], fwhm[1],Dfwhm[1],vsn[1],Dvsn[1], DelT[1],DDelT[1],//33
   murel[1],Dmurel[1],Map[1],DMap[1],nbl[1],Dnbl[1],Ext[1],DExt[1],ros[1],Dros[1],Extv[1],DExtv[1],opdB[1],opdP[1],Dopd[1],//48
   catg[1], Struc[1],Eru0, ErtE, Erfb, ErpiE, ErtetE, Erml, Erdl, Ermul, Ermus, Eros,EffiL, ntot, gtot, //63
   dchi1, dchi2,fra1,fra2,fra3,Npeak,mag1,mag2,Sepa, int(numd[1]), int(nerr),//74
   fraBD[1],fraMS[1],fraBH[1],color[0][1],Dcol[0][1],color[1][1],Dcol[1][1],color[2][1],Dcol[2][1],mag[0][1],Dmag[0][1],//85
   mag[1][1],Dmag[1][1],mag[2][1],Dmag[2][1],mag[3][1],Dmag[3][1],mag[4][1],Dmag[4][1],mag[5][1],Dmag[5][1],//95
   EffiD,fra[0],fra[1],fra[2],fra[3],fra[4],fra[5],fra[6],fra[7],fra[8],fra[9],ls.argm[0],ls.argm[1], ls.flag);//109
   fclose(fil3C); 
   

   }//end of loop galactic longitude
   return(0);
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//               tE & ML & DL & U0 Functions                          //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
int FunctE(lens & l){
   int gg=-1;  
   if(l.tE<=l.tEs[0]){  
   gg=0; 
   }
   else if(l.tE>=l.tEs[GG])  gg=int(GG); 
   else{
   for(int i=1;  i<=GG; ++i){ 
   if(double((l.tE-l.tEs[i-1])*(l.tE-l.tEs[i]))<0.0 or l.tE==l.tEs[i-1]){gg=int(i-1);  break;}}}
   if(gg<0 or gg>GG or l.tE<=0.0 or l.tEs[1]<=0.0 or l.tEs[10]<=0.0){
   cout<<"Error gg:  "<<gg<<"\t tE:   "<<l.tE<<endl; 
   cout<<"tE_min:  "<<tE_min<<"\t tE_max:   "<<tE_max<<endl; 
   for(int i=0;i<GG;++i) cout<<"tEs[i]:  "<<l.tEs[i]<<endl;
   int iis; cin>>iis;} 

   return(gg); 
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH 
int FuncMl(lens & l){

   int gg=-1;  
   if(l.Ml<=Ml_min)      gg=0; 
   else if(l.Ml>=Ml_max) gg=int(GG); 
   else{
   for(int i=1;  i<=GG; ++i){ 
   if(double((l.Ml-l.Mls[i-1])*(l.Ml-l.Mls[i]))<0.0 or l.Ml==l.Mls[i-1]){gg=int(i-1);  break;}}}
   
   if(gg<0 or gg>GG or l.Ml<0.0 or l.Ml==0.0 or l.Mls[10]<=0.0){
   cout<<"Error gg:  "<<gg<<"\t Ml:   "<<l.Ml<<endl; 
   cout<<"Ml_min:  "<<Ml_min<<"\t Ml_max:   "<<Ml_max<<endl; 
   for(int i=0;i<GG;++i)   cout<<"Mls[i]:  "<<l.Mls[i]<<endl;
   int iis; cin>>iis;} 

   return(gg); 
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH 
int FuncPi(lens & l){

   int gg=-1;  
   double lpi=log10(fabs(l.pirel));  
   if(lpi<=pi_min)      gg=0; 
   else if(lpi>=pi_max) gg=int(GG); 
   else{
   for(int i=1;  i<=GG; ++i){ 
   if(double((lpi-l.pis[i-1])*(lpi-l.pis[i]))<0.0 or lpi==l.pis[i-1]){gg=int(i-1);  break;}}}
   
   if(gg<0 or gg>GG or l.pirel==0.0 or l.pirel<0.0){
   cout<<"Error gg:  "<<gg<<"\t log(pirel):   "<<lpi<<endl; 
   cout<<"pi_min:  "<<pi_min<<"\t pi_max:   "<<pi_max<<endl; 
   for(int i=0;i<GG;++i) cout<<"pis[i]:  "<<l.pis[i]<<endl;
   int iis; cin>>iis;} 

   return(gg); 
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH 
int Funcu0(lens & l){

   int gg=-1;  
   if(l.u0<=u0_min)      gg=0; 
   else if(l.u0>=u0_max) gg=int(GG); 
   else{
   for(int i=1;  i<=GG; ++i){ 
   if(double((l.u0-l.u0s[i-1])*(l.u0-l.u0s[i]))<0.0 or l.u0==l.u0s[i-1]){gg=int(i-1);  break;}}}
   
   if(gg<0 or gg>GG or l.u0<0.0 or l.u0>u0m or l.u0s[10]<=0.0){
   cout<<"Error gg:  "<<gg<<"\t u0:   "<<l.u0<<endl; 
   cout<<"u0_min:  "<<u0_min<<"\t u0_max:   "<<u0_max<<endl; 
   for(int i=0; i<GG; ++i) cout<<"u0s[i]:  "<<l.u0s[i]<<endl;
   int iis; cin>>iis;} 

   return(gg); 
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
int FuncMu(lens & l){

   double mur=double(l.murel*year);//mas/years 
   int gg=-1;  
   if(mur<=mu_min)      gg=0; 
   else if(mur>=mu_max) gg=int(GG); 
   else{
   for(int i=1;  i<=GG; ++i){ 
   if(double((mur-l.mus[i-1])*(mur-l.mus[i]))<0.0 or mur==l.mus[i-1]){gg=int(i-1);  break;}}}
   
   if(gg<0 or gg>GG or mur<0.0 or l.murel==0.0 or l.mus[10]<=0.0){
   cout<<"Error gg:  "<<gg<<"\t        mu_rel:   "<<mur<<endl; 
   cout<<"mu_min:  "<<mu_min<<"\t   mu_max:   "<<mu_max<<endl; 
   for(int i=0; i<GG; ++i) cout<<"mus[i]:  "<<l.mus[i]<<endl;
   int iis; cin>>iis;} 

   return(gg); 
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  
int FuncMb(lens & l, double gadr){

   int gg=-1;  
   if(gadr<=mb_min)       gg=0; 
   else if(gadr>=mb_max)  gg=int(GG); 
   else{
   for(int i=1;  i<=GG; ++i){ 
   if(double((gadr-l.mbs[i-1])*(gadr-l.mbs[i]))<0.0 or gadr==l.mbs[i-1]){gg=int(i-1);  break;}}}
   
   if(gg<0 or gg>GG or gadr<=0.0 or l.mbs[10]<=0.0){
   cout<<"Error gg:  "<<gg<<"\t mbase[2]:   "<<gadr<<endl; 
   cout<<"mb_min:  "<<mb_min<<"\t mb_max:   "<<mb_max<<endl; 
   for(int i=0; i<GG; ++i) cout<<"u0s[i]:  "<<l.mbs[i]<<endl;
   int iis; cin>>iis;} 

   return(gg); 
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH 
int FuncFb(lens & l, double blen){

   int gg=-1;  
   if(blen<=fb_min)       gg=0; 
   else if(blen>=fb_max)  gg=int(GG); 
   else{
   for(int i=1;  i<=GG; ++i){ 
   if(double((blen-l.fbs[i-1])*(blen-l.fbs[i]))<0.0 or blen==l.fbs[i-1]){gg=int(i-1);  break;}}}
   if(gg<0 or gg>GG or blen<0.0 or blen>1.0 or l.fbs[10]<0.0){
   cout<<"Error gg:  "<<gg<<"\t fblend[2]:   "<<blen<<endl; 
   cout<<"fb_min:  "<<fb_min<<"\t fb_max:   "<<fb_max<<endl; 
   for(int i=0; i<GG; ++i) cout<<"fbs[i]:  "<<l.fbs[i]<<endl;
   int iis; cin>>iis;} 

   return(gg); 
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH 
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Fisher and Covariance matrixes             //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void FisherM(source & s, lens & l, astromet & as,  covarian & co, int ndw,  VBBinaryLensing & vbb)
{

   //VBBinaryLensing vbb;
   vbb.Tol=1.e-4;
   vbb.a1 =0.0;  
   vbb.LoadESPLTable("./ESPL.tbl");


   co.flagi=+1; 
   int tt, tev; 
   if(s.fb[0]<0.15)       {co.bb[0]=+0.07; co.bb[1]=+0.15;} 
   else if(s.fb[0]<0.85)  {co.bb[0]=-0.07; co.bb[1]=+0.07;}
   else                   {co.bb[0]=-0.07; co.bb[1]=-0.15;}
   double sum=0.0 ;  

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///        Photometry                                  ///
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH   
   co.Delta1[0]= 0.1507586576;//u0
   co.Delta1[1]= l.tE*0.254674;//tE [days]
   co.Delta1[3]= 0.2509463534656*l.piE;//piE 
   co.Delta1[4]= 3.0*M_PI/180.0; //for xi [radian]
   co.Delta1[5]=fabs(s.ros*0.708745475); //for rho* 
   
   for(int j=0; j<Nx; ++j){
   for(int k=0; k<Nx; ++k){
   co.inputA[j][k]=0.0; co.inverA[j][k]=0.0; 
   co.adjunA[j][k]=0.0; co.tempoA[j][k]=0.0; co.summ[j][k]=0.0;}} 

      
      
   for(int i=0; i<ndw; ++i){//data
   tt=int(l.tele[i]);//telescope[0,1] LSST, ELT
   
   for(int j=0; j<Nx; ++j){   
   for(int h=0; h<2;  ++h){
   
   if(j==0){co.diff=double(+co.Delta1[j]*sig[h] );      l.u0+=co.diff;}
   if(j==1){co.diff=double(+co.Delta1[j]*sig2[h]);      l.tE+=co.diff;}//+
   if(j==2){co.diff=double(co.bb[h]);               s.fb[tt]+=co.diff;}  
   if(j==3){co.diff=double(+co.Delta1[j]*sig2[h]);     l.piE+=co.diff;}//+
   if(j==4){co.diff=double(+co.Delta1[j]*sig[h] );      s.xi+=co.diff;}  
   if(j==5){co.diff=double(+co.Delta1[j]*sig2[h]);     s.ros+=co.diff;}//+
  
  
   lightcurve(s,l,as, l.timn[i]);
   
   if(s.ut<double(30.0*s.ros))  s.Astar=vbb.ESPLMag2(s.ut, s.ros);
   else                         s.Astar=double(s.ut*s.ut+2.0)/sqrt(s.ut*s.ut*(s.ut*s.ut+4.0));
   
   
   
   co.magw= s.mbs[tt]-2.5*log10(s.Astar*s.fb[tt] + 1.0 - s.fb[tt]); 
   co.derm1[h]=double(co.magw-l.magn[i])/co.diff;
   
   if(l.tE<=0.0 or s.fb[tt]<=0.0 or s.fb[tt]>1.0 or l.piE<=0.0 or l.u0==0.0 or co.diff==0.0  or l.errm[i]<=0.0 or s.ros<0.0){
   cout<<"Big ERROR(1) l.tE:"<<l.tE<<"\t s.fb: "<<s.fb[tt]<<endl;  
   cout<<"l.piE: "<<l.piE<<"\t u0:"<<l.u0<<"\t diff:  "<<co.diff<<"\t ros:  "<<s.ros<<endl;
   cout<<"magw: "<<co.magw<<"\t l.magn[i]"<<l.magn[i]<<"\t errM:  "<<l.errm[i]<<endl; cin>>tev;}
  
  
   if(j==0)      l.u0-=co.diff;  
   if(j==1)      l.tE-=co.diff;  
   if(j==2)  s.fb[tt]-=co.diff;   
   if(j==3)     l.piE-=co.diff;
   if(j==4)      s.xi-=co.diff;
   if(j==5)     s.ros-=co.diff;}
   
   co.derm1f=double(co.derm1[0]+co.derm1[1])*0.5; 
   
   for(int k=0; k<=j; ++k){
   for(int h=0; h<2;  ++h){
   
   if(k==0){co.diff=double(+co.Delta1[k]*sig[h] );     l.u0+=co.diff;}
   if(k==1){co.diff=double(+co.Delta1[k]*sig2[h]);     l.tE+=co.diff;} 
   if(k==2){co.diff=double(co.bb[h]);              s.fb[tt]+=co.diff;}  
   if(k==3){co.diff=double(+co.Delta1[k]*sig2[h]);    l.piE+=co.diff;}
   if(k==4){co.diff=double(+co.Delta1[k]*sig[h] );     s.xi+=co.diff;}   
   if(k==5){co.diff=double(+co.Delta1[k]*sig2[h]);    s.ros+=co.diff;}//+
   
   
   lightcurve(s,l,as, l.timn[i]);
   if(s.ut<double(30.0*s.ros))    s.Astar=vbb.ESPLMag2(s.ut, s.ros);
   else                           s.Astar=double(s.ut*s.ut+2.0)/sqrt(s.ut*s.ut*(s.ut*s.ut+4.0));
   
   co.magw=s.mbs[tt]-2.5*log10(s.Astar*s.fb[tt]+1.0-s.fb[tt]); 
   co.derm2[h]=double(co.magw-l.magn[i])/co.diff;
   
   if(l.tE<=0.0 or s.fb[tt]<0.0 or s.fb[tt]>1.0 or l.piE<=0.0 or l.u0==0.0 or co.diff==0.0  or l.errm[i]<=0.0 or s.ros<0.0){
   cout<<"Big ERROR(2) l.tE:"<<l.tE<<"\t s.fb: "<<s.fb[tt]<<endl;  
   cout<<"l.piE: "<<l.piE<<"\t u0:"<<l.u0<<"\t diff:  "<<co.diff<<endl;
   cout<<"magw: "<<co.magw<<"\t l.magn[i]"<<l.magn[i]<<"\t errM:  "<<l.errm[i]<<endl; cin>>tev;}
  
   if(k==0)      l.u0-=co.diff;  
   if(k==1)      l.tE-=co.diff;  
   if(k==2)  s.fb[tt]-=co.diff;   
   if(k==3)     l.piE-=co.diff;
   if(k==4)      s.xi-=co.diff;
   if(k==5)     s.ros-=co.diff;}
   
   co.derm2f=double(co.derm2[0]+co.derm2[1])*0.5;
   co.inputA[j][k] += co.derm1f*co.derm2f/(l.errm[i]*l.errm[i]);}
   }//end of for J
   }//end of data for
   
   for(int j=0; j<Nx;   ++j){
   for(int k=(j+1);k<Nx;++k) co.inputA[j][k]=co.inputA[k][j];}

 
   //cout<<"Covariance matrix Photometry"<<endl; 
   for(int i=0; i<Nx; ++i){ 
   for(int j=0; j<Nx; ++j){
   //cout<<co.inputA[i][j]<<"\t   ";  
   if(co.inputA[i][j]==0.0) co.inputA[i][j]=RandR(0.0,1.0)*pow(10.0,-10.0);  }}
   //cout<<"\n"; 
   //cout<<"*************Input*******************"<<endl;
   
   
   
   inverse(co,0);   
   //for(int i=0; i<Nx; ++i){ for(int j=0; j<Nx; ++j){cout<<co.inverA[i][j]<<"\t   ";}  cout<<"\n"; } 
   //cout<<"************Inverse*****************"<<endl;  
   
   for(int i=0; i<Nx; ++i){
   for(int j=0; j<Nx; ++j){
   co.summ[i][j]=0.0;  
   for(int k=0; k<Nx;  ++k){co.summ[i][j]+= co.inputA[i][k] *co.inverA[k][j]; }}}
   ///if(co.summ[i][j]<1.0e-4) co.summ[i][j]=0.0;}}
   
   
   
   for(int i=0; i<Nx; ++i){
   for(int j=0; j<Nx; ++j){
   if((i==j and fabs(co.summ[i][j]-1.0)>0.2) or (i!=j and fabs(co.summ[i][j]-0.0)>0.2)){
   //cout<<"Error_Inverse sum:  "<<co.summ[i][j]<<"\t i:  "<<i<<"\t j:  "<<j<<endl; 
   //cout<<""<<l.t0<<"\t"<<l.u0<<"\t"<<l.tE<<"\t"<<s.xi*M_PI/180.0<<"\t"<<s.fb[tt]<<"\t"<<s.mbs[tt]<<"\t"<<l.piE<<endl; 
   co.flagi=-1; 
   }}}
 
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///                   Astrometry                                            ///
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

   co.Delta2[0]= l.tetE*0.2500467345692345;    
   co.Delta2[1]= s.mus1*0.2500465923465493;  
   co.Delta2[2]= s.mus2*0.25000465936443; 
   co.Delta2[3]= 0.25000375423543264*l.piE; 
   for(int j=0; j<Ny; ++j){
   for(int k=0; k<Ny; ++k){
   co.inputB[j][k]=0.0; co.inverB[j][k]=0.0; 
   co.adjunB[j][k]=0.0; co.tempoB[j][k]=0.0; co.summ[j][k]=0.0;}}
   
   for(int i=0; i<ndw; ++i){
   for(int j=0; j<Ny;  ++j){
   
   for(int h=0; h<2; ++h)  {
   if(j==0){co.diff=double(co.Delta2[j]*sig2[h]); l.tetE+=co.diff;}
   if(j==1){co.diff=double(co.Delta2[j]*sig[h]);  s.mus1+=co.diff;}
   if(j==2){co.diff=double(co.Delta2[j]*sig[h]);  s.mus2+=co.diff;}
   if(j==3){co.diff=double(co.Delta2[j]*sig2[h]); l.piE +=co.diff;}
   
   
   lightcurve(s,l,as, l.timn[i]);
   co.dera1[h]= double(s.pos1c-l.soux[i])/co.diff;  
   co.derb1[h]= double(s.pos2c-l.souy[i])/co.diff;
   
   if(l.tetE<=0.0 or l.piE<=0.0 or co.diff==0.0 or (s.pos1c==l.soux[i] and s.pos2c==l.souy[i]) or l.erra[i]<=0.0){
   cout<<"BIG ERROR (3), tetE: "<<l.tetE<<"\t l.piE: "<<l.piE<<"\t diff:  "<<co.diff<<"\t l.erra[i]: "<<l.erra[i]<<endl;
   cout<<"pos1c:  "<<s.pos1c<<"\t soux:  "<<l.soux[i]<<"\t pos2c:  "<<s.pos2c<<"\t souy: "<<l.souy[i]<<endl;  cin>>tev; }
   
   if(j==0) l.tetE -=co.diff;
   if(j==1) s.mus1 -=co.diff; 
   if(j==2) s.mus2 -=co.diff;
   if(j==3)  l.piE -=co.diff;}
   
   co.dera1f=(co.dera1[0]+co.dera1[1])*0.5;  
   co.derb1f=(co.derb1[0]+co.derb1[1])*0.5;  
  
  
   for(int k=0; k<=j; ++k){
      
   for(int h=0; h<2; ++h){
   if(k==0){co.diff=double(co.Delta2[k]*sig2[h]);  l.tetE+=co.diff;}
   if(k==1){co.diff=double(co.Delta2[k]*sig[h] );  s.mus1+=co.diff;}
   if(k==2){co.diff=double(co.Delta2[k]*sig[h] );  s.mus2+=co.diff;}
   if(k==3){co.diff=double(co.Delta2[k]*sig2[h]);  l.piE +=co.diff;}
   
   lightcurve(s,l,as, l.timn[i]);
   co.dera2[h]= double(s.pos1c-l.soux[i])/co.diff;  
   co.derb2[h]= double(s.pos2c-l.souy[i])/co.diff; 
   
   if(l.tetE<=0.0 or l.piE<=0.0 or co.diff==0.0 or (s.pos1c==l.soux[i] and s.pos2c==l.souy[i]) ){
   cout<<"BIG ERROR (4), tetE: "<<l.tetE<<"\t l.piE: "<<l.piE<<"\t diff: "<<co.diff<<endl;
   cout<<"pos1c:  "<<s.pos1c<<"\t soux:  "<<l.soux[i]<<"\t pos2c:  "<<s.pos2c<<"\t souy: "<<l.souy[i]<<endl;  cin>>tev; }
   
   if(k==0) l.tetE -=co.diff;
   if(k==1) s.mus1 -=co.diff; 
   if(k==2) s.mus2 -=co.diff;
   if(k==3)  l.piE -=co.diff;}
   
   co.dera2f=(co.dera2[0]+ co.dera2[1])*0.5;  
   co.derb2f=(co.derb2[0]+ co.derb2[1])*0.5;  
   co.inputB[j][k] += (co.dera1f*co.dera2f + co.derb1f*co.derb2f)/(l.erra[i]*l.erra[i]*2.0);}
  
   }//end of loop J
   }//end of loop data
   for(int j=0; j<Ny; ++j){
   for(int k=(j+1);k<Ny; ++k) co.inputB[j][k]= co.inputB[k][j];}



   //cout<<"Covariance matrix Astrometry"<<endl;  
   for(int i=0; i<Ny; ++i){ 
   for(int j=0; j<Ny; ++j){
   if(co.inputB[i][j]==0.0) co.inputB[i][j]=RandR(0.0,1.0)*pow(10.0,-10.0);}}
   //}  cout<<"\n"; }
   //cout<<"*************Input*******************"<<endl;
   //for(int i=0; i<Nx; ++i){ for(int j=0; j<Nx; ++j){cout<<co.inverA[i][j]<<"\t   ";}  cout<<"\n"; } 
   //cout<<"************Inverse*****************"<<endl;  */
   
   inverse(co,1);   
   
   
   for(int i=0; i<Ny; ++i){
   for(int j=0; j<Ny; ++j){
   co.summ[i][j]=0.0;  
   for(int k=0; k<Ny; ++k){co.summ[i][j]+=co.inputB[i][k] *co.inverB[k][j];}}}
   ///if(co.summ[i][j]<1.0e-8) co.summ[i][j]=0.0;  }}
   //cout<<co.summ[i][j]<<"\t";} cout<<"\n";}
   //cout<<"************************************"<<endl;
   
   for(int i=0; i<Ny; ++i){
   for(int j=0; j<Ny; ++j){
   if( (i==j and fabs(co.summ[i][j]-1.0)>0.2) or (i!=j and fabs(co.summ[i][j]-0.0)>0.2)){
   //cout<<"Error_Inverse sum:  "<<co.summ[i][j]<<"\t i:  "<<i<<"\t j:  "<<j<<endl; //cin>>tev;
   //cout<<""<<l.t0<<"\t"<<l.u0<<"\t"<<l.tE<<"\t"<<s.xi*M_PI/180.0<<"\t"<<s.fb[tt]<<"\t"<<s.mbs[tt]<<"\t"<<l.piE<<endl; 
   co.flagi=-1; 
   }}}
} 
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Error parameters Calculation               //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&// 
void ErrorCal(covarian & co, lens & l , source & s){

  double corr1;  
  double epiE, etetE;  
  
  for(int k=0;k<Nx; ++k)  co.Era[k]=sqrt(fabs(co.inverA[k][k]));    
  for(int k=0;k<Ny; ++k)  co.Erb[k]=sqrt(fabs(co.inverB[k][k]));      
  corr1=double(co.inverA[1][4]/fabs(co.Era[1]*co.Era[4]));//correlation tE && xi  

  co.resu[0]= double(co.Era[0]/(fabs(l.u0)   +eps));
  co.resu[1]= double(co.Era[1]/(fabs(l.tE)   +eps));  
  co.resu[2]= double(co.Era[2]/(fabs(s.fb[0])+eps));   
  co.resu[3]= double(co.Era[3]/(fabs(l.piE)  +eps));
  co.resu[4]= double(co.Era[4]/(fabs(s.xi)   +eps));
  co.resu[5]= double(co.Erb[0]/(fabs(l.tetE) +eps));
  co.resu[6]= double(co.Erb[1]/(fabs(s.mus1) +eps));
  co.resu[7]= double(co.Erb[2]/(fabs(s.mus2) +eps));
  co.resu[8]= double(co.Erb[3]/(fabs(l.piE)  +eps));
  co.resu[15]=double(co.Era[5]/(fabs(s.ros) +eps));///delta_ros/ros==delta_tetE/tetE
  
  epiE= MIN(co.resu[3], co.resu[8]); 
  etetE=MIN(co.resu[5], co.resu[15]);   

  co.resu[9] =sqrt(epiE*epiE + etetE*etetE);//Lens Mass
  co.resu[10]=fabs(co.resu[9]*(s.Ds-l.Dl)/s.Ds);//Dl 
    
  co.f1=etetE*etetE+co.resu[1]*co.resu[1]+(co.Era[4]*tan(s.xi))*(co.Era[4]*tan(s.xi))-2.0*co.resu[1]*co.Era[4]*tan(s.xi)*corr1;
  co.f2=etetE*etetE+co.resu[1]*co.resu[1]+(co.Era[4]/tan(s.xi))*(co.Era[4]/tan(s.xi))-2.0*co.resu[1]*co.Era[4]/tan(s.xi)*corr1;
  co.sigmul1=sqrt(co.Erb[1]*co.Erb[1] + l.murel*l.murel*cos(s.xi)*cos(s.xi)*fabs(co.f1));
  co.sigmul2=sqrt(co.Erb[2]*co.Erb[2] + l.murel*l.murel*sin(s.xi)*sin(s.xi)*fabs(co.f2));  
  
  co.resu[11]= double(co.sigmul1/(fabs(l.mul1)+eps));//mu_lens_1
  co.resu[12]= double(co.sigmul2/(fabs(l.mul2)+eps));//mu_lens_2
  co.resu[13]= sqrt(co.sigmul1*co.sigmul1+ co.sigmul2*co.sigmul2)/(fabs(l.mul)+eps);//mu_lens
  co.resu[14]= sqrt(co.Erb[1]*co.Erb[1] +  co.Erb[2]*co.Erb[2]  )/(fabs(s.mus)+eps);// mu_source
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Light curves and Astrometry                //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void lightcurve(source & s, lens & l, astromet & as, double timh)
{ 
   double dvex, dvey, Ve_x, tt;
   double pis= double(1.0/s.Ds);//[mas] source parallax 
   double pil= double(1.0/l.Dl);//[mas] lens parallax  
   double int1=0.0, int2=0.0; 
   
   
   for(int ig=0; ig<2; ++ig){
   if(ig==0) tt=timh; 
   if(ig==1) tt=0.0;    
   dvex= +vearth*sin(omegae*tt + M_PI/2.0)/omegae;//km/s
   dvey= -vearth*cos(omegae*tt + M_PI/2.0)/omegae;//km/s
   as.Ve_n1= cos(tetp)*dvex*sin(l.deltao)-dvey*cos(l.deltao);
   Ve_x=    -cos(tetp)*dvex*cos(l.deltao)-dvey*sin(l.deltao);
   as.Ve_n2=-sin(s.FI)*Ve_x+cos(s.FI)*sin(tetp)*dvex;
   if(ig==0){int1 =as.Ve_n1; int2 =as.Ve_n2;}
   if(ig==1){int1-=as.Ve_n1; int2-=as.Ve_n2;}}
   as.ue_n1=int1;//[radian] without dimention 
   as.ue_n2=int2;//[radian] without dimention 
   
      
   s.ux=-l.u0*sin(s.xi) + (timh-l.t0)*cos(s.xi)/l.tE + l.piE*as.ue_n1;//[] +parallax
   s.uy= l.u0*cos(s.xi) + (timh-l.t0)*sin(s.xi)/l.tE + l.piE*as.ue_n2;//[] +parallax
   
   s.ut0=sqrt((s.ux-l.piE*as.ue_n1)*(s.ux-l.piE*as.ue_n1) + (s.uy-l.piE*as.ue_n2)*(s.uy-l.piE*as.ue_n2));
   s.ut= sqrt(s.ux*s.ux+ s.uy*s.uy);  
   
   if(s.ut ==0.0)  s.ut=1.0e-50; 
   if(s.ut0==0.0)  s.ut0=1.0e-50; 
   
   
   s.def1a= (s.ux-l.piE*as.ue_n1)*l.tetE/(s.ut0*s.ut0 + 2.0);//x-deflection, without parallax[mas]
   s.def2a= (s.uy-l.piE*as.ue_n2)*l.tetE/(s.ut0*s.ut0 + 2.0);//y-deflection, without parallax[mas] 
   
   s.def1c=  s.ux*l.tetE/(s.ut*s.ut + 2.0); //x-deflection[mas]
   s.def2c=  s.uy*l.tetE/(s.ut*s.ut + 2.0); //y-deflection[mas] 
  
   s.pos1b= -l.u0*l.tetE*sin(s.xi) + s.mus1*(timh-l.t0) - as.ue_n1*pis;//x-source trajectory+parallax[mas]
   s.pos2b= +l.u0*l.tetE*cos(s.xi) + s.mus2*(timh-l.t0) - as.ue_n2*pis;//y-source trajectory+parallax[mas]  
  
   s.pos1c= -l.u0*l.tetE*sin(s.xi) + s.mus1*(timh-l.t0) - as.ue_n1*pis +s.def1c;//x-source trajectory+parallax[mas]+lensing
   s.pos2c= +l.u0*l.tetE*cos(s.xi) + s.mus2*(timh-l.t0) - as.ue_n2*pis +s.def2c;//y-source trajectory+parallax[mas]+lensing      
   
   l.pos1= l.mul1*(timh-l.t0) - as.ue_n1*pil ;//x-lens trajectory && parallax[mas]
   l.pos2= l.mul2*(timh-l.t0) - as.ue_n2*pil ;//y-lens trajectory && parallax[mas] 
   
   //if(s.ut<double(5.0*s.ros) or s.ut0<double(5.0*s.ros)){
   //cout<<"Error finite ros: "<<s.ros<<"\t ut:  "<<s.ut<<"\t ut0:  "<<s.ut0<<endl;  
   //int uue; cin>>uue; 
  //}
   
         
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Error LSST calculations                    //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
double errlsstM(double maga, int fi, double sig5){//LSST Photometric Error

   double x,Delta1=0.0;    
   x=pow(10.0 , float(0.4*(maga-sig5)) );
   Delta1=double(sqrt(fabs((0.04-gama[fi])*x+gama[fi]*x*x))); 
   //if(Delta1<0.0001)   Delta1=0.0001; 
   if(sig5<0.0 or sig5>40.0  or Delta1<0.00001 or x<=0.0){
   cout<<"Error errLsstM:  sig5:    "<<sig5<<"\t Delta1:  "<<Delta1<<"\t maga:  "<<maga<<endl;
   int uue; cin>>uue;}
   return( sqrt(delta2*delta2+Delta1*Delta1) );

}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double errlsstA(lsst & ls, double ghadr){//LSST Astrometric Error
    
    double error=-1.0, shib=0.0; 
    
    if(ghadr<ls.maga[0] or  ghadr==ls.maga[0])  error=double(ls.erra[0]); 
    
    else if(ghadr>ls.maga[Na-1] or ghadr==ls.maga[Na-1]){
    shib=(ls.erra[Na-1]-ls.erra[Na-2])/(ls.maga[Na-1]-ls.maga[Na-2]); 
    error=double(ls.erra[Na-1]+shib*(ghadr-ls.maga[Na-1]));}   
    
    else{
    for(int i=1; i<Na; ++i){
    if(double((ghadr-ls.maga[i])*(ghadr-ls.maga[i-1]))<0.0 or ghadr==ls.maga[i-1]){
    shib=(ls.erra[i]-ls.erra[i-1])/(ls.maga[i]-ls.maga[i-1]);
    error=double(ls.erra[i-1]+shib*(ghadr-ls.maga[i-1])); 
    break;}}}
    
    if(error<=0.0 or error<ls.erra[0] or ghadr<0.0 or Na>6){
    cout<<"Error errLSSTA: ghadr:  "<<ghadr<<"\t error: "<<error<<endl;  int dds; cin>>dds;} 
    
    return(error);     
}     
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double errlsstA2(double deltam, double sigast, double sigmag){//LSST Astrometric ability to resolve two close stars
   
    double snr=double(1.0/sigmag);  
    double Factor;  

    if(snr<=5.0)                     Factor=5.0; 
    else if(snr>=5.0 and snr<10.0)   Factor=5.0+(2.0/5.0)*(snr-5.0);
    else if(snr>=10.0 and snr<100.0) Factor=7.0+(13.0/90.0)*(snr-10.0); 
    else                             Factor=20.0+(80.0/900.0)*(snr-100.0);
    
    if(deltam<3.0 and deltam>=2.0)   Factor=1.4*Factor;//Deltam~2.0
    if(deltam>=3.0)                  Factor=3.0*Factor;//Deltam~3

    return( double(Factor*sigast) );  

}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Error ELT  calculations                    //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
/*double errELT(lsst & ls, double ghadr, int flag){
     double error=-1.0, shib=0.0; 
     
     if(flag==0){//photometry-ELT
     if(ghadr<ls.mage[0] or ghadr==ls.mage[0])   error= ls.ere1[0];
     
     else if(ghadr>ls.mage[Nel-1] or ghadr==ls.mage[Nel-1]){
     shib=double(ls.ere1[Nel-1]-ls.ere1[Nel-2])/(ls.mage[Nel-1]-ls.mage[Nel-2]); 
     error=ls.ere1[Nel-1]+shib*(ghadr-ls.mage[Nel-1]);     }
     
     else{
     for(int i=1; i<Nel; ++i){
     if(double((ghadr-ls.mage[i])*(ghadr-ls.mage[i-1]))<0.0 or  ghadr==ls.mage[i-1]){
     shib=double(ls.ere1[i]-ls.ere1[i-1])/(ls.mage[i]-ls.mage[i-1]); 
     error=ls.ere1[i-1]+shib*(ghadr-ls.mage[i-1]); 
     break;}}} 
     }
     
     
     if(flag==1){//astrometry-ELT
     if(ghadr<ls.mage[0] or  ghadr==ls.mage[0])  error= ls.ere2[0];
      
     else if(ghadr>ls.mage[Nel-1] or ghadr==ls.mage[Nel-1]){
     shib=double(ls.ere2[Nel-1]-ls.ere2[Nel-2])/(ls.mage[Nel-1]-ls.mage[Nel-2]); 
     error=ls.ere2[Nel-1]+shib*(ghadr-ls.mage[Nel-1]); }   

     else{
     for(int i=1; i<Nel; ++i){
     if(double((ghadr-ls.mage[i])*(ghadr-ls.mage[i-1]))<0.0 or ghadr==ls.mage[i-1]){
     shib=double(ls.ere2[i]-ls.ere2[i-1])/(ls.mage[i]-ls.mage[i-1]);
     error=ls.ere2[i-1] +  shib*(ghadr-ls.mage[i-1]); 
     break;}}}}
     
     if(error<=0.0 or flag<0 or flag>1 or ghadr<0.0 or ghadr>35.0){
     cout<<"ERROR ELT error: "<<error<<"\t flag:  "<<flag<<"\t ghadr: "<<ghadr<<endl;  int tte; cin>>tte;}
 
     return(error);     
}
*/
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Optical Depth calculations                 //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void optical_depth(source & s)
{
    double ds=(double)s.nums*step;///kpc
    double CC=4.0*GN*M_PI*ds*ds*pow(10.0,9.0)*Msun/(velocity*velocity*KP);
    double dl,x,dx;
    
    s.od_disk=s.od_ThD=s.od_bulge=s.opt=0.0;
    s.od_halo[0]=s.od_halo[1]=0.0; 
    
    for(int k=1; k<s.nums; ++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo[0] +=  s.rho_halo[k][0]*x*(1.0-x)*dx*CC;
    s.od_halo[1] +=  s.rho_halo[k][1]*x*(1.0-x)*dx*CC;}
    
    s.opt=fabs(s.od_disk+ s.od_ThD+ s.od_bulge +s.od_halo[0] +s.od_halo[1]);///total
    //cout<<"total_opticalD: "<<s.opt<<"\t od_disk: "<<s.od_disk<<endl;
    //cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH   
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Func source calculations                   //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void func_source(source & s, CMD & cm, extinc & ex){

    int struc,nums,num;
    double rho,rf;
    double Ds,Av, Alv;
    double Map[M+1],Mab[M+1],Ai[M+1];        
    double test,frand;
    double maxnb=0.0;
    double Avk, Avks, Akv; 
    
    
    for(int i=0; i<(M+1); ++i){
    s.nsbl[i]=  fabs(s.Nstart*pow(FWHM[i]*0.5,2)*M_PI/(3600.0*3600.0));
    s.nsbl[i]+= RandN(sqrt(s.nsbl[i]),2.0);
    if(s.nsbl[i]<=1.0)   s.nsbl[i]=1.0;
    if(s.nsbl[i]>maxnb)  maxnb=s.nsbl[i];}
    
    
    for(int i=0; i<int(M+1); ++i){
    s.Fluxb[i]=0.0; Ai[i]=0.0;  Mab[i]=0.0; Map[i]=0.0; 
    s.magb[i]=0.0; s.Ai[i]=0.0; s.Map[i]=0.0; s.Mab[i]=0.0;}
   
   
  
    for(int k=1; k<=int(maxnb+eps); ++k){  
    
    
    do{
    nums=int(RandR(5.0 , Num-2.0));
    rho=RandR(s.Romins[0] , s.Romaxs[0]); 
    Ds=double(nums*step);
    }while(rho>s.Rostari[nums][0] or Ds<0.0 or Ds>MaxD);///distance larger than 50.0
    //cout<<"k:  "<<k<<"\t Ds:  "<<Ds<<"\t nums:  "<<nums<<endl;
    
    struc=0; 
    rf=RandR(0.0, s.Rostar0[nums][0]); 
     if (rf<= s.rho_disk[nums])                    struc=0;///thin disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums][0])) struc=3;///halo
   // cout<<"struc :  "<<struc<<endl; 
 

    if(struc==0){///thin disk
    num=int(RandR(0.0 , N1-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_d[i][num];}
    Mab[M]=cm.Mk_d[num];
    if(k==1){ 
    s.Rstar=cm.radi_d[num]; 
    s.mass= cm.mass_d[num]; 
    s.age =float(cm.age_d[num]*1.0); 
    s.logT =cm.logT_d[num];
    s.cl =    cm.cl_d[num];
    s.typ=   cm.typ_d[num];}}


    if(struc==1 ){///bulge
    num=int(RandR(0.0,N2-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_b[i][num];}
    Mab[M]=cm.Mk_b[num];
    if(k==1){ 
    s.Rstar=cm.radi_b[num]; 
    s.mass= cm.mass_b[num];
    s.age=float(cm.age_b[num]*1.0); 
    s.logT= cm.logT_b[num];
    s.cl=     cm.cl_b[num];
    s.typ=   cm.typ_b[num];}}


    if(struc==2){///thick disk
    num=int(RandR(0.0,N3-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_t[i][num];}
    Mab[M]=cm.Mk_t[num];
    if(k==1){ 
    s.Rstar=cm.radi_t[num]; 
    s.mass= cm.mass_t[num];
    s.age=float(cm.age_t[num]*1.0); 
    s.logT= cm.logT_t[num];
    s.cl=     cm.cl_t[num];
    s.typ=   cm.typ_t[num];}}


    if(struc==3){///stellar halo
    num=int(RandR(0.0,N4-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_h[i][num];}
    Mab[M]=cm.Mk_h[num];
    if(k==1){ 
    s.Rstar=cm.radi_h[num]; 
    s.mass= cm.mass_h[num];
    s.age=float(cm.age_h[num]*1.0); 
    s.logT= cm.logT_h[num];
    s.cl=     cm.cl_h[num]; 
    s.typ=   cm.typ_h[num];}}
    
    
    
    
    
    double EBV=Interpol(Ds,ex)+eps;///extinction in Ks-band   ???
    //Akv=fabs(a0[M]+ b0[M]/Rv[struc] ) ;///Cardelli et al. 1989 A(k-band)/A(v-band)
    //Avk=1.0/Akv; ////A(v)/A(k)
    //Avks=Avk*0.95;///Marshal 2006 A(v)/A(ks)= Av/Ak Ak/Aks  [8, 9]
    //Av=double(ex.Aks*Avks);
    
    Av=fabs(EBV*Rv[struc]+RandN(sigma[2],0.5));   
    if(Rv[struc]<=0.0  or  Av<0.0  or  Av==0.0 or struc>3){
    cout<<"Error(1) Avks:  "<<Avks<<"\t Av:  "<<Av<<"\t Ds:  "<<Ds<<"\t Aks:  "<<ex.Aks<<endl;  int iis;   cin>>iis; }
   
    
    /*
    extI=-3.0; 
    Alv=fabs( a0[M+1]+ b0[M+1]/Rv[struc]);///AIc/AV
    extI=double(Av*Alv+RandN(sigma[M+1],1) ); //extinction in Ic-band
    if(extI<0.0) extI=0.0;
    if(Av<0.0 or Av==0.0 or extI<0.0 or fabs(Alv-0.6)>0.1){
    cout<<"Error(2) Av:  "<<Av<<"\t extI:   "<<Av*Alv<<"\t Alv:  "<<Alv<<"\t Ds:   "<<Ds<<endl;  int uue;  cin>>uue;}
    */


 
    for(int i=0; i<(M+1); ++i){
    Alv=fabs(a0[i]+b0[i]/Rv[struc]);///Alambda/AV
    Ai[i]= Av*Alv+RandN(sigma[i],1); //extinction in other bands
    if(Ai[i]<0.0) Ai[i]=0.0; 
    Map[i]=double(Mab[i]+5.0*log10(Ds*100.0)+Ai[i]);
    if(s.nsbl[i]>=k) s.Fluxb[i]+=pow(10.0,-0.4*Map[i]);
    //cout<<"filter:  "<<i<<"\t Alv:  "<<Alv<<"\t Ext:  "<<Ai[i]<<"\t Mab:  "<<Mab[i]<<"\t Map:  "<<Map[i]<<endl;
    }
    //cout<<"*****************************"<<endl;  

   
    if(k==1){
    //cout<<"mass:  "<<s.mass<<"\t age:  "<<s.age<<"\t type:  "<<s.typ<<endl;
   // cout<<"cl:  "<<s.cl<<"\t age:  "<<s.age<<"\t tef:  "<<s.logT<<endl;
    s.struc=struc;   
    s.Ds=Ds;  
    s.nums=nums;
    s.Av=Av; 
    for(int i=0; i<(M+1); ++i){
    s.Ai[i] = Ai[i];   
    s.Map[i]=Map[i];   
    s.Mab[i]=Mab[i];}}    
    //cout<<"=========================="<<endl; 
    
    }///loop over the stars 


    for(int i=0; i<int(M+1); ++i){
    s.magb[i]=-2.5*log10(s.Fluxb[i]);
    s.blend[i]=double(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);
    if(s.Fluxb[i]<=0.0 or s.nsbl[i]<1.0 or s.blend[i]>1.00001 or s.blend[i]<=0.0 or(s.nsbl[i]==1.0 and s.blend[i]<1.0)or 
    fabs(s.Map[i]-s.Mab[i]-s.Ai[i]-5.0*log10(s.Ds*100.0))>0.1 or s.Ds<0.0 or s.Ds>MaxD or Av<=0.0){ 
    cout<<"nsbl[i]: "<<s.nsbl[M]<<"\t Ds:  "<<s.Ds<<"\t Av: "<<Av<<endl; 
    cout<<"BIg ERROR nsbl: "<<s.nsbl[i]<<"\t nsbl[i]: "<<s.nsbl[i]<<"\t blend: "<<s.blend[i]<<endl; 
    int uue; cin>>uue;}}
    
    
    s.fb[0]=s.blend[2];//r-LSST 
    s.fb[1]=s.blend[6];//K-band 
    s.mbs[0]=s.magb[2];//r-LSST      
    s.mbs[1]=s.magb[6];//K-band 
    //cout<<">>>>>>>>>>>>  End of Func_source <<<<<<<<<<<<<<<<<<<<<"<<endl;
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Func lens  calculations                    //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void func_lens(lens & l, source & s){
    
    double f,test, tt, Am, DD;
    double randflag; 
    l.rhomaxl=0.0;
     
    for(int k=1; k<int(s.nums-1); ++k){
    l.Dl=double(k*step);
    tt=sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k][1];
    if(tt>l.rhomaxl)  l.rhomaxl=tt;}




    do{
    l.numl=int(RandR(1.0,s.nums-1.0));
    l.Dl  =double(l.numl*step);
    test  =RandR(0.0,l.rhomaxl);
    tt=sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[l.numl][1];
    if(tt>l.rhomaxl or l.Dl<=0 or s.Ds<=0.0 or l.Dl>=s.Ds){
    cout<<"ERROR: rholens[numl]: "<<tt<<""<<l.rhomaxl<<endl;   cout<<"Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;int ue; cin>>ue;}
    }while(test>tt);
    l.Dl=double(l.numl*step);
    //cout<<"Dl :  "<<l.Dl<<"\t numl:  "<<l.numl<<endl;    
   
   
   
   
   
     randflag=RandR(0.0,s.Rostar0[l.numl][1]);
     if(randflag<=(s.rho_disk[l.numl]) ) l.struc=0;//disk
else if(randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1;//bulge
else if(randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2;//thick
else if(randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl][0]))l.struc=3;//halo
else if(randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl][0]+s.rho_halo[l.numl][1]))l.struc=4;
else {cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl][1]<<endl; int ye;  cin>>ye;}



    test= RandR(0.0,1000.01);  
    if(test<898.0){//main-sequence
    l.catg=2;  
    if(l.struc==0){///thin disk
    do{
    l.Ml=RandR(0.08,4.5);
    test=RandR(0.0,57.0);
    if(l.Ml<=1.0) f=pow(l.Ml,-1.6);
    if(l.Ml>1.0)  f=pow(l.Ml,-3.0);
    }while(test>f);}


    if(l.struc==1){//bulge
    do{
    l.Ml=RandR(0.08, 1.4);
    test=RandR(0.3,379.0);
    f=pow(l.Ml, -2.35);
    }while(test>f);}

    if(l.struc==2 or l.struc==3){///thick_disk  & Halo 
    do{
    l.Ml=RandR(0.08, 1.4);
    test=RandR(0.8, 3.6);
    f=pow(l.Ml,-0.5);
    }while(test>f);}

    if(l.struc==4){///Galactic dark-matter halo from stars
    do{
    l.Ml=RandR(0.08,0.8);
    test=RandR(1.5,122.0);
    f=pow(l.Ml,-1.9);
    }while(test>f);}}
  
  
    else if(test>=898 and test<998.0){//Brown dwarfs 
    l.catg=1;  
    do{
    l.Ml=RandR(13.0,  80.0);
    test=RandR(46.0, 167.0);
    f=pow(l.Ml,-0.7)*1000.0;
    }while(test>f);
    l.Ml=double(l.Ml*Mjupiter/Msun);}
  
  
    else{///Stellar mass black holes 
    l.catg=3;///log uniform mass function   
   // l.Ml=RandR(3.0,100.0);
    do{
    l.Ml=RandR(3.0, 100.0);
    test=RandR(0.01, 1.0/3.0)*100.0;
    f=pow(l.Ml,-1.0)*100.0;
    }while(test>f);
    }
 



    l.xls   =double(l.Dl/s.Ds);
    l.RE    =sqrt(4.0*GN*l.Ml*Msun*s.Ds*KP)/velocity;
    l.RE    =double(l.RE*sqrt(l.xls*(1.0-l.xls)));///meter
    l.tetE  =double(l.RE/AU/l.Dl);//[mas]
    s.ros   =double(s.Rstar*Rsun*l.xls/l.RE); 
    l.pirel =double(1.0/l.Dl- 1.0/s.Ds);//[mas]    
    l.piE   =double(l.pirel/l.tetE);//[]   
  
    l.u0    = RandR(0.0, u0m);
    if(l.u0==0.0) l.u0=1.0e-50;  
  
  
  
    l.t0    = RandR(10.0,Tobs-10.0); 
    l.DeltaT= sqrt(4.0+l.u0*l.u0)*l.tetE;//[mas]
    vrel(s,l );
    l.tE=fabs(l.RE/(l.Vt*1000.0*3600.0*24.0));//[days]
    l.A0= (l.u0*l.u0+2.0)/sqrt(l.u0*l.u0*(l.u0*l.u0+4.0)); // 
    l.mi1=s.magb[2]-2.5*log10(fabs(l.A0+1.0)*0.5*s.blend[2] + 1.0 - s.blend[2]);
    l.mi2=s.magb[2]-2.5*log10(fabs(l.A0-1.0)*0.5*s.blend[2] + 1.0 - s.blend[2]);   
    Am=l.A0+1.0; 
    DD=double(4.0-Am*Am + Am*sqrt(Am*Am-4.0))/(Am*Am*0.5-2.0); 
    s.FWHM=2.0*l.tE*sqrt(fabs(DD-l.u0*l.u0));//days 
    // if(DD<double(l.u0*l.u0))  s.FWHM=l.tE;  
    //cout<<"Ml:  "<<l.Ml<<"\t RE:  "<<l.RE/AU<<"\t tE:  "<<l.tE<<endl;
    //cout<<"ros:  "<<s.ros<<"\t murel:  "<<l.murel<<"\t tetE:  "<<l.tetE<<endl;
    //cout<<"u0:  "<<l.u0<<"\t pirel:  "<<l.pirel<<"\t piE:  "<<l.piE<<endl;
    //cout<<"mus1:  "<<s.mus1<<"\t mus2:  "<<s.mus2<<"\t mul1:  "<<l.mul1<<"\t mul2:  "<<l.mul2<<endl;
    //cout<<"DD:  "<<DD<<"\t FWHM:  "<<s.FWHM<<endl;
 

    if(l.tE<=0.0 or l.Dl>s.Ds or l.Vt<=0.0 or l.t0>Tobs or l.t0<0.0 or l.Ml>100.0 or l.Ml<0.0 or l.pirel<0.0 or l.pirel==0.0 or
    s.ros<0.0 or l.RE<0.0 or l.xls>=1.0 or l.Dl<0.0 or s.Ds<0.0 or s.Ds>MaxD or l.u0<0.0 or s.FWHM<0.0 or s.ros<=0.0){
    cout<<"Error_(fun_lens)   Vt:  "<<l.Vt<<"\t ros:  "<<s.ros<<endl;
    cout<<"Error Am:  "<<Am<<"\t DD:  "<<DD<<"\t FWHM:  "<<s.FWHM<<endl;
    cout<<"RE: "<<l.RE/AU<<"\t xls:  "<<l.xls<<"\t tE: "<<l.tE<<endl;
    cout<<"Ml:  "<<l.Ml<<"\t Dl:  "<<l.Dl<<"\t Ds:  "<<s.Ds<<"\t pirel:  "<<l.pirel<<endl;
    cout<<"numl:  "<<l.numl<<"\t Vt:  "<<l.Vt<<"\t murel:  "<<l.murel<<endl;  
    cout<<"Dl:  "<<l.Dl<<"\t Ds:   "<<s.Ds<<"\t xls:  "<<l.xls<<"\t Ml:   "<<l.Ml<<endl;
    cout<<"BIG ERROR te: "<<l.tE<<"\t RE: "<<l.RE<<"\t V_rel: "<<l.Vt<<"\t u0: "<<l.u0<<endl; 
    int iie; cin>>iie;}
  //cout<<">>>>>>>>>>>> End of func_lens <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
}
//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Read CMD                                   //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void read_cmd(CMD & cm ,  yfilter & yf)
{
//  mass  logT    Age  logL logG  Z  u_LSST  g_LSST  r_LSST i_LSST z_LSST  B  V  R    I   J    H    K    CL    Type//1402/04/06
//  0      1       2     3    4   5   6       7       8      9      10     11 12  13  14  15   16   17   18     19
    int yye;     
    double logL, gravity,metal,B,V,R,I,J, H; 
    char filename[40];
    FILE *fp2;
   
////=================================== THIN DISK ============================== 
    int j=0; 
    sprintf(filename,"./files/CMDBESANCON14020406/%c%c%c%c%c%c%d.dat",'C','M','D','T','i','b',3);///changed 1402/04/06
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTi.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    
    fscanf(fp2,"%lf %lf %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %lf %d %lf\n",
    &cm.mass_d[j],&cm.logT_d[j],&cm.age_d[j],&logL,&gravity,&metal,
    &cm.Mab_d[0][j],&cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.Mab_d[3][j],&cm.Mab_d[4][j],&B,&V,&R,&I,&J,&H,&cm.Mk_d[j],&cm.cl_d[j],&cm.typ_d[j]);
    
    gravity=pow(10.0,gravity)*0.01;//[m/s^2] 
    cm.radi_d[j]= sqrt(GN*cm.mass_d[j]*Msun/gravity)/Rsun; ///in sun radius [ Rsun]
    

    //cout<<"metal:   "<<metal<<"\t age:  "<<age<<"\t z_LSST:   "<<cm.Mab_d[4][j]<<endl;
    cm.Mab_d[5][j]=ylsst(yf, metal, float(cm.age_d[j]*1.00), cm.Mab_d[4][j]);
    if(cm.mass_d[j]<0.0 or cm.radi_d[j]<0.0 or cm.radi_d[j]>10000.0 or cm.logT_d[j]<0.0 or cm.Mab_d[2][j]>20.0 or metal>0.12 or cm.age_d[j]>10 or cm.cl_d[j]>7 or cm.typ_d[j]>9.0){
    cout<<"radius "<<cm.radi_d[j]<<endl;
    cout<<"ERROR(cmd thin_disk) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} 
    fclose(fp2);
    if(j!=N1){cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}




////=================================== BULGE ================================== 
    j=0; 
    sprintf(filename,"./files/CMDBESANCON14020406/%c%c%c%c%c%d.dat",'C','M','D','B','b',3);
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDB.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %lf %d %lf\n",
    &cm.mass_b[j],&cm.logT_b[j],&cm.age_b[j],&logL,&gravity,&metal,
    &cm.Mab_b[0][j],&cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.Mab_b[3][j],&cm.Mab_b[4][j],&B,&V,&R,&I,&J,&H,&cm.Mk_b[j],&cm.cl_b[j],&cm.typ_b[j]);
    cm.Mab_b[5][j]= ylsst(yf, metal, float(cm.age_b[j]*1.00), cm.Mab_b[4][j]);
    
    gravity=pow(10.0,gravity)*0.01;//[m/s^2] 
    cm.radi_b[j]= sqrt(GN*cm.mass_b[j]*Msun/gravity)/Rsun; ///in sun radius [ Rsun]
    
    
    if(cm.mass_b[j]<0.0 or cm.logT_b[j]<0.0 or cm.radi_b[j]<0.0 or cm.radi_b[j]>1000.0 or cm.Mab_b[2][j]>18.0 or metal>2.0 or cm.age_b[j]>10 or cm.cl_b[j]>7 or cm.typ_b[j]>9.0){
        cout<<"radius "<<cm.radi_b[j]<<endl;
    cout<<"ERROR(cmd bulge) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} 
    fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N2<<endl;  cin>>yye;}


////=================================== THICK DISK =============================
    j=0; 
    sprintf(filename,"./files/CMDBESANCON14020406/%c%c%c%c%c%c%d.dat",'C','M','D','T','k','b',3);
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTk.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %lf %d %lf\n",
    &cm.mass_t[j],&cm.logT_t[j],&cm.age_t[j],&logL,&gravity,&metal,
    &cm.Mab_t[0][j],&cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.Mab_t[3][j],&cm.Mab_t[4][j],&B,&V,&R,&I,&J,&H,&cm.Mk_t[j],&cm.cl_t[j],&cm.typ_t[j]);
    cm.Mab_t[5][j]= ylsst( yf, metal, float(cm.age_t[j]*1.00), cm.Mab_t[4][j]);
   
    gravity=pow(10.0,gravity)*0.01;//[m/s^2] 
    cm.radi_t[j]=sqrt(GN*cm.mass_t[j]*Msun/gravity)/Rsun; ///in sun radius [ Rsun]
   
    if(cm.mass_t[j]<0.0 or cm.logT_t[j]<0.0 or cm.radi_t[j]<0.0  or cm.radi_t[j]>1000.0 or cm.Mab_t[2][j]>20.0 or metal>0.03 or cm.age_t[j]>8 or cm.cl_t[j]>7 or cm.typ_t[j]>9.0){
            cout<<"radius "<<cm.radi_t[j]<<endl;
    cout<<"ERROR(cmd Thick disk) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} 
    fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N3<<endl;  cin>>yye;}


////=================================== STELLAR HALO =========================== 
    j=0; 
    sprintf(filename,"./files/CMDBESANCON14020406/%c%c%c%c%c%d.dat",'C','M','D','H','b',3);
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDH.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.logT_h[j],&cm.age_h[j],&logL,&gravity,&metal,
    &cm.Mab_h[0][j],&cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.Mab_h[3][j],&cm.Mab_h[4][j],&B,&V,&R,&I,&J,&H,&cm.Mk_h[j],&cm.cl_h[j],&cm.typ_h[j]);
    cm.Mab_h[5][j]= ylsst( yf, metal, float(cm.age_h[j]*1.00), cm.Mab_h[4][j]);
    
    gravity=pow(10.0,gravity)*0.01;//[m/s^2] 
    cm.radi_h[j]=sqrt(GN*cm.mass_h[j]*Msun/gravity)/Rsun; ///in sun radius [ Rsun]
    
    
    if(cm.mass_h[j]<0.0 or cm.logT_h[j]<0.0 or cm.radi_h[j]<0.0  or cm.radi_h[j]>1000.0 or cm.Mab_h[2][j]>20.0 or metal>0.02 or cm.age_h[j]>9 or cm.cl_h[j]>7 or cm.typ_h[j]>9.0){
    cout<<"radius "<<cm.radi_h[j]<<endl;
    cout<<"ERROR(cmd Halo) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} 
    fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N4<<endl;  cin>>yye;}

   cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
   
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double ylsst(yfilter & yf, double metal, double age, double magz)
{
    int h, g, k1,k2; 
    double magy;  
    //cout<<"metal_0:  "<<yf.Metal[0]<<"\t metal_1:  "<<yf.Metal[1]<<"\t metal_3:   "<<yf.Metal[3]<<endl;

    h=-1; 
    if(metal<=yf.Metal[0])          h=0; 
    else if(metal>=yf.Metal[met-1]) h=met-1;
    else{ 
    for(int i=0; i<int(met-1); ++i){
    if(float((metal-yf.Metal[i])*(metal-yf.Metal[i+1]))<0.0 or metal==yf.Metal[i]) {h=i; break;}}}
       
       
    g=-1;    
    k1=int(yf.count[h]);   
    k2=int(yf.count[h]+yf.number[h]);
    if(k1<0)           k1=0;
    if(k2>int(YZ-1))   k2=int(YZ-1);
    
    if(age<=yf.Age[k1])          g=int(k1); 
    else if(age>=yf.Age[k2-1])   g=int(k2-1); 
    
    else{
    for(int k=k1;  k<int(k2-1);  ++k){
    if(float((age-yf.Age[k])*(age-yf.Age[k+1]))<0.0  or  age==yf.Age[k]){g=k; break;}}}
    magy = double(yf.B[g] + yf.M[g]*magz);   
  
    if(g<0 or g>YZ or g<k1 or g>k2 or h<0 or h>met or k1>k2 or k2>YZ or fabs(age-yf.Age[g])>5.0 ){
    cout<<"ERROR2 g:"<<g<<"\t h:  "<<h<<"\t age:"<<age<<"\t metal: "<<metal<<"\t Age[g]:  "<<yf.Age[g]<<endl;
    cout<<"k1: "<<k1<<"\t k2:  "<<k2<<"\t metal[h]:  "<<yf.Metal[h]<<endl;  int uui;  cin>>uui;}
    return( magy ); 
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Galactic model calculations                //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void Disk_model(source & s, int numt)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nnf=0.4/0.8;
   double mBarre;//stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RAa;
   double xf,yf,zf,rho;
   double fd=1.0;//see the program mass_averaged.cpp.  we do not apply any limitation
   double fb=1.0;//0.657066;////just stars brighter than V=11.5, but we change to consider all stars  
   double fh=1.0;//No limitation 
   double Rdd=2.17;//2.53;///2.17;
   double Rhh=1.33;//1.32;//1.33;
   double xb0, yb0, zb0;//kPC
   double frac=0.05;//fraction of halo in the form of compact objects
   double qq=0.688; 
   double Rd0=1.54;
   double xol, yol, zol, x0, y0, z0; 
   double r0,  pos1, inc1, Rlm; 
   double rb, rscale, rhodmh, scale;  
   double almb, rlm_0b;  
   double alm,  rlm_0;  
   
   char filename[40];
   FILE *filj;
   FILE *fill;
   
        
        
        
   int flagf=0;
   if(numt<10){
   flagf=1; 
   filj=fopen("./files/density/right.txt", "a+"); 
   sprintf(filename,"./files/density/%c%d.dat",'E',numt);
   fill=fopen(filename,"w");}
   
   
   s.Nstart=0.0; 
   for(int j=0; j<2; ++j){
   s.Romaxs[j]=s.Rostart[j]=0.0;   
   s.Romins[j]=10000000000.0;}




for(int i=1;i<Num;++i){
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_ThD[i]=0.0;
   s.Nstari[i]=0.0; 
   for(int j=0; j<2; ++j){s.Rostari[i][j]=s.Rostar0[i][j]=s.rho_halo[i][j]=0.0;}




   x=double(i*step);//kpc
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = x*cos(s.FI)*cos(s.TET)-Dsun;
   Rb=sqrt(xb*xb+yb*yb);



///========== Galactic Thin Disk =================================================================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   rho=rho*1.2;  ///totalmass= 4.25e10
   s.rho_disk[i]=fabs(s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii]);  } ///Msun/pc^3  
///===============================================================================================






///========== Galactic Thick Disk ================================================================
  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=fabs((rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nnf))));
  else             s.rho_ThD[i]=fabs((rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*exp(nnf)*exp(-fabs(zb)/0.8)/(1.0+0.5*nnf));///Msun/pc^3
  s.rho_ThD[i]=s.rho_ThD[i]*2.67;//total_mass=0.8e10
///===============================================================================================




///========== Galactic Stellar Halo (spheroid)====================================================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi<=0.5)  s.rho_halo[i][0]=fabs(1.0*(0.932*0.00001/867.067)*pow(0.5/Dsun,-2.44));
   else           s.rho_halo[i][0]=fabs(1.0*(0.932*0.00001/867.067)*pow(rdi/Dsun,-2.44));///Msun/pc^3 totalmass= 1.2e9
   s.rho_halo[i][0]=s.rho_halo[i][0]*5281.0;///Total_mass=1.2e9
///===============================================================================================





///========== Galactic dark-matter Halo ==========================================================
   rb=sqrt(zb*zb+Rb*Rb); //radial distance from the Galactic center [kpc]
   rscale=2.697; //[kpc] 
   rhodmh=0.1079;//Solar mass /pc^3
   scale=double(8.4*0.001/(rhodmh/(1.0+Dsun*Dsun/rscale/rscale))); 
   s.rho_halo[i][1]=frac*scale*rhodmh/(1.0+pow(rb/rscale,2.0));  
   s.rho_halo[i][1]=s.rho_halo[i][1]*0.25; ///total mass=1.8e11
///===============================================================================================




///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));

   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]=  (fabs(rhoS)+fabs(rhoE));///Msun/pc^3 
  s.rho_bulge[i]=s.rho_bulge[i]*0.45;////total mass= 1.7e10
///==================================================================




s.Rostar0[i][0]=fabs(s.rho_disk[i]+s.rho_ThD[i]+s.rho_bulge[i]+s.rho_halo[i][0]);///[Msun/pc^3]///for choosing source stars 
s.Rostar0[i][1]=s.Rostar0[i][0]+fabs(s.rho_halo[i][1]);///for choosing lens objects+ dark-matter halo 




s.Nstari[i]=BF*((s.rho_disk[i])*fd/0.403445+s.rho_ThD[i]*fh/0.4542+(s.rho_halo[i][0])*fh/0.4542+(s.rho_bulge[i])*fb/0.308571);//[Nt/pc^3] 




for(int j=0; j<2; ++j){
s.Rostari[i][j]=s.Rostar0[i][j]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Msun/deg^2]
s.Rostart[j]  +=s.Rostari[i][j];///[Mt/deg^2]
if(s.Rostari[i][j]>s.Romaxs[j]) s.Romaxs[j]=s.Rostari[i][j];///source selection
if(s.Rostari[i][j]<s.Romins[j]) s.Romins[j]=s.Rostari[i][j];///source selection
}

s.Nstari[i]=s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]
s.Nstart  +=s.Nstari[i];///[Nt/deg^2]


  if(flagf>0)
  fprintf(fill,"%e   %e   %e   %e   %e   %e  %e  %e\n",
  x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i][0], s.rho_halo[i][1], s.Rostar0[i][1],s.Nstari[i]);}
  if(flagf>0){
  fprintf(filj,"%.5lf  %.5lf  %.5lf  %.5lf\n",s.right, s.decli, s.TET*RAa, s.FI*RAa);
  fclose(fill); fclose(filj);}

   //cout<<"Nstart [Nt/deg^2]: "<<s.Nstart<<"\t Ro_star [Mass/deg^2]: "<<s.Rostart<<endl;
   //cout<<"Romaxs:  "<<s.Romaxs<<"\t Romins:  "<<s.Romins<<endl;
   //cout<<">>>>>>>>>>>>>>>>>>>>>> END OF DISK MODLE <<<<<<<<<<<<<<<<<<<<"<<endl;    
   //exit(0);
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double RandN(double sigma, double nnd){
   double rr,f,frand;
   do{
   rr=RandR(-sigma*nnd, sigma*nnd);//[-N sigma:N sigma]
   f =exp(-0.5*rr*rr/sigma/sigma);
   frand=RandR(0.0,1.0);
   }while(frand>f);
   return(rr);
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double RandR(double down, double up){
   double p=(double)rand()/((double)(RAND_MAX)+(double)(1.0));
   return(p*(up-down)+down );
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///==============================================================//
///                                                              //
///                  Linear interpolarion                        //
///                                                              //
///==============================================================//
double Interpol(double ds, extinc & ex)
{
  double F=-1.0;
  if(ds<ex.dis[0])       F=ex.Extks[0];
  else if(ds>=ex.dis[99]) F=ex.Extks[99];
  else{ 
  for(int i=0; i<99; ++i){
  if(ex.dis[i]>=ex.dis[i+1]){
  cout<<"ERROR dis[i]: "<<ex.dis[i]<<"\t disi+1: "<<ex.dis[i+1]<<endl;  int yye; cin>>yye; }
  if(ds>=ex.dis[i] and ds<ex.dis[i+1]){
  F = ex.Extks[i]+(double)(ds-ex.dis[i])*(ex.Extks[i+1]-ex.Extks[i])/(ex.dis[i+1]-ex.dis[i]);
  break;}}}
  if(F==-1.0||F<0.0){cout<<"ERROR big Extinction(ds): "<<F<<"\t ds: "<<ds<<endl; int uut; cin>>uut;}
  return(F);
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Relative Velocity calculations             //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void vrel(source & s, lens & l){

  if (l.Dl==0.0) l.Dl=0.00034735;
  double pi=M_PI;
  double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*l.Dl*cos(s.TET)*cos(s.FI));
  double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*s.Ds*cos(s.TET)*cos(s.FI));
  if(Rlc==0.0) Rlc=0.0000000000034346123;
  if(Rsc==0.0) Rsc=0.000000000004762654134; 
 
 
 
 
  double LVx, SVx, vt2;
  double SVT, SVR, SVZ, LVT, LVR, LVZ;
  double fv, testfv, test, age;
  double  VSunx, vls2, vls1;
  double  tetd ;
  double sigma_R_halo[2], sigma_T_halo[2], sigma_Z_halo[2]; 


  double NN=3.5;
  double sigma_R_Disk,        sigma_T_Disk,       sigma_Z_Disk;
  double sigma_R_DiskL,       sigma_T_DiskL,      sigma_Z_DiskL;
  double sigma_R_DiskS,       sigma_T_DiskS,      sigma_Z_DiskS;
  double sigma_R_TDisk=67.0,  sigma_T_TDisk=51.0, sigma_Z_TDisk=42.0;
  double sigma_R_Bulge=113.0, sigma_T_Bulge=115.0,sigma_Z_Bulge=100.0;
  sigma_R_halo[0]= 131.0, sigma_T_halo[0]=106.0, sigma_Z_halo[0]=85.0;//stellar-origin 
  sigma_R_halo[1]= 120.0, sigma_T_halo[1]=120.0, sigma_Z_halo[1]=120.0;//dark matter   Moniez paper+parallax 
  
  
  double Rho[8]={00.0}; 
  double maxr=0.0;
  for(int i=0; i<8; ++i){ 
  Rho[i]=rho0[i]*corr[i]/d0[i]; maxr+=Rho[i];}

 
   for (int i=0;i<2; ++i){
   test= RandR(0.0,maxr); ///total ages
   if(test<=Rho[0])                                          {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0; age= 0.075;}
   else if(test<=(Rho[0]+Rho[1]))                            {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0; age=0.575;}
   else if(test<=(Rho[0]+Rho[1]+Rho[2]))                     {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;age=1.5;  }
   else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]))              {sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2; age=2.5; }
   else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]))       {sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8; age=4.0; }
   else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5])){sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4; age=6.0; }
   else if(test<=maxr)                                       {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5; age=8.5; }
   if(i==0) {
   sigma_R_DiskS= sigma_R_Disk;
   sigma_T_DiskS= sigma_T_Disk;
   sigma_Z_DiskS= sigma_Z_Disk;}
   if(i==1){
   sigma_R_DiskL= sigma_R_Disk;
   sigma_T_DiskL= sigma_T_Disk;
   sigma_Z_DiskL= sigma_Z_Disk;}}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH   SOURCE
    if(s.struc==0){///Galactic disk   source velocity 
    SVR= RandN(sigma_R_DiskS, NN);
    SVT= RandN(sigma_T_DiskS, NN);
    SVZ= RandN(sigma_Z_DiskS, NN); }

    else if(s.struc==1){///Galactic bulge
    SVR= RandN(sigma_R_Bulge, NN);
    SVT= RandN(sigma_T_Bulge, NN);
    SVZ= RandN(sigma_Z_Bulge, NN); }

    else if(s.struc==2){///thick disk
    SVR= RandN(sigma_R_TDisk, NN);
    SVT= RandN(sigma_T_TDisk, NN);
    SVZ= RandN(sigma_Z_TDisk, NN); }

    else if(s.struc==3){///stellar halo
    SVR= RandN(sigma_R_halo[0], NN);
    SVT= RandN(sigma_T_halo[0], NN);
    SVZ= RandN(sigma_Z_halo[0], NN); }
    
        
    if(s.struc==0 or s.struc==2) SVT =SVT+ vro_sun*(1.00762*pow(Rsc/Dsun,0.0394) + 0.00712);//Thin and Thick diskes
    
    s.vs=sqrt( SVR*SVR + SVT*SVT + SVZ*SVZ );
    
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH LENZ

   if(l.struc==0){///Galactic disk
   LVR= RandN(sigma_R_DiskL, NN) ;
   LVT= RandN(sigma_T_DiskL, NN) ;
   LVZ= RandN(sigma_Z_DiskL, NN) ; }

   else if(l.struc==1){///Galactic bulge
   LVR= RandN(sigma_R_Bulge, NN) ;
   LVT= RandN(sigma_T_Bulge, NN) ;
   LVZ= RandN(sigma_Z_Bulge, NN) ; }

   else if(l.struc==2){///thick disk
   LVR= RandN(sigma_R_TDisk, NN) ;
   LVT= RandN(sigma_T_TDisk, NN) ;
   LVZ= RandN(sigma_Z_TDisk, NN) ; }

   else if(l.struc==3){///stellar halo  Lens 
   LVR= RandN(sigma_R_halo[0], NN);
   LVT= RandN(sigma_T_halo[0], NN);
   LVZ= RandN(sigma_Z_halo[0], NN); }
   
   else if(l.struc==4){///stellar halo  Lens ***********
   LVR= RandN(sigma_R_halo[1], NN);
   LVT= RandN(sigma_T_halo[1], NN);
   LVZ= RandN(sigma_Z_halo[1], NN);}
  
 
   
   if(l.struc==0 or l.struc==2) LVT = LVT + vro_sun*(1.00762*pow(Rlc/Dsun,0.0394)+0.00712);
   l.vl=sqrt(LVT*LVT + LVZ*LVZ + LVR*LVR);
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  BETA  HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    l.betal=0.0; l.betas=0.0;
    tetd= s.TET;
    test= double(l.Dl*cos(s.FI)*sin(tetd)/Rlc);
    if(fabs(test-1.0)<0.00000000001745647563)      l.betal= pi/2.0;
    else if(fabs(test+1.0)<0.0000000000012354369)  l.betal=-pi/2.0;
    else                          l.betal=asin(test);
    
    test= double(s.Ds*cos(s.FI)*sin(tetd)/Rsc); 
    if( fabs(test-1.0)<0.01)     l.betas=pi/2.0;
    else if(fabs(test+1.0)<0.01) l.betas=-pi/2.0;
    else                         l.betas=asin(test);
    
    if(Dsun < fabs(l.Dl*cos(s.FI)*cos(tetd)) and l.betal>=0.0) l.betal= pi-l.betal; 
    if(Dsun < fabs(s.Ds*cos(s.FI)*cos(tetd)) and l.betas>=0.0) l.betas= pi-l.betas; 
    if(Dsun < fabs(l.Dl*cos(s.FI)*cos(tetd)) and l.betal<0.0) l.betal=-(pi-fabs(l.betal)); 
    if(Dsun < fabs(s.Ds*cos(s.FI)*cos(tetd)) and l.betas<0.0) l.betas=-(pi-fabs(l.betas)); 


    if(fabs(l.Dl*cos(s.FI)*sin(tetd))>Rlc or fabs(test)>1.0){
    cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
    cout<<"FI: "<<s.FI<<"\t TET: "<<tetd<<"\t betal:  "<<l.betal<<endl;
    cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<"\t betas:   "<<l.betas<<endl;
    cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(tetd)/Rlc<<"\t sin(s): "<<test<<endl;}
    //int ew; cin>>ew;}
///HHHHHHHHHHHHHHHHHHHHHHHHHH  DELTA   HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH



    if(s.TET>pi)  tetd=s.TET-2.0*pi; 
    l.deltal= pi - fabs(tetd) -fabs(l.betal);
    l.deltas= pi - fabs(tetd) -fabs(l.betas);  
    if(l.betal<0.0)  l.deltal= -1.0*l.deltal;
    if(l.betas<0.0)  l.deltas= -1.0*l.deltas;
    l.deltao= pi-fabs(tetd);
    if(tetd<0.0)  l.deltao=-1.0*l.deltao; 

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    s.SV_n1 =+SVR * sin(l.deltas)- SVT * cos(l.deltas);
    s.LV_n1 =+LVR * sin(l.deltal)- LVT * cos(l.deltal);
    s.VSun_n1=+VSunR*sin(l.deltao)-VSunT*cos(l.deltao);
    
    SVx= -SVR*cos(l.deltas)- SVT*sin(l.deltas);
    LVx= -LVR*cos(l.deltal)- LVT*sin(l.deltal);
    VSunx= -VSunR*cos(l.deltao) -VSunT*sin(l.deltao);
    
    s.SV_n2=-sin(s.FI)*(SVx) + cos(s.FI)*SVZ;
    s.LV_n2=-sin(s.FI)*(LVx) + cos(s.FI)*LVZ;
    s.VSun_n2=-sin(s.FI)*(VSunx)+cos(s.FI)*(VSunZ);
 
    
    vls1= l.xls*s.SV_n1 - s.LV_n1 +(1.0-l.xls)*s.VSun_n1;  ///Source - lens 
    vls2= l.xls*s.SV_n2 - s.LV_n2 +(1.0-l.xls)*s.VSun_n2;  /// Source -lens
    l.Vt=sqrt(fabs( vls1*vls1 + vls2*vls2 ) );
    
    
    s.mus1= double(s.SV_n1- s.VSun_n1)*1000.0*3600.0*24.0/(s.Ds*AU);//[mas/days]
    s.mus2= double(s.SV_n2- s.VSun_n2)*1000.0*3600.0*24.0/(s.Ds*AU);//[mas/days]
    l.mul1= double(s.LV_n1- s.VSun_n1)*1000.0*3600.0*24.0/(l.Dl*AU);//[mas/days]
    l.mul2= double(s.LV_n2- s.VSun_n2)*1000.0*3600.0*24.0/(l.Dl*AU);//[mas/days]
    s.mus =sqrt(s.mus1*s.mus1+s.mus2*s.mus2);//constant 
    l.mul= sqrt(l.mul1*l.mul1+l.mul2*l.mul2);
     
    l.murel= sqrt( (s.mus1-l.mul1)*(s.mus1-l.mul1) + (s.mus2-l.mul2)*(s.mus2-l.mul2) );//[mas/days]  
    vt2=l.murel*l.Dl*AU/(1000.0*24.0*3600.0);//[km/s]

    l.murel1= double(s.mus1-l.mul1);
    l.murel2= double(s.mus2-l.mul2);

    if(vls1>0.0 and vls2>=0.0)         s.xi=atan(fabs(vls2)/fabs(vls1));//in radian
    else if(vls1<=0.0 and vls2>0.0)    s.xi=atan(fabs(vls1)/fabs(vls2)) + M_PI/2.0;
    else if(vls1<0.0 and vls2<=0.0)    s.xi=atan(fabs(vls2)/fabs(vls1)) + M_PI;
    else if(vls1>=0.0 and vls2<0.0)    s.xi=atan(fabs(vls1)/fabs(vls2)) + 3.0*M_PI/2.0;
    else if(s.xi > 2.0*M_PI)           s.xi-= 2.0*M_PI; 
    
       
    if(vls1==0.0 and vls2==0.0 and fabs(vt2-l.Vt)>0.1 ){
    cout<<"Big error both zeros:  "<<vls1<<"\t vls_b:  "<<vls2<<endl;
    cout<<"Error l.Vt:  "<<l.Vt<<"\t   vt2:  "<<vt2<<endl;  int iee;  cin>>iee;}

    if (l.Vt<0.0 or l.Vt>1.0e6){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<s.vs<<endl;
    cout<<"source type:  "<<s.cl<<"\t s.struc:  "<<s.struc<<endl;}
   // int yee; cin>>yee;}
    //cout<<"(vrel_func):   l.deltao:  "<<l.deltao<<"FI: "<<s.FI<<endl;   
}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Inverse matrix    calculations             //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void getCofactorA(double input[Nx][Nx], double temp[Nx][Nx], int p, int q, int n){
    int i = 0, j = 0;
    for(int row=0;  row<n;  row++){
    for(int col=0;  col<n;  col++){
    if (row!=p and col!=q){
    temp[i][j++] = input[row][col];
    if (j==int(n-1)){ j=0;  i++;}
    }
    }
    } 
}
//#############################################################################
void getCofactorB(double input[Ny][Ny], double temp[Ny][Ny], int p, int q, int n){
    int i = 0, j = 0;
    for (int row = 0; row<n; row++) {
    for (int col = 0; col<n; col++) {
    if (row != p && col != q) {
    temp[i][j++] = input[row][col];
    if(j ==int(n-1)){ j=0; i++; }
    }
    }
    } 
}
//#############################################################################
double determinantA(double input[Nx][Nx], int n){
    double deter=0.0;
    double A11, A12, A13, A14, A15, A21, A22, A23, A24, A25, A31, A32, A33, A34, A35, A41, A42, A43, A44, A45; 
    double A51, A52, A53, A54, A55;
    double A[5][5];  
    double detr1, detr2, detr3, detr4, detr5, detr6;     
    
   
    if(n==1)      {deter= double(input[0][0]); }
    else if(n==2) {deter= double(input[0][0]*input[1][1]-input[1][0]*input[0][1]);}
    else if(n==3) {deter= double(input[0][0]*(input[1][1]*input[2][2]-input[1][2]*input[2][1])+
                                 input[0][1]*(input[1][2]*input[2][0]-input[1][0]*input[2][2])+
                                 input[0][2]*(input[1][0]*input[2][1]-input[1][1]*input[2][0]));}
          
    else if(n==4){
    A11=input[0][0];     A12=input[0][1];  A13=input[0][2];  A14=input[0][3];
    A21=input[1][0];     A22=input[1][1];  A23=input[1][2];  A24=input[1][3];
    A31=input[2][0];     A32=input[2][1];  A33=input[2][2];  A34=input[2][3];
    A41=input[3][0];     A42=input[3][1];  A43=input[3][2];  A44=input[3][3];
    deter=double(A11*(A22*A33*A44 + A23*A34*A42 + A24*A32*A43 - A24*A33*A42 - A23*A32*A44 - A22*A34*A43 )-
                 A21*(A12*A33*A44 + A13*A34*A42 + A14*A32*A43 - A14*A33*A42 - A13*A32*A44 - A12*A34*A43 )+
                 A31*(A12*A23*A44 + A13*A24*A42 + A14*A22*A43 - A14*A23*A42 - A13*A22*A44 - A12*A24*A43 )-
                 A41*(A12*A23*A34 + A13*A24*A32 + A14*A22*A33 - A14*A23*A32 - A13*A22*A34 - A12*A24*A33 ));} 

    else if(n==5){//we had to set this because argument of deter55 takes only a 5x5 matrix
    A[0][0]=input[0][0]; A[0][1]=input[0][1];  A[0][2]=input[0][2];  A[0][3]=input[0][3];  A[0][4]=input[0][4];
    A[1][0]=input[1][0]; A[1][1]=input[1][1];  A[1][2]=input[1][2];  A[1][3]=input[1][3];  A[1][4]=input[1][4];
    A[2][0]=input[2][0]; A[2][1]=input[2][1];  A[2][2]=input[2][2];  A[2][3]=input[2][3];  A[2][4]=input[2][4];
    A[3][0]=input[3][0]; A[3][1]=input[3][1];  A[3][2]=input[3][2];  A[3][3]=input[3][3];  A[3][4]=input[3][4];
    A[4][0]=input[4][0]; A[4][1]=input[4][1];  A[4][2]=input[4][2];  A[4][3]=input[4][3];  A[4][4]=input[4][4];
    deter=deter55(A);}
    
    
    else if(n==6){
    A[0][0]=input[1][1]; A[0][1]=input[1][2];  A[0][2]=input[1][3];  A[0][3]=input[1][4]; A[0][4]=input[1][5];
    A[1][0]=input[2][1]; A[1][1]=input[2][2];  A[1][2]=input[2][3];  A[1][3]=input[2][4]; A[1][4]=input[2][5];
    A[2][0]=input[3][1]; A[2][1]=input[3][2];  A[2][2]=input[3][3];  A[2][3]=input[3][4]; A[2][4]=input[3][5];
    A[3][0]=input[4][1]; A[3][1]=input[4][2];  A[3][2]=input[4][3];  A[3][3]=input[4][4]; A[3][4]=input[4][5];
    A[4][0]=input[5][1]; A[4][1]=input[5][2];  A[4][2]=input[5][3];  A[4][3]=input[5][4]; A[4][4]=input[5][5];
    detr1=deter55(A);  
    
    A[0][0]=input[1][0]; A[0][1]=input[1][2];  A[0][2]=input[1][3];  A[0][3]=input[1][4]; A[0][4]=input[1][5];
    A[1][0]=input[2][0]; A[1][1]=input[2][2];  A[1][2]=input[2][3];  A[1][3]=input[2][4]; A[1][4]=input[2][5];
    A[2][0]=input[3][0]; A[2][1]=input[3][2];  A[2][2]=input[3][3];  A[2][3]=input[3][4]; A[2][4]=input[3][5];
    A[3][0]=input[4][0]; A[3][1]=input[4][2];  A[3][2]=input[4][3];  A[3][3]=input[4][4]; A[3][4]=input[4][5];
    A[4][0]=input[5][0]; A[4][1]=input[5][2];  A[4][2]=input[5][3];  A[4][3]=input[5][4]; A[4][4]=input[5][5];
    detr2=deter55(A);  
    
    A[0][0]=input[1][0]; A[0][1]=input[1][1];  A[0][2]=input[1][3];  A[0][3]=input[1][4]; A[0][4]=input[1][5];
    A[1][0]=input[2][0]; A[1][1]=input[2][1];  A[1][2]=input[2][3];  A[1][3]=input[2][4]; A[1][4]=input[2][5];
    A[2][0]=input[3][0]; A[2][1]=input[3][1];  A[2][2]=input[3][3];  A[2][3]=input[3][4]; A[2][4]=input[3][5];
    A[3][0]=input[4][0]; A[3][1]=input[4][1];  A[3][2]=input[4][3];  A[3][3]=input[4][4]; A[3][4]=input[4][5];
    A[4][0]=input[5][0]; A[4][1]=input[5][1];  A[4][2]=input[5][3];  A[4][3]=input[5][4]; A[4][4]=input[5][5];
    detr3=deter55(A); 
    
    A[0][0]=input[1][0]; A[0][1]=input[1][1];  A[0][2]=input[1][2];  A[0][3]=input[1][4]; A[0][4]=input[1][5];
    A[1][0]=input[2][0]; A[1][1]=input[2][1];  A[1][2]=input[2][2];  A[1][3]=input[2][4]; A[1][4]=input[2][5];
    A[2][0]=input[3][0]; A[2][1]=input[3][1];  A[2][2]=input[3][2];  A[2][3]=input[3][4]; A[2][4]=input[3][5];
    A[3][0]=input[4][0]; A[3][1]=input[4][1];  A[3][2]=input[4][2];  A[3][3]=input[4][4]; A[3][4]=input[4][5];
    A[4][0]=input[5][0]; A[4][1]=input[5][1];  A[4][2]=input[5][2];  A[4][3]=input[5][4]; A[4][4]=input[5][5];
    detr4=deter55(A); 
    
    A[0][0]=input[1][0]; A[0][1]=input[1][1];  A[0][2]=input[1][2];  A[0][3]=input[1][3]; A[0][4]=input[1][5];
    A[1][0]=input[2][0]; A[1][1]=input[2][1];  A[1][2]=input[2][2];  A[1][3]=input[2][3]; A[1][4]=input[2][5];
    A[2][0]=input[3][0]; A[2][1]=input[3][1];  A[2][2]=input[3][2];  A[2][3]=input[3][3]; A[2][4]=input[3][5];
    A[3][0]=input[4][0]; A[3][1]=input[4][1];  A[3][2]=input[4][2];  A[3][3]=input[4][3]; A[3][4]=input[4][5];
    A[4][0]=input[5][0]; A[4][1]=input[5][1];  A[4][2]=input[5][2];  A[4][3]=input[5][3]; A[4][4]=input[5][5];
    detr5=deter55(A); 
    
    A[0][0]=input[1][0]; A[0][1]=input[1][1];  A[0][2]=input[1][2];  A[0][3]=input[1][3]; A[0][5]=input[1][4];
    A[1][0]=input[2][0]; A[1][1]=input[2][1];  A[1][2]=input[2][2];  A[1][3]=input[2][3]; A[1][4]=input[2][4];
    A[2][0]=input[3][0]; A[2][1]=input[3][1];  A[2][2]=input[3][2];  A[2][3]=input[3][3]; A[2][4]=input[3][4];
    A[3][0]=input[4][0]; A[3][1]=input[4][1];  A[3][2]=input[4][2];  A[3][3]=input[4][3]; A[3][4]=input[4][4];
    A[4][0]=input[5][0]; A[4][1]=input[5][1];  A[4][2]=input[5][2];  A[4][3]=input[5][3]; A[4][4]=input[5][4];
    detr6= deter55(A);
    deter=double(detr1*input[0][0]-detr2*input[0][1]+detr3*input[0][2]-detr4*input[0][3]+detr5*input[0][4]-detr6*input[0][5]);}      
    else{cout<<"Error n is larger than 6:  "<<n<<"\t deter "<<deter<<endl;  int yye; cin>>yye;} 
    return(deter);
}
//#################################################################################################
double deter55(double AA[5][5]){

    double A11, A12, A13, A14, A15, A21, A22, A23, A24, A25, A31, A32, A33, A34, A35, A41, A42, A43, A44, A45; 
    double A51, A52, A53, A54, A55, deter=0.0;  
    
    A11=AA[0][0];     A12=AA[0][1];  A13=AA[0][2];  A14=AA[0][3]; A15=AA[0][4];
    A21=AA[1][0];     A22=AA[1][1];  A23=AA[1][2];  A24=AA[1][3]; A25=AA[1][4];
    A31=AA[2][0];     A32=AA[2][1];  A33=AA[2][2];  A34=AA[2][3]; A35=AA[2][4];
    A41=AA[3][0];     A42=AA[3][1];  A43=AA[3][2];  A44=AA[3][3]; A45=AA[3][4];
    A51=AA[4][0];     A52=AA[4][1];  A53=AA[4][2];  A54=AA[4][3]; A55=AA[4][4];
    
    //for (int row = 0; row<5; row++) {
    //for (int col = 0; col<5; col++) {
    //cout<<"A:  "<<AA[row][col]<<endl; }}
    
    deter=+A11*A22*(A33*A44*A55 + A34*A45*A53 + A35*A43*A54 - A35*A44*A53 - A34*A43*A55 - A33*A45*A54 )
          -A11*A32*(A23*A44*A55 + A24*A45*A53 + A25*A43*A54 - A25*A44*A53 - A24*A43*A55 - A23*A45*A54 )
          +A11*A42*(A23*A34*A55 + A24*A35*A53 + A25*A33*A54 - A25*A34*A53 - A24*A33*A55 - A23*A35*A54 )
          -A11*A52*(A23*A34*A45 + A24*A35*A43 + A25*A33*A44 - A25*A34*A43 - A24*A33*A45 - A23*A35*A44 )
          
          -A21*A12*(A33*A44*A55 + A34*A45*A53 + A35*A43*A54 - A35*A44*A53 - A34*A43*A55 - A33*A45*A54 )
          +A21*A32*(A13*A44*A55 + A14*A45*A53 + A15*A43*A54 - A15*A44*A53 - A14*A43*A55 - A13*A45*A54 )
          -A21*A42*(A13*A34*A55 + A14*A35*A53 + A15*A33*A54 - A15*A34*A53 - A14*A33*A55 - A13*A35*A54 )
          +A21*A52*(A13*A34*A45 + A14*A35*A43 + A15*A33*A44 - A15*A34*A43 - A14*A33*A45 - A13*A35*A44 )
          
          +A31*A12*(A23*A44*A55 + A24*A45*A53 + A25*A43*A54 - A25*A44*A53 - A24*A43*A55 - A23*A45*A54 )
          -A31*A22*(A13*A44*A55 + A14*A45*A53 + A15*A43*A54 - A15*A44*A53 - A14*A43*A55 - A13*A45*A54 )
          +A31*A42*(A13*A24*A55 + A14*A25*A53 + A15*A23*A54 - A15*A24*A53 - A14*A23*A55 - A13*A25*A54 )
          -A31*A52*(A13*A24*A45 + A14*A25*A43 + A15*A23*A44 - A15*A24*A43 - A14*A23*A45 - A13*A25*A44 )
          
          -A41*A12*(A23*A34*A55 + A24*A35*A53 + A25*A33*A54 - A25*A34*A53 - A24*A33*A55 - A23*A35*A54 )
          +A41*A22*(A13*A34*A55 + A14*A35*A53 + A15*A33*A54 - A15*A34*A53 - A14*A33*A55 - A13*A35*A54 )
          -A41*A32*(A13*A24*A55 + A14*A25*A53 + A15*A23*A54 - A15*A24*A53 - A14*A23*A55 - A13*A25*A54 )
          +A41*A52*(A13*A24*A35 + A14*A25*A33 + A15*A23*A34 - A15*A24*A33 - A14*A23*A35 - A13*A25*A34 )
          
          +A51*A12*(A23*A34*A45 + A24*A35*A43 + A25*A33*A44 - A25*A34*A43 - A24*A33*A45 - A23*A35*A44 )
          -A51*A22*(A13*A34*A45 + A14*A35*A43 + A15*A33*A44 - A15*A34*A43 - A14*A33*A45 - A13*A35*A44 )
          +A51*A32*(A13*A24*A45 + A14*A25*A43 + A15*A23*A44 - A15*A24*A43 - A14*A23*A45 - A13*A25*A44 )
          -A51*A42*(A13*A24*A35 + A14*A25*A33 + A15*A23*A34 - A15*A24*A33 - A14*A23*A35 - A13*A25*A34 );
          
   //cout<<"deter: "<<deter <<endl;        
   return(deter); 
}
//#################################################################################################          
double determinantB(double input[Ny][Ny], int n){
    double deter=0.0;

    
    if(n==1) deter= double(input[0][0]);
    
    else if(n==2) deter= double(input[0][0]*input[1][1]-input[1][0]*input[0][1]);
    
    else if(n==3){
    deter=input[0][0]*(input[1][1]*input[2][2]-input[1][2]*input[2][1])+
          input[0][1]*(input[1][2]*input[2][0]-input[1][0]*input[2][2])+
          input[0][2]*(input[1][0]*input[2][1]-input[1][1]*input[2][0]); } 
    
    
    else if(n==4){
    double A11, A12, A13, A14, A21, A22, A23, A24, A31, A32, A33, A34, A41, A42, A43, A44;  
    A11=input[0][0];     A12=input[0][1];  A13=input[0][2];  A14=input[0][3];
    A21=input[1][0];     A22=input[1][1];  A23=input[1][2];  A24=input[1][3];
    A31=input[2][0];     A32=input[2][1];  A33=input[2][2];  A34=input[2][3];
    A41=input[3][0];     A42=input[3][1];  A43=input[3][2];  A44=input[3][3];
    deter=A11*(A22*A33*A44 + A23*A34*A42 + A24*A32*A43 - A24*A33*A42 - A23*A32*A44 - A22*A34*A43 )-
          A21*(A12*A33*A44 + A13*A34*A42 + A14*A32*A43 - A14*A33*A42 - A13*A32*A44 - A12*A34*A43 )+
          A31*(A12*A23*A44 + A13*A24*A42 + A14*A22*A43 - A14*A23*A42 - A13*A22*A44 - A12*A24*A43 )-
          A41*(A12*A23*A34 + A13*A24*A32 + A14*A22*A33 - A14*A23*A32 - A13*A22*A34 - A12*A24*A33 ); }
     
    else {cout<<"Error n is larger than 4:  "<<n<<"\t deter "<<deter<<endl;  int yye; cin>>yye;} 
     
     
    /*
    double temd[Ny][Ny]; 
    int sign = 1; 
    for(int f=0; f<n; f++){
    getCofactorB(input, temd, 0, f, n);
    deter+=sign * input[0][f]*determinantB(temd, int(n-1));
    sign  =-sign;}*/
        
    return(deter);
}
//#############################################################################
void inverse(covarian & co, int flag){

    co.sign=0; 
    co.deter=0.0; 
     
 
    if(flag==0){//Nx
    //cout<<"Adjuint matrix \n"<<endl;
    //for (int i = 0; i<Nx; i++){
    //for (int j = 0; j<Nx; j++){co.tempoA[i][j]=0.0;  co.adjunA[j][i]
    for(int i = 0; i<Nx; i++){
    for(int j = 0; j<Nx; j++){
    getCofactorA(co.inputA, co.tempoA, i, j, Nx);
    co.sign =((i + j) % 2 == 0) ? 1 : -1;
    co.adjunA[j][i] = co.sign * determinantA(co.tempoA , int(Nx-1) );}}
    
    //for (int i = 0; i<Nx; i++){
    //for (int j = 0; j<Nx; j++){
    //cout<<co.adjunA[i][j]<<"\t   ";//}
    //cout<<"\n"; //}
    
    co.deter=determinantA(co.inputA , Nx);
    //cout<<"Determinate: "<< co.deter <<endl;
    if(co.deter==0.0){cout << "Singular matrix, can't find its inverse." << endl; co.deter=1e-50; }// int uue; cin>>uue;} 
    for(int i=0; i<Nx; i++){
    for(int j=0; j<Nx; j++){co.inverA[i][j] =co.adjunA[i][j]/double(co.deter);}}}
    
    
    
    if(flag==1){//Ny
    //cout<<"Adjuint matrix \n"<<endl;
    for (int i = 0; i<Ny; i++){
    for (int j = 0; j<Ny; j++){
    getCofactorB(co.inputB, co.tempoB, i, j, Ny);
    co.sign =((i + j) % 2 == 0) ? 1 : -1;
    co.adjunB[j][i] = co.sign * determinantB(co.tempoB ,  int(Ny-1) );}}
    for (int i = 0; i<Ny; i++){
    for (int j = 0; j<Ny; j++){
    //cout<<co.adjunB[i][j]<<"\t   ";
    }
    //cout<<"\n";
    }
    co.deter=determinantB(co.inputB , Ny);
    //cout<<"Determinate: "<< co.deter <<endl;
    if (co.deter==0.0){cout << "Singular matrix, can't find its inverse." << endl; co.deter=1e-50; }// int uue; cin>>uue;} 
    for(int i=0; i<Ny; i++){
    for(int j=0; j<Ny; j++){ co.inverB[i][j] =co.adjunB[i][j]/double(co.deter);}} }


}
//###############################################################################################################
//###############################################################################################################
///==============================================================//
///                                                              //
///                  EXtinction                                 //
///                                                              //
///==============================================================//
int Extinction(extinc & ex,source & s)
{
     double sig,Lon,Lat;
     int uue, flag=0;
     
     if(s.lon<0.0){sig=-1.0;cout<<"Strange!!!!longtitude is negative:  s.lon "<<s.lon<<endl; cin>>uue;}
     else sig=1.0;
     double delt=fabs(s.lon)-floor(fabs(s.lon));
     if(delt>1.0 or delt<0.0) {cout<<"ERROR longtitude: delt: "<<delt<<"\t s.lon: "<<s.lon<<endl;  cin>>uue; }
     else if(delt<0.25) Lon=(floor(fabs(s.lon))+0.00)*sig;
     else if(delt<0.50) Lon=(floor(fabs(s.lon))+0.25)*sig;
     else if(delt<0.75) Lon=(floor(fabs(s.lon))+0.50)*sig;
     else               Lon=(floor(fabs(s.lon))+0.75)*sig;
     //if(s.lon==0.0)     Lon=360.00;

     if(s.lat<0.0) sig=-1.0;
     else sig=1.0;
     delt=fabs(s.lat)-floor(fabs(s.lat));
     if(delt>1.0 or delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lon: "<<s.lat<<endl;  cin>>uue;}
     else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
     else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
     else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
     else                Lat=(floor(fabs(s.lat))+0.75)*sig;
     if(Lon>360.000 or Lon<0.25 or fabs(Lat)>20.0 or (Lon>100 and Lon<260)){
     cout<<"BIG error (stopped program) s.lon: "<<Lon<<"\t s.lat: "<<Lat<<endl;   cin>>uue;}

     
     if(Lon==360.0)      Lon=0.00;
     char filename[40];
     FILE *fpd;
     sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',Lat,'_',Lon);
     fpd=fopen(filename,"r");
     cout<<"Extnction file name:  Lat:  "<<Lat<<"\t Lon:    "<<Lon<<endl;
   
    


     double lonti,latit;
     if(!fpd){
     cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
     int yye; cin>>yye;  
     FILE *SD;
     SD=fopen("./files/Ext/saved_direction.txt","r");
     for(int i=0; i<128961; ++i){
     fscanf(SD,"%lf %lf \n",&latit,&lonti);
     if(fabs(Lat-latit)<0.1 and fabs(Lon-lonti)<0.1){
     cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
     cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
     flag=-1;}
     else{
     flag=1;
     for(int i=0; i<100; ++i){
     fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
     //cout<<"distance: "<<ex.dis[i]<<"\t Extks: "<<ex.Extks[i]<<endl;
     if(ex.dis[i]<0.2  or ex.dis[i]>50.0 or ex.Extks[i]<0.0){
     cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
     cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0;  }
    /// if(ex.dis[i]==8.25) ex.exks8= ex.Extks[i]; 
     }}
     //cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<endl;
     fclose(fpd);
     return(flag);
}
//#########################################################################################################################
void sortarray(lsst &ls){

    int sum=0, pos, max;   
    for(int i=0; i<6; i++){
    sum+=ls.Nfilt[i];
    ls.Nfilt2[i]=0;
    ls.argm[i]  =0;}
 
    for(int i=0; i<6; i++){
    max=0;
    pos=0;
    for(int q=0; q<6; q++){
    if(ls.Nfilt[q]>max){
    max=ls.Nfilt[q];  pos=q;}}
    ls.Nfilt2[i]=max;
    ls.argm[i]  =pos;
    cout<<"Sorted_N_filters:  "<<ls.Nfilt2[i]<<"\t its_argument:  "<<ls.argm[i]<<endl;
    ls.Nfilt[pos]=0;}
    
 }
//#########################################################################################################################
