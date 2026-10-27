import numpy as np

################################################################################
Avks=float(8.35816659571)##8.20922) ##for Rv=3.1 for disk stars, because all directions with one array was toward the disk 
distance=np.zeros(100)
extK=np.zeros(100)
disp=np.zeros(100) 
extp=np.zeros(100)
fill=open("saved_direction.txt","w")
fill.close()
sss=np.zeros((1,2))
dff=open("Marshalmap.txt","w")
dff.close()
nno=int(0)
row=int(0)
Rv=0.0;  
with open("3DExtMapMarshal.dat","r") as fd:
    for line in fd:
        row=row+1
        text=line.split(' ')
        NN=len(text) 
        if(int(NN)%2!=0):
            print("NN is ODD !!!", NN)
        ll=int((NN-6)/2)
        if(row>0):
            disp=distance; 
            extp=extK;    
        print ("*************************************")   
        print ("number of distances:  ", ll )
        print ("b: ", text[0] , "l: ", text[1] )
        
        dis=np.zeros(ll);  
        ext=np.zeros(ll)
        for i in range(ll):
            j=int(i*2+4)           
            if(j>=NN):
                print (j, "An error " )
            dis[i]=float(text[j])
            ext[i]=float(text[j+1])
            if(dis[i]<dis[i-1] and i>=1):
                print(">>>>>>>>>>>>>BIG ERROR 1:  dist:", dist, extin)
                print ("Length of row: ", NN, "  number of row: ", row )
                print ("dis, ext: ", dis, ext )
                print (text[0], text[1] )
################################################################################    
        if(ll>2):
            for i in range(100):
                dist=float(0.25+i*0.5)
                extin=2000.0;
                if dist<=float(dis[0]):
                    m=float((ext[1]-ext[0])/(dis[1]-dis[0]));
                    extin=float(m*(dist-dis[0])+ext[0])
                    if(extin<0.0): extin=0.0
                #######################################################    
                elif(dist>float(dis[ll-1]) ):
                    if(dist<15.0 and i<30):### up to distance 15kpc the extinction is increasing 
                        m=(ext[ll-1]-ext[ll-2])/(dis[ll-1]-dis[ll-2])*abs(30.0-i)/30.0
                    else:  ### for dist>15 extinction does not change
                        m= 0.0      
                    extin=float(m*0.5 + extK[i-1])                                    
                    if(extin>0.7 and abs(float(text[0]))>2.5): 
                        extin=extK[i-1]
                         
                #######################################################          
                else:
                    for j in range(ll-1):
                        if(dist>dis[j] and dist<=dis[j+1]):
                            m=float((ext[j+1]-ext[j])/(dis[j+1]-dis[j]))
                            extin=float((dist-dis[j])*m+ext[j])
                            break 
                #######################################################
                if(extin==2000.0):
                    print ("BIG ERROR extin: ", extin )
                    input("Enter a number") 
                if(extin<0.0 or m<0.0):   
                    print ("EXtinction is negative:  ", extin, m, dist, dis[0], dis[ll-1] )
                    print ("b(degree), l(0<l<=360.0, degree): ",  text[0],   text[1] )
                    extin=0.0
                    input("Enter a number ")
                

                    #print "Length of row: ", NN, "  number of row: ", row
                    #print "dis, ext: ", dis, ext
                    #print text[0], text[1]
                    #input("Enter a number ")
                distance[i]=float(dist)
                extK[i]=extin
          
               
################################################################################
        else:
            print("one entry: ",   text[0], text[1] )
            nno+=1
            #input ("Enter a number ")
            '''
            for i in range(100):
                dist=###float(0.25+i*0.5)
                mm=###float(0.7/Avks/1.0);
                extK[i]=float((dist-dis[0])*mm+ext[0])
                distance[i]=float(dist) 
                if(float(extK[i])<0.0):
                    extK[i]=float(0.0)
            '''
            distance= disp; 
            extK= extp; 
################################################################################
        '''
        plt.clf()
        plt.plot(distance,extK,"ro", markersize=2.5) 
        plt.xticks(fontsize=13, rotation=0)
        plt.yticks(fontsize=13, rotation=0)
        plt.xlabel(r"$\rm{distance(kpc)}$",fontsize=16.0)
        plt.ylabel(r"$\rm{Ext(Ks)}$",fontsize=16.0)
        plt.xlim([0.0,50.0])
        plt.grid("True")
        plt.grid(linestyle='dashed')
        if(row%1==0):
            fig3=plt.gcf()
            fig3.savefig("./figs/{0:d}_{1:.2f}_{2:.2f}.jpg".format(row,float(text[0]),float(text[1])),dpi=200)
            py.clf()
        '''            
################################################################################
        '''
        for i in range(100):
            #print distance[i], extK[i]      
            if(float(extK[i]*1.00083647)<extK[i-1] and float(extK[i])!=float(extK[i-1]) and i>0):
                print(">>>>>>>>>>>>>BIG ERROR 2:  dist:", distance ,  extK)
                print ("Length of row: ", NN, "  number of row: ", row )
                print "dis, ext: ", dis, ext
                print text[0], text[1]
                input("Enter a number  2 ")
        '''        
        
        lat, longi=float(text[0]), float(text[1])
        if(abs(lat)<2.5 and (abs(longi)<2.5 or longi>357.5)):  Rv=2.5## Galactic bulge 
        else:                                                  Rv=3.1## Galactic disks  
        extK=np.abs(extK*0.95/(Rv*(0.1615-0.1483/Rv)))##A(Ks)->E(B-V) reddening  2 Mehr 1405  last updated 
                
        filew=open("./ExtMarshal/Ext{0:.2f}_{1:.2f}.dat".format(float(lat),float(longi)),"w") 
        np.savetxt(filew,np.concatenate((distance.reshape((-1,1)),extK.reshape((-1,1))),axis=1),fmt="%.2f  %.4f")
        filew.close()
        
        if(longi==360.0): 
            filew=open("./ExtMarshal/Ext{0:.2f}_{1:.2f}.dat".format(float(lat),float(0.00)),"w") 
            np.savetxt(filew,np.concatenate((distance.reshape((-1,1)),extK.reshape((-1,1))),axis=1),fmt="%.2f  %.4f")
            filew.close() 
        
        sss[0,:]=[float(text[0]),float(text[1])]
        fill=open("saved_direction.txt","a+")
        np.savetxt(fill,sss,fmt='%.2f   %.2f')
        fill.close()
        #input("Enter a number")
################################################################################
        if(int(row)%500==0):
            print ("Step: ", row )
if(row!=int(64881)):
    print (row, "ERROR in the number of rows !!!!", "nno:  ", nno )
print (" ***************** END OF CODE *********************************") 






