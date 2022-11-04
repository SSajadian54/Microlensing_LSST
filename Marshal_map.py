import numpy as np
import matplotlib.pyplot as plt


dff=open("Marshalmap.txt","w")
dff.close()

row=int(0)
point=np.zeros((2))
with open("3D-ExtMap.dat","r") as fd:
    for line in fd:
        row=row+1
        text=line.split(' ')
        #NN=len(text) 
        #print text
        point[0]=float(text[0])
        point[1]=float(text[1])
        dff=open("Marshalmap.txt","a+")
        np.savetxt(dff,point.reshape((1,2)),fmt="%.3f   %.3f")
        #print point 
        dff.close()
        #print "latitude: ", text[0],  "longtitude: ", text[1], "row: ", row
        #print "*********************************************************" 
        #input("Enter a number ")

print "number of rows: ", row 
################################################################################


a=360.0; 

point2=np.zeros((row,2))
point2=np.loadtxt("Marshalmap.txt")
plt.clf()
for i in range(row):
    if (float(point2[i,1])>100.0):
        point2[i,1]=float(point2[i,1]-360.0)

plt.plot(point2[:,1],point2[:,0],'rd',ms=0.05 ,label="Marshal map")
plt.xlabel(r"$Galactic longtitude (deg)$")
plt.ylabel(r"$GAlactic Latitude (deg)$")
fig3=plt.gcf()
fig3.savefig("Marshal map")
print ">>>>>>>>>>> Marshal map was made <<<<<<<<<<<<<<<"


flag=np.zeros((row))
point3=np.zeros((row,2))
point3=np.loadtxt("Marshalmap.txt")

ii=int(-1)
for i in range(80):
    b=float(-10.0+i*0.25)
    for j in range(801):
        if(j<400):
            l= 0.25+j*0.25
        else: 
            l= 260.0+(j-400.0)*0.25
        ii+=1
        test=int(-1)
        for k in range (row):
            if(abs(float(point3[k,0])-b)<0.05 and abs(float(point3[k,1])-l)<0.05):
                if(flag[k]==0):
                    flag[k]=int(ii)
                else:
                    print flag[k], "ERROR !!!!!"
                    print "ERROR latitude: ", b, "longtitude: ", l
                    input ("Enter a number ")
                #print "FOUND latitude: ", b, "longtitude: ", l, point3[k,:], int(flag[k])
                test=1
                break
        if(test==-1):
            print "ERROR latitude: ", b, "longtitude: ", l
            input ("Enter a number ")
       
 






