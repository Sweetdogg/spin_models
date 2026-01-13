import numpy as np
import pathlib
import random as rd
isboot=0
index=6
boot0sample=10
boot0bins=10
filepath = sorted(pathlib.Path('.').rglob('out.dat'))
data0list=[]
for i in filepath:
    data=np.loadtxt(i,usecols=(0,1,3,5))
    boot0samples=[]
    data0samples=[]
    if(not isboot) :
        bins=len(data)
        data0samples=data
        data0samples[:,2]=data[:,3] * data[:,0]**2
        # data0samples[:,2]=(data[:,2]+data[:,3]/data[:,5]**2)*(1-data[:,4]**2/data[:,5]**2)
        data0list.append([data[0][0],data[0][1],data0samples.mean(0)[2],np.sqrt(data0samples.var(0)[2]/(bins-1))])
    else:
        for j in range(boot0sample):
            boot0val=np.array(rd.choices(data,k=boot0bins))
            boot0samples.append(boot0val.mean(0))  #boot0samples.append(boot0val.mean(0)[7]/boot0val.mean(0)[6]**2)
        boot0samples=np.array(boot0samples)
        data0list.append([data[0][0],data[0][1],boot0samples.mean(0)[2],boot0samples.std(0)[2]])
np.savetxt('processed_data/xx/out00.txt',sorted(data0list,key=lambda x:x[0]),fmt='%d'+'   %.5f'+2*'   %.12f')
