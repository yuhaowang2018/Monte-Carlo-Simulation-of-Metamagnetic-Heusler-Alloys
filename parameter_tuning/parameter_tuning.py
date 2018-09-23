#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 09:38:23 2017

@author: yuhao
"""
import numpy as np
import scipy as sp


from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from sklearn.gaussian_process.kernels import (Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct)

from skopt import Optimizer                                             
from matplotlib import pyplot as pl
import subprocess
import sys

from pexpect import pxssh
import getpass
import time

def generateinput(inputfilename,parameterlist):
    templatefile=open("template.inp","r")
    inputfile=open(inputfilename,"w")
    for line in templatefile:
        if("STRUCTPARA" in line):
            line="STRUCTPARA "+' '.join(map(str, parameterlist[0:6]))+"\n"
        elif("MagField" in line):
            line="MagField "+str(parameterlist[6])+"\n"
        elif("Kaniscub " in line):
            line="Kaniscub "+str(parameterlist[7])+"\n"
        elif("Kanistet " in line):
            line="Kanistet "+str(parameterlist[8])+"\n"
        inputfile.write(line)
    templatefile.close()
    inputfile.close()
    return

def readinputfile(filename):
    inputfile=open(filename,"r")
    for line in inputfile:
        if("STRUCTPARA" in line):
           para=line[:-1].split(" ")
           parameterlist=para[1:]
        elif("MagField" in line):
            para=line[:-1].split(" ")
            parameterlist.append(para[1])
        elif("Kaniscub " in line):
            para=line[:-1].split(" ")
            parameterlist.append(para[1])
        elif("Kanistet " in line):
            para=line[:-1].split(" ")
            parameterlist.append(para[1])
    return parameterlist

def readoutputfile(filename):
    outputfile=open(filename,"r")
    heat=[]
    cool=[]
    flag=0
    count=0
    for line in outputfile:
        if(flag==1):
          strains=line.split("\t")
          strain=[float(i) for i in strains]
          if(strain[0]>=200):
             heat.append(strain[1]) 
          if(strain[0]==400):
              flag=2
        elif (flag==2):
            strains=line.split("\t")
            strain=[float(i) for i in strains]
            if(strain[0]==410):
                count+=1
            if(count==2):
                flag=3
        elif (flag==3):
            strains=line.split("\t")
            strain=[float(i) for i in strains]
            cool.append(strain[1])
            if(strain[0]==200):
                flag=4
        elif ("Entropy" in line):
            flag=1
    cool.reverse()
    return [heat,cool]  

def checkstrain(X,y,ii):
    T=[]
    
    name=["NiCoMnIn","NiCoMnInH1","NiCoMnInH7"]
    epsheat=[0,0,0]
    epscool=[0,0,0]
    tt=[[260,240],[250,220],[0,0]]
    for i in range(len(name)):
        epsheat[i]=[]
        epscool[i]=[]
    for jj in range(len(name)):
        for i in range(17):
            T.append(i*5+200)
            t=i*5+200
            if (t<tt[jj][0]):
                epsheat[jj].append(1)
            else:
                epsheat[jj].append(0)
            if(t<tt[jj][1]):
                epscool[jj].append(1)
            else:
                epscool[jj].append(0)
   
    
    for i in range(ii):
        yy=0
        for ij in range(len(name)):
            inputfilename=name[ij]+".inp."+str(i+1)
            outputfilename=name[ij]+"."+str(i+1)+".out"
            parameters=readinputfile(inputfilename)
            parameterlist=[float(i) for i in parameters]
            if(ij==0): X.append([parameterlist[0],parameterlist[1],parameterlist[4],parameterlist[5]])
            cal=readoutputfile(outputfilename)
            heat=cal[0]

            cool=cal[1]

            yyy=0

            for j in range(len(heat)-1):
                yyy+=abs(heat[j]+heat[j+1]-epsheat[ij][j]-epsheat[ij][j+1])*5/2+abs(cool[j]+cool[j+1]-epscool[ij][j]-epscool[ij][j+1])*5/2

            yyy=yyy#*abs(deleps-delhc)
            yy+=yyy
        y.append(yy)
            
    return



def checkmag(X,y,ii):
    factor=35
    nameexp=["0.05T.txt","1T.txt","7T.txt"]
    epsheat=[[],[]]
    epsmeanheat=[[],[]]
    epscool=[[],[],[]]
    epsmeancool=[[],[],[]]
    pos=0
    for i in nameexp:
        expfile=open(i,"r")
        for line in expfile:
            spl=line[:-1].split(" ")
            if(i=="7T.txt"):
                data=[float(i) for i in spl[0:2]]
            else:
                 data=[float(i) for i in spl]
            if (len(data))==4:
                epscool[pos].append([data[0],data[1]])
                epsheat[pos].append([data[2],data[3]])
            else:
                epscool[pos].append([data[0],data[1]])
        expfile.close()
        epscool[pos].reverse()
        pos+=1
    temp=[10*i+200 for i in range(21)]
    epsheat=np.array(epsheat)
    epscool=np.array(epscool)
    for nn,i in enumerate(epscool):
        start=0;
        end=0;
        iters=0;
        for num,d in enumerate(i):
            tempdata=temp[iters]
            if (d[0]>=tempdata-0.6 and start==0):
                start=num
            if (d[0]>tempdata-0.1 and end==0):
                end=num
            if(start!=0 and end !=0):
                epsmeancool[nn].append(np.mean(i[start:end,1]))
                iters+=1
                start=0
                end=0
                if(abs(d[0]-temp[-1])<0.1): break
    
    for nn,i in enumerate(epsheat):
        start=0;
        end=0;
        iters=0;
        for num,d in enumerate(i):
            tempdata=temp[iters]
            if (d[0]>=tempdata-0.6 and start==0):
                start=num
            if (d[0]>tempdata-0.1 and end==0):
                end=num
            if(start!=0 and end !=0):
                epsmeanheat[nn].append(np.mean(i[start:end,1]))
                iters+=1
                start=0
                end=0
                if(abs(d[0]-temp[-1])<0.1): break
            
            
        
    name=["NiCoMnIn","NiCoMnInH1","NiCoMnInH7"]
    for i in range(ii):
        yy=0
        for ij in range(len(name)):
            inputfilename=name[ij]+".inp."+str(i+1)
            outputfilename=name[ij]+"."+str(i+1)+".out"
            # parameterlist: J, K, Uc, Ut, K1(U1),U2, MagField,Kaniscub,Kanistet
            parameters=readinputfile(inputfilename)
            parameterlist=[float(i) for i in parameters]
            if(ij==0): X.append(parameterlist[0:6]+parameterlist[7:9])
            cal=readoutputfile(outputfilename)
            heat=[i*factor for i in cal[0]]

            cool=[i*factor for i in cal[1]]
            yyy=0

            if (name[ij]!="NiCoMnInH7"):
                for j in range(len(heat)-1):
                    yyy+=abs(heat[j]+heat[j+1]-epsmeanheat[ij][j]-epsmeanheat[ij][j+1])*10/2+abs(cool[j]+cool[j+1]-epsmeancool[ij][j]-epsmeancool[ij][j+1])*10/2
                    
                yyy=yyy
                yy+=yyy
            else:
                for j in range(len(heat)-1):
                    yyy+=abs(cool[j]+cool[j+1]-epsmeancool[ij][j]-epsmeancool[ij][j+1])*10/2
                    
                yyy=yyy
                yy+=yyy
        y.append(yy/(2*len(heat)))
                
    
        

ii=1 #number of outputfiles at first
X=[]
y=[]
checkmag(X,y,ii)
       
np.random.seed(1)


kernel=C(1.0, (1e-5, 1e5))*RBF([1,1,1,1,1,1,1,1], (0.01, 9))+C(1.0, (1e-5, 1e5))  * ExpSineSquared(length_scale=1000, periodicity=50,
                                length_scale_bounds=(0.01, 10),
                                periodicity_bounds=(0.01, 1e5))
gp = GaussianProcessRegressor(kernel=kernel,n_restarts_optimizer=100,normalize_y=True)
bounds=[(4.5, 4.5),(2.35,2.35),(0,2),(0,2),(-2.5,-2.5),(-0.3,-0.3),(-0.012,-0.005),(-0.2,-0.1)]
# J,K,Uc,Ut,U1,U2,anisc,anist
opt = Optimizer(bounds, gp,n_random_starts=100,acq_func="EI",acq_optimizer="sampling",acq_optimizer_kwargs={"n_points":9999999999,"n_restarts_optimizer":100})
#opt.n_points=999999999
for i in range(ii):
    opt.tell(X[i],y[i])

maxEI=[]
for count in range(500):
    X=np.array(X)
    print('X combo: ')

    print(X)
   

    y=np.array(y)
    print ('ymin combo')
    print (X[np.argmin(y)])
    print (np.argmin(y))

    if(count>0):opt.tell(X[-1],y[-1])
    print(X[-1])
    next_x=opt.ask()
    print("next_x")
    print (next_x)

    X=np.append(X,next_x)
    X=np.reshape(X,(int(len(X)/8),8))

    pp=next_x
    parameterlist=[[pp[0], pp[1],  pp[2],  pp[3],  pp[4],pp[5], 0.05, pp[6],pp[7]],
                   [pp[0], pp[1],  pp[2],  pp[3],  pp[4],pp[5], 1, pp[6], pp[7]],
                   [pp[0], pp[1],  pp[2],  pp[3],  pp[4],pp[5], 7, pp[6], pp[7]]]
    inputfilename=["NiCoMnIn.inp."+str(ii+1),"NiCoMnInH1.inp."+str(ii+1),"NiCoMnInH7.inp."+str(ii+1)]
    
    for iii in range(len(inputfilename)):
        generateinput(inputfilename[iii],parameterlist[iii])

    commands = [
        "./yuhao_new_Magstr "+inputfilename[0],
        "./yuhao_new_Magstr "+inputfilename[1],
        "./yuhao_new_Magstr "+inputfilename[2],
        
    ]
    processes = [subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) for cmd in commands]
# do other things here..
# wait for completion
    for p in processes: p.wait()

    outputfilename=["NiCoMnIn."+str(ii+1)+".out","NiCoMnInH1."+str(ii+1)+".out","NiCoMnInH7."+str(ii+1)+".out"]
    while True:
        flag=0
        for n in outputfilename:
            while True:
               try:
                  checkoutput=open(n,'r')
                  break
               except FileNotFoundError:
                 print("Oops!  No file keep trying")
                 time.sleep(5)
            
            lines=checkoutput.readlines()
            #print("read "+n)
            if ('END\n' in lines):flag+=1
        if (flag==3):break
        else:  time.sleep(1000)

    ii+=1
    X=[]
    y=[]
    checkmag(X,y,ii)