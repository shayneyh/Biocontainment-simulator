#this program implements the Gillespie algorithm for a coupled logistic system
#ideally it would have input of the parameters R, a1, a2, K1, K2
# and of the number of steps per run and number of runs
#ideally it would have output of the mean time to extinction, or all extinction times
#created 3 September, 2013

#from pylab import *
import random
from math import log
#from scipy.integrate import odeint
import sys

#NTS make these input-able upon call of program
constants=[0]*6;
inputs=sys.argv;
constants[0]=1;constants[1]=1; #growth rate, r1,r2
constants[2]=8; constants[3]=constants[2]; #carrying capacity k1, k2
constants[4]=0.2; constants[5]=constants[4]; #a1,a2
#constants[2]=int(inputs[1]); constants[3]=constants[2];
#constants[4]=int(inputs[2]); constants[5]=constants[4];

t0 = 0;
nsteps = 10**5;
nruns = 1000;
nK = 1+1; #does powers of 2 up to 2**(nK-1)

def initialconditions(constants):
  rconstN=constants[0];rconstM=constants[1];
  KconstN=constants[2];KconstM=constants[3];
  alpha=constants[4];beta=constants[5];
  if(1-alpha*beta==0):
      npreyN0 = int(KconstN/2);
      npreyM0 = int(KconstM/2);
  else:
      npreyN0 = int((KconstN-alpha*KconstM)/(1-alpha*beta));
      npreyM0 = int((KconstM-beta*KconstN)/(1-alpha*beta));
  return [npreyN0,npreyM0]


#one step of the Gillespie algorithm
def stepVerhulst2(time, nprey1, nprey2, constants):
  n1 = random.random(); n2 = random.random();
  rconstN=constants[0];rconstM=constants[1];
  KconstN=constants[2];KconstM=constants[3];
  alpha=constants[4];beta=constants[5];
  #(*defining instantaneous rates*)
  ratePrey1Growth = rconstN*nprey1; 
  ratePrey1Death = nprey1*rconstN*(nprey1 + alpha*nprey2)/KconstN; 
  ratePrey2Growth = rconstM*nprey2; 
  ratePrey2Death = rconstM*nprey2*(beta*nprey1 + nprey2)/KconstM;
  rateTotal = ratePrey1Growth + ratePrey1Death + ratePrey2Growth + ratePrey2Death;
  r1 = ratePrey1Growth/rateTotal; r2 = ratePrey1Death/rateTotal; 
  r3 = ratePrey2Growth/rateTotal; r4 = ratePrey2Death/rateTotal;
  outPrey = [nprey1, nprey2]; outTime = time;
  #(*update population, then update time*)
  if(n1 < r1):
      outPrey[0] += 1
  elif(n1 < r1 + r2):
      outPrey[0] -= 1
  elif(n1 < r1 + r2 + r3):
      outPrey[1] +=1
  else:
     outPrey[1] -= 1
  outTime = time - log(n2)/rateTotal; #normalized against rate total since rate is an indication of how fast an event happens
  #(*output*)
  return [outTime, outPrey[0], outPrey[1]]


#one run of the Gillespie, to extinction
def runVerhulst2(t0,npreyN0,npreyM0):
    outV = [t0,npreyN0,npreyM0];
    counter=0;
    while((counter<nsteps)&(outV[1]>0)&(outV[2]>0)):
        outV = stepVerhulst2(outV[0],outV[1],outV[2],constants);  
        counter+=1
    return outV


#many runs of Gillespie, presumably at one value of parameters
def manyruns(nruns, t0,npreyN0,npreyM0):
    extinction = [0,0];
    etime = [0,0];
    lived=0;
    for j in range(nruns):
        tempoutV = runVerhulst2(t0, npreyN0, npreyM0);
        if(tempoutV[1]==0): #first population dies, then add 1 to extinction event, add up all the extinction time for all extinction
            extinction[0]+=1; etime[0]+=tempoutV[0];
        elif(tempoutV[2]==0):# same as previous one except for another species
            extinction[1]+=1; etime[1]+=tempoutV[0];
        else:#no population died within specified runs (steps)
            lived+=1;
#    print extinction;#
    for i in range(len(etime)): # i is number of species
        if(extinction[i]!=0):
            etime[i]=etime[i]/extinction[i]; #average extinction time only for runs that encoutered extinction
            extinction[i]=float(extinction[i])/float(nruns); # % of times one species becomes extinct
##        else:
##            print("none of species %s died",%(i+1))
    print lived;#
    return [etime, extinction]


#convert integer to parameter
def intconvert(k):
    return 2**k

#many parameter values
def main():
    output1=open("tvK.txt","w");output2=open("evK.txt","w");
    output1.write("#extinction time versus K in the symmetric case, case 3, with a1=a2=0.1\n");
    output2.write("#extinction chance versus K in the symmetric case, case 3, with a1=a2=0.1\n");
    etimevK=[];
    ehappvK=[];
    for k in range(1,nK): 
        
        param=intconvert(k);
        constants[2]=param;constants[3]=param; #vary carrying capacity in each iteration by a power of 2
        [npreyN0,npreyM0]=initialconditions(constants);
        temper=manyruns(nruns, t0,npreyN0,npreyM0); #outputs average extinction time and extinction probability given a constant k1, k2
        if(sum(temper[1])!=0): #if either species has a non-zero extinction probability
            etimevK.append([param,temper[0]]);
            ehappvK.append([param,temper[1]]);
            output1.write(str([param,temper[0]])+"\n");
            output2.write(str([param,temper[1]])+"\n");
        else:
            output1.write("none died\n");
            output2.write("none died\n");
        print("done that one")
    print(etimevK)
    print(ehappvK)
    output1.close();output2.close()


if __name__ == "__main__":
    main()

