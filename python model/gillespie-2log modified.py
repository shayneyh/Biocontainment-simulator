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
constants=[0]*8;
inputs=sys.argv;
constants[0]=1.0;constants[1]=0.1; #growth rate, r1,r2
constants[2]= 100000.0; constants[3]=constants[2]; #carrying capacity k1, k2
constants[4]=1.0; constants[5]=constants[4]; #a1,a2
constants[6]=99000.0; constants[7]= 1000.0; #initial population n1, n2



t0 = 0;
nsteps = 10**3;
nruns = 1;

def initialconditions(constants):
  rconstWT=constants[0];rconstMut=constants[1];
  KconstWT=constants[2];KconstMut=constants[3];
  alpha=constants[4];beta=constants[5];
  nWT0 = constants[6];
  nMut0 = constants[7];
  
  return [nWT0,nMut0]


#one step of the Gillespie algorithm
def stepVerhulst2(time, nWT, nMut, constants):
  event = random.random(); timestep = random.random();
  rconstWT=constants[0];rconstMut=constants[1];
  KconstWT=constants[2];KconstMut=constants[3];
  alpha=constants[4];beta=constants[5];
  #(*defining instantaneous rates*)
  rateWTGrowth = rconstWT*nWT; #990
  rateWTDeath = nWT*rconstWT*(nWT + alpha*nMut)/KconstWT; #990
  rateMutGrowth = rconstMut*nMut; #10 
  rateMutDeath = rconstMut*nMut*(beta*nWT + nMut)/KconstMut;  #0.1
  rateTotal = rateWTGrowth + rateWTDeath + rateMutGrowth + rateMutDeath;
  r1 = rateWTGrowth/rateTotal; r2 = rateWTDeath/rateTotal; 
  r3 = rateMutGrowth/rateTotal; r4 = rateMutDeath/rateTotal;
  outPrey = [nWT, nMut]; outTime = time;
  #(*update population, then update time*)
  ##print event, r1, r2, r3, r4
  if(event < r1):
      outPrey[0] += 1
  elif(event < r1 + r2):
      outPrey[0] -= 1
  elif(event < r1 + r2 + r3):
      outPrey[1] +=1
  else:
     outPrey[1] -= 1
  outTime = time - log(timestep)/rateTotal; #normalized against rate total since rate is an indication of how fast an event happens
  ##print("time: %d", outTime);
  print outPrey[0], outPrey[1]
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
    for i in range(len(etime)): # i is number of species
        if(extinction[i]!=0):
            etime[i]=etime[i]/extinction[i]; #average extinction time only for runs that encoutered extinction
            extinction[i]=float(extinction[i])/float(nruns); # % of times one species becomes extinct
        else:
            print("none of species %s died",(i+1))
    return [etime, extinction, lived]




#many parameter values
def main():

   #change the program to determine if the mutant population stays under the init population instead of growing
    output1=open("summary2.txt","w");
    output1.write("Parameters:r1 = " + str(constants[0]) + " r2 = " + str(constants[1]) + " k1 = " + str(constants[2]) + " k2 = " + str(constants[3]) + " a1 = " + str(constants[4]) + " a2 = " + str(constants[5]) + " p1 = " + str(constants[6]) + " p2 = " + str(constants[7]) + "\n");
    output1.write("#average extinction time\n");
    output1.write("#extinction probability\n");
    output1.write("#lived\n");
  ##  etimevK=[];
  ##  ehappvK=[];
    
    
    [npreyN0,npreyM0]=initialconditions(constants);
    temper=manyruns(nruns, t0,npreyN0,npreyM0); #manyruns(nruns, t0,npreyN0,npreyM0):outputs average extinction time and extinction probability given a constant k1, k2
    if(sum(temper[1])!=0): #if either species has a non-zero extinction probability
   ##         etimevK.append([param,temper[0]]);
     ##       ehappvK.append([param,temper[1]]);
          output1.write(str(temper[0])+"\n");
          output1.write(str(temper[1])+"\n");
          output1.write(str(temper[2])+"\n");
    else:
            output1.write("none died\n");
            
    
    ##print(etimevK)
    ##print(ehappvK)
    output1.close();


if __name__ == "__main__":
    main()

