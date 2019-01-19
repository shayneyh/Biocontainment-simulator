% this program implements the Gillespie algorithm for a coupled logistic system
% ideally it would have input of the parameters R, a1, a2, K1, K2
% and of the number of steps per run and number of runs
% ideally it would have output of the mean time to extinction, or all extinction times
% created 3 September, 2013
%
% modified August 14, 2014

%coupled of notes about Matlab:
%1. All arrays starts at index 1 instead of 0
%2. Use () instead of [] to access individual elements in an array
%3. Use disp instead of print and use ' instead of "
%4. Must set directory on the workspace to the same directory as the file
% in order to run file
%5. To change the filename, you must also change the name of the function
%to be able to run the file
%6. uint64 can support up to 2^64-1 which is about 1.84e19

%plot graphs,

%default function call:
%coupled_logistics_poisson_simulation(2.5076557,3.5e10, 3.5e10, 8, 20,24,10)
%plasmid_size in kb, each step represents an hour

function output_data = coupled_logistics_poisson_simulation(rconstWT, KconstWT, KconstMut, plasmid_size, copy_number, nsteps, nruns)

rconstMut = rconstWT + log(1-(7.2*10e-5)*plasmid_size*copy_number); % slope is normalized by e^(Rmut-Rwt),Rwt is normalized to one


% %initial population
% KconstWT = uint64(KconstWT);
% KconstMut = uint64(KconstMut);

filename = 'testdata.xlsx';
output_data = {};
for a =1:10
    alpha = 1-0.01*a; %a2 to be 0.0 - 1.0
    beta = alpha;
    
    
    for p = 1:10
        efficiency = 1-0.002*p; %varying efficiency
        nWT0 = KconstMut*efficiency;
        nMut0 = KconstMut*(1-efficiency); %initial population n1, n2
        t0 = 0;
        
        %run simulation
        [etime_average, extinction_prob, lived, prey1, prey2] = manyruns(nruns, t0,nWT0,nMut0); %manyruns(nruns, t0,npreyN0,npreyM0):outputs average extinction time and extinction probability given a constant k1, k2
        %Summarize results
        [avg_prey1, avg_prey2, popratio, stdev_prey1, stdev_prey2] =popcal(nruns, prey1, prey2);
        
        %output do excel
        
        %Moved the next line to be inside the second for loop because
        %initial population depends on efficiency, which is varied in this
        %loop
        constants = {'r1', rconstWT, 'r2', rconstMut, 'k1', KconstWT, 'k2', KconstMut; 'a1', alpha, 'a2', beta, 'p1', nWT0, 'p2',nMut0};
        avgpop = {'avg_population1', avg_prey1,'avg_population2', avg_prey2, 'ratio', popratio,'',''};
        stdev = {'std_pop1', stdev_prey1, 'std_pop2', stdev_prey2,'','','',''};
        
        
        if isempty(output_data)
           
            output_data = [constants; avgpop; stdev];
        else
            
            %append data from each loop to existing data
            output_data = [output_data; constants; avgpop; stdev];
        end
        
    end
end

xlswrite(filename,output_data);





%functions implementation--------------------------------------------------
    function [etime_average, extinction_prob, lived, prey1, prey2] = manyruns(nruns, t0,nWT0,nMut0)
        extinction = zeros(1,2);
        extinction_prob = zeros(1,2);
        etime = zeros(1,2);
        etime_average = zeros(1,2);
        lived=0;
        time=zeros(1,nruns);
        prey1=zeros(1,nruns);
        prey2=zeros(1,nruns);
        
        for j  = 1:nruns
            [temp_nWTtemp, temp_nMut, temp_time] = runVerhulst2(nWT0, nMut0,t0);
            prey1(j) = temp_nWTtemp;
            prey2(j) = temp_nMut;
            time(j) = temp_time;
            if temp_nWTtemp == 0 %first population dies, then add 1 to extinction event, add up all the extinction time for all extinction
                extinction(1) = extinction(1) + 1;
                etime(1) = etime(1) + temp_time;
            elseif temp_nMut==0 %same as previous one except for another species
                extinction(2) = extinction(2) + 1;
                etime(2)= etime(2) + temp_time;
                
            else   %no population died within specified runs (steps)
                lived = lived + 1;
            end
        end
        
        
        
        for i = 1:length(extinction) % i is number of species
            if extinction(i)~= 0
                etime_average(i) = etime(i)/extinction(i); %average extinction time only for runs that encoutered extinction
                extinction_prob(i) = float(extinction(i))/float(nruns); % of times one species becomes extinct
            else
                disp('none of the species died');
            end
        end
    end

    function [nWT,nMut,time] = runVerhulst2(nWT0,nMut0,t0) % combined stepVerhulst into runVerhulst
        
        
        %%one step of the Poisson algorithm
        nWT = nWT0;
        nMut = nMut0;
        time = t0;
        %defining instananeous growth rates
        %I took out the /1000 and *1000 because the integer type in matlab can
        %support integer value up to 2^64
        while time<nsteps && nWT>0 && nMut>0
            
            rateWTGrowth = rconstWT*nWT;
            rateWTDeath = nWT*rconstWT*(nWT + alpha*nMut)/KconstWT;
            rateMutGrowth = rconstMut*nMut;
            rateMutDeath = rconstMut*nMut*(beta*nWT + nMut)/KconstMut;
            
            WTgrowth = poissrnd(rateWTGrowth);
            WTdeath = poissrnd(rateWTDeath);
            Mutgrowth = poissrnd(rateMutGrowth);
            Mutdeath = poissrnd(rateMutDeath);
            
            nWT = nWT + WTgrowth - WTdeath;
            nMut = nMut + Mutgrowth - Mutdeath;
            time = time + 1;
        end
        
    end


    function [avg_prey1, avg_prey2, popratio, stdev_prey1, stdev_prey2] = popcal(nruns, prey1, prey2)
        
        stdev_prey1 = std(prey1,1);
        stdev_prey2 = std(prey2,1);
        avg_prey1 = sum(prey1)/nruns;
        avg_prey2 = sum(prey2)/nruns;
        totalPop = avg_prey1+avg_prey2;
        popratio = avg_prey1/totalPop;
        
    end





end



