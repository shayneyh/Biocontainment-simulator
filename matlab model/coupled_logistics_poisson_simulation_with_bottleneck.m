%% Poisson Process Model for Coupled Logistic System
%% Introduction
% This program models the population dynamics of two species using poisson 
% stochastic modelling. The growth profile are 
% modelled using coupled logistic system
% Governing equations for growth and death rate:
% dx/dt = r1x(1-(x+ay)/k1)
% dy/dt = r2y(1-(y+bx)/k2)
% The two species modelled in this program represenst wild type (WT) 
% MG1655 bacteria and MG1655 mutant(contains plasmids).
% r1,k1,k2 are found by experiment while r2 depends on several factors
% and is implemented at the begninning of the program
%
% The program varies a,b in the coupled logistics eqautions as well as
% CRSIPR efficiency, which will change the initial WT and mutant
% population.
%
% The program run the poisson stochstics model several times for each
% condition to obtain statistically significant data. For each condition,
% it calculates 5 output parameters: WT and mutant population at the end 
% of smulation, WT population/total population(ratio) and standard 
% deviation of WT and mutant population. The result is summarized in a
% table and exported to .xlsx file. The program can also generate 3-D 
% graphs for any input and output variables
%
% The input parameters are as follows:
% rconstWT (r1): growth rate constant for WT MG1655
% KconstWT (k1): carrying capacity for WT MG1655 per cm^2
% KconstMut (k2): carrying capacity for mutant
% plasmid_size: size of plasmid transformed into mutant (unit in kb)
% copy_number: nubmer of plasmid in each mutant
% nsteps: each step of the poisson simulation represents 10 mins, nsteps 
%         specifies how many step the model will run before the first
%         bottleneck occurs
% nruns: number of trials for each condition
%
% The output variable of this program, output_data, is a cell matrix
% containing values of input and output parameters for each condition with
% a header in the first row. The columns are as follows:
% 1. r1, 2. r2, 3. k1, 4. k2, 5. a1, 6. a2, 7. p1, 8. p2, 9. avg_prey1
% 10. avg_prey2, 11. popratio, 12. stdev_prey1, 13. stdev_prey2
%
% Published on August 15, 2014
% Credit to Matthew Badali, Department of Physics, University of Toronto
%
% default function call:
% coupled_logistics_poisson_simulation_with_bottleneck(2.5076557,4.46e8,4.46e8,8,20,24,10) Dont use this
% coupled_logistics_poisson_simulation_with_bottleneck(1.5188,4.46e8,4.46e8,8,20,5,10)
%
%% Implementation
% *main function*

function output_data = coupled_logistics_poisson_simulation_with_bottleneck(rconstWT,...
    KconstWT, KconstMut, plasmid_size, copy_number, nsteps, nruns)

%modify r2 based on r1 and total plasmid DNA
%slope is normalized by e^(Rmut-Rwt),Rwt is normalized to one
rconstMut = rconstWT + log(1-(7.2*10e-5)*plasmid_size*copy_number); 

%header of the output file
output_data = {'r1','r2','k1','k2','a1','a2','p1','p2','avg_population1'...
    ,'avg_population2','ratio', 'std_pop1', 'std_pop2'};

%indicate range of conditions before running the simulation
% interval_a = 0.05; %decrement by 0.01
% upper_bound_a = 1.5;
% lower_bound_a = 1/1.5; %a ranges from 1~0.9
%total_points_a = round((1-lower_bound_a)/interval_a)+1;

% interval_b = 0.01;
% upper_bound_b = 0.99;
% lower_bound_b = 0.8;
% total_points_b = round((1-lower_bound_b)/interval_b)+1;

%vary a1, a2
alpha_range = [1/1.5 1/1.45 1/1.4 1/1.35 1/1.3 1/1.25 1/1.2 1/1.15 1/1.1 1/1.05 1 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5];
efficiency_range = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99];
for a = alpha_range %assume efficiency = 50 and dilution factor = 0.01
    alpha = a;
    beta = 1/alpha;
    efficiency = 0.9;
    nWT0 = KconstMut*efficiency;
    nMut0 = KconstMut*(1-efficiency); %initial population n1, n2
    t0 = 0;
    dilution_factor = 0.01;
    %run simulation
    [etime_average, extinction_prob, lived, prey1, prey2] = ...
        manyruns(nruns, t0,nWT0,nMut0, dilution_factor);
    
    %Summarize results
    [avg_prey1, avg_prey2, popratio, stdev_prey1, stdev_prey2]...
        = popcal(nruns, prey1, prey2);
    
    %Moved the next line to be inside the second for loop because
    %initial population depends on efficiency, which is varied in this
    %loop
    %Add rows to existing data after each loop
    output_data = [output_data; {rconstWT, rconstMut, KconstWT,...
        KconstMut, alpha, beta, nWT0, nMut0, avg_prey1, avg_prey2,...
        popratio, stdev_prey1, stdev_prey2}];
    
end

%vary CRISPR efficiency, assuming alpha = beta = 1 and dilution factor =
%0.01
for b = efficiency_range
    alpha = 1;
    beta = 1;
    efficiency = b;
    nWT0 = KconstMut*efficiency;
    nMut0 = KconstMut*(1-efficiency); %initial population n1, n2
    t0 = 0;
    dilution_factor = 0.01;
    
    %run simulation
    [etime_average, extinction_prob, lived, prey1, prey2] = ...
        manyruns(nruns, t0,nWT0,nMut0, dilution_factor);
    
    %Summarize results
    [avg_prey1, avg_prey2, popratio, stdev_prey1, stdev_prey2]...
        = popcal(nruns, prey1, prey2);
    
    %Moved the next line to be inside the second for loop because
    %initial population depends on efficiency, which is varied in this
    %loop
    %Add rows to existing data after each loop
    output_data = [output_data; {rconstWT, rconstMut, KconstWT,...
        KconstMut, alpha, beta, nWT0, nMut0, avg_prey1, avg_prey2,...
        popratio, stdev_prey1, stdev_prey2}];
end


%vary dilution factors
for c = [0.1 0.01 0.001]
    [etime_average, extinction_prob, lived, prey1, prey2] = ...
        manyruns(nruns, t0,nWT0,nMut0, dilution_factor);
    %Summarize results
    [avg_prey1, avg_prey2, popratio, stdev_prey1, stdev_prey2]...
        = popcal(nruns, prey1, prey2);
    
    %Moved the next line to be inside the second for loop because
    %initial population depends on efficiency, which is varied in this
    %loop
    %Add rows to existing data after each loop
    output_data = [output_data; {rconstWT, rconstMut, KconstWT,...
        KconstMut, alpha, beta, nWT0, nMut0, avg_prey1, avg_prey2,...
        popratio, stdev_prey1, stdev_prey2}];
    
end
%Output do excel
filename = 'testdata2.xlsx';
xlswrite(filename,output_data);


%Plot graphs
%1. r1, 2. r2, 3. K1, 4. K2, 5. a1, 6. a2, 7. p1, 8. p2, 9. avg_prey1
%10. avg_prey2, 11. popratio, 12. stdev_prey1, 13. stdev_prey2

%function definition: plotgraph(x,y,z,x_label, y_label, z_label)
plotgraph(5,7,9, 'a1', 'p1', 'WT population')
title('Graph of average WT population vs a1 and initial WT population')
%Graph of Average Mutant Population vs a1 and Initial WT Population
plotgraph(5,7,10, 'a1', 'p1', 'Mutant population')
title('Graph of average mutant population vs a1 and initial WT population')
%Graph of WT Population/Total Population vs a1 and Initial WT Population
plotgraph(5,7,11,'a1', 'p1', 'ratio')
title...
('Graph of WT Population/Total Population vs a1 and Initial WT Population')

%%    
% *Sub functions*

%runs simulation for "nruns" times given initial populations 
function [etime_average, extinction_prob, lived, prey1, prey2] =...
        manyruns(nruns, t0,nWT0,nMut0, dilution_factor)
    extinction = zeros(1,2);
    extinction_prob = zeros(1,2);
    etime = zeros(1,2);
    etime_average = zeros(1,2);
    lived=0;
    time=zeros(1,nruns);
    prey1=zeros(1,nruns);
    prey2=zeros(1,nruns);
    temp_nWT = nWT0;
    temp_nMut = nMut0;
    temp_time = t0;

    for i  = 1:nruns
        
        [temp_nWT,temp_nMut,temp_time] = runVerhulst2(temp_nWT,temp_nMut,temp_time);
%         plot(temp_nWT)
%         plot(temp_nMut)
            
        cycles = 0; %define 100 cycles as survival threshold
       while temp_nWT(end) > 0 &&temp_nMut(end)> 0 && cycles <= 100
        %bottleneck
        [temp_nWT,temp_nMut,temp_time] = bottle_neck(temp_nWT,temp_nMut,temp_time, dilution_factor);
        
        cycles = cycles + 1;
       end
       
        prey1(i) = temp_nWT;
        prey2(i) = temp_nMut;
        time(i) = temp_time;

        %if first population dies, then add 1 to extinction event,
        %then add up all the extinction time for all extinction
        if temp_nWT == 0 
            extinction(1) = extinction(1) + 1;
            etime(1) = etime(1) + temp_time;
        %same as previous one except for another species
        elseif temp_nMut==0 
            extinction(2) = extinction(2) + 1;
            etime(2)= etime(2) + temp_time;
        %no population died within specified runs (steps)
        else   
            lived = lived + 1;
        end
    end           

    for j = 1:length(extinction) % i is number of species
        if extinction(j)~= 0
            %average extinction time only for runs that 
            %encoutered extinction
            etime_average(j) = etime(j)/extinction(j);
            % # of times one species becomes extinct
            extinction_prob(j) = extinction(j)/nruns; 
        end
    end
end

%runs one step of the Poisson algorithm
%simulate population growth of both species at carrying capacity
function [nWT,nMut,time] = runVerhulst2(nWT0,nMut0,t0)

    nWT = nWT0;
    nMut = nMut0;
    time = t0;
    %defining instananeous growth rates       
    while time(end)-t0(end)<nsteps && nWT(end)>0 && nMut(end)>0

        rateWTGrowth = rconstWT*nWT(end);
        rateWTDeath = nWT(end)*rconstWT*(nWT(end) + alpha*nMut(end))/KconstWT;
        rateMutGrowth = rconstMut*nMut(end);
        rateMutDeath = rconstMut*nMut(end)*(beta*nWT(end) + nMut(end))/KconstMut;

        WTgrowth = poissrnd(rateWTGrowth);
        WTdeath = poissrnd(rateWTDeath);
        Mutgrowth = poissrnd(rateMutGrowth);
        Mutdeath = poissrnd(rateMutDeath);

        nWT = [nWT nWT(end) + WTgrowth - WTdeath];
        nMut = [nMut nMut(end) + Mutgrowth - Mutdeath];
        time = [time time(end)+1/6];
    end    

end

%simulate effect of bottleneck at 1% (GenDrift.pdf,Evaluating the Impact of Population Bottlenecks
%in Experimental Evolution http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1462272/pdf/12399403.pdf), ends when total population reaches carrying
%capacity (assume growth rate stays the same after bottleneck, http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2902570/#!po=22.7273 fig.2)
function [nWT,nMut,time] = bottle_neck(nWT0,nMut0,t0, dilution_factor)
    
    pop_after_bottleneck = round((nWT0(end)+nMut0(end))*dilution_factor);
    WT_fraction = nWT0(end)/(nWT0(end)+nMut0(end));
    
    WT_after_bottleneck = sum(rand(1,round(pop_after_bottleneck))<WT_fraction);
    nWT = [nWT0 WT_after_bottleneck];
    nMut = [nMut0 pop_after_bottleneck-WT_after_bottleneck];
    time = [t0 t0(end)];
    while nMut(end)+nWT(end) < KconstWT
        rateWTGrowth = rconstWT*nWT(end);
        rateWTDeath = nWT(end)*rconstWT*(nWT(end) + alpha*nMut(end))/KconstWT;
        rateMutGrowth = rconstMut*nMut(end);
        rateMutDeath = rconstMut*nMut(end)*(beta*nWT(end) + nMut(end))/KconstMut;

        WTgrowth = poissrnd(rateWTGrowth);
        WTdeath = poissrnd(rateWTDeath);
        Mutgrowth = poissrnd(rateMutGrowth);
        Mutdeath = poissrnd(rateMutDeath);

        nWT = [nWT nWT(end) + WTgrowth - WTdeath];
        nMut = [nMut nMut(end) + Mutgrowth - Mutdeath];
        time = [time time(end) + 1/6]; 
        n = n+1;
    end
    
    
end
    
    

%summarize results for all trials for a particular condition
function [avg_prey1, avg_prey2, popratio, stdev_prey1, stdev_prey2] =...
        popcal(nruns, prey1, prey2)

    stdev_prey1 = std(prey1,1);
    stdev_prey2 = std(prey2,1);
    avg_prey1 = sum(prey1)/nruns;
    avg_prey2 = sum(prey2)/nruns;
    totalPop = avg_prey1+avg_prey2;
    popratio = avg_prey1/totalPop;

end

%plot 3-D graphs, x,y are indep variables (usually input variables) 
%and z is the dependent variable (usually output variables)
function plotgraph(x,y,z,x_label, y_label, z_label)
    %turn output_data from cell array to a matrix for plotting purposes
    data = cell2mat(output_data(2:size(output_data,1),:));
    %5 denotes a1, 7 denotes p1 and 11 denotes popratio. See
    %introduction for references
    x = reshape(data(:,x),length(data(:,x))/total_points_a,total_points_a);
    y = reshape(data(:,y),length(data(:,y))/total_points_a,total_points_a);
    z = reshape(data(:,z),length(data(:,z))/total_points_a,total_points_a);
    %eg, if the array is 333222111 then it is reshaped into 333
    %222
    %111 in order to create a meshgrid for 3-D plotting

    %plot graph
    figure
    mesh(x,y,z)
    xlabel(x_label);
    ylabel(y_label);
    zlabel(z_label);

end

end



