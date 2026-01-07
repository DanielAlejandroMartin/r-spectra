% Code to generate data for Fig. SM4 and SM6 in Gabaldon et al, 
% "Data-driven inference of brain dynamical states from the r-spectrum of correlation matrices"

clear all
close all
clc
%% Parameters
total=50000; %Number of iterations (time steps)
%rvals=[0:0.001:1]; %r values used to perform the r-spectrum
Tvec=unique(sort([1.0:0.2:2 2.05:0.05:2.5 2.6:0.2:3.4 ])); % Vector of values of T: from 1.2 to 3.4 . Zoomed about 2.25
%Tvec=unique(sort([2.35 2.4:0.2:3.4 ])); % Vector of values of T: from 1.2 to 3.4 . Zoomed about 2.25

ds=5; % Every how many steps you take a snapshot for computing Pearson Correlation.
L=64; % Side lenght of the LxL grid
addpath('../functions')
%%%%%%%
% Periodic Boundary Conditions
up_nei=2:L;
up_nei(L)=1;
dwn_nei=0:L-1;
dwn_nei(1)=L;


mkdir("Simulation")
%% Simulation
for run=1:5
    run
    rng(run*1234) % Choose random number generator for reproducibility  
    for T=Tvec
        Tinv=1/T;
        T
        tic
        if T==Tvec(1) % % For lower T start with half spins up, for following T start from pre-equilibrated state
                S=ones(L);
                S(rand(L) > 0.5) = -1;
        end        
        for it=2:total % Iterations start (MonteCarlo steps)
            for i=1:L^2
                 x=randi(L);
                 y=randi(L);
                 neigh = S(up_nei(x),y) + S(dwn_nei(x),y) + S(x,up_nei(y)) + S(x,dwn_nei(y));
                 dE= 2*S(x,y) * neigh;
                 if  rand() < exp(-dE*Tinv) % Turn the spin if it is energy-efficient
                     S(x,y) = -S(x,y);
                 end                 
            end % one MonteCarlo step
            if mod(it,ds)==0
                S_td(fix(it/ds),:)=reshape(S, [1, L^2]); % Store an snapshot for computing the correlation matrix
            end
            
            mag(it)=sum(sum(S)); % Magnetization for all the iterations
        end % iterations
         
        toc
        time=size(S_td,1);
        S_td=S_td(fix(time/10)+1:time,:); % Delete the first 10% steps (let the system reach stationary state)
        mag=mag(fix(time/10)+1:end); 
        ACx=autocorr(mag);  
        AC1=ACx(2); % First autocorrelation coefficient 

        CorrMat=corr(S_td); % Pearson correlation matrix

        CorrMat(isnan(CorrMat))=0; % Avoid NANs due to  spins that don't flip over simulation time at low T and short times. 
     
        Cmt=triu(CorrMat);
         filename = sprintf('Simulation/CorrIsing%d_N%d_T%f.mat',run,L^2,T);
         save(filename,'Cmt','mag','S','-v7.3')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end  % end Tvec 
end % Move to another network
