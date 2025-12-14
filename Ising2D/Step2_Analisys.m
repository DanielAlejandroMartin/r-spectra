% Code to generate data for Fig. SM3 in Gabaldon et al, 
% "Heuristic Inferencence of the brain dynamical state from the r-spectrum of the correlation matrix"

clear all
close all
clc
%% Parameters
total=50000; %Number of iterations (time steps)
rvals=[0:0.001:1]; %r values used to perform the r-spectrum
Tvec=unique(sort([1.0:0.2:2 2.05:0.05:2.5 2.6:0.2:3.4 ])); % Vector of values of T: from 1.2 to 3.4 . Zoomed about 2.25


L=64; % Side lenght of the LxL grid
addpath('../functions')


mkdir("Analysis")


for run=1:5
run
    for T=Tvec(2:end)
        T
tic


filename = sprintf('Simulation/CorrIsing%d_N%d_T%f.mat',run,L^2,T);
         load(filename)

         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Cmt=Cmt+Cmt';
        for i=1:L^2
           Cmt(i,i)=1; 
        end
        CorrMat=Cmt;
        % From the CorrMat, compute the diferent r-spectrum and r_c values
        [rc,results,MMat,NCL,AllSizesCell] = rc_function(CorrMat,rvals,1);
       
       
       
               rc_S2=rc.S2;
        rc_sigm_S1=rc.sigm_S1;
        rc_LG=rc.LG;
        rc_Skew=rc.Skew;
        
        S1=results.S1;
        S2=results.S2;
        Skew=results.Skew;


        
               

 % Save values, for results visualization run the routine 'PlotIsing.m' 
             file=sprintf('Analysis/ResultsR%d_L%d_T%f.mat',run,L,T);
        
  save(file,'mag','rc_S2','rc_sigm_S1','rc_LG','AllSizesCell','Skew','rc_Skew')
        toc
    end  % end Tvec 
end % Move to another network