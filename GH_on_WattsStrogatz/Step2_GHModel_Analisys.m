% Code edited from  "Haimovici et al, Brain Organization into Resting State Networks Emerges at
% Criticality on a Model of the Human Connectome PRL 110, 178101 (2013)"


clear all
close all
%% Parameters
N_nodes=5000;
Tvec=unique(sort([0.10:0.01:0.22 0.15:0.002:0.18])); % Vector of values of T: from 0.1 to 0.22. Zoomed about 0.175
rvals=[0:0.001:1]; %r values used to perform the r-spectrum
addpath('../functions')
mkdir("Analysis")

%Tvec=Tvec(15:end)
for run=[1:40]

    for T=Tvec(1:end)
      
        tic
 file=sprintf('Analysis/ResultsR%d_L%d_T%f.mat',run,N_nodes,T);
if ~isfile(file)
        filename = sprintf('Simulation/CorrWS%d_N%d_T%f.mat',run,N_nodes,T)
        load(filename)
        ACx=autocorr(Act);  
        AC1=ACx(2); 
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Cmt=Cmt+Cmt';
        for i=1:N_nodes
           Cmt(i,i)=1; 
        end
        CorrMat=Cmt;
        
        % From the CorrMat, compute the diferent r-spectrum and r_c values
        [rc,results,MMat,NCL,AllSizesCell] = rc_function(CorrMat,rvals,1);
        
       
        % Save values, for results visualization run the routine 'PlotIsing.m' 
        
        file=sprintf('Analysis/ResultsR%d_L%d_T%f.mat',run,N_nodes,T);
        save(file,'rc','AC1','Act')
        toc

end %if file exists
    end  % end Tvec 
end % Move to another network