% Code to generate data for Fig. SM3 in Gabaldon et al, 
% "Data-driven inference of brain dynamical states from the r-spectrum of correlation matrices"

% Code edited from  "Haimovici et al, Brain Organization into Resting State Networks Emerges at
% Criticality on a Model of the Human Connectome PRL 110, 178101 (2013)"


clear all
close all
%% Parameters
N_nodes=5000;
lambda=12.5; % parameter for the exponential distribuion of weights of the WS network
p1=0.001;   % Spontaneous activation probability (it set a very small rate of background activity) 
p2=0.2;     %  Refractory to Quiescent probability (it determine the lengh of time in which is active)
total=50000; % Number of iterations (time steps)
Tvec=unique(sort([0.10:0.01:0.22 0.15:0.002:0.18])); % Vector of values of T: from 0.1 to 0.22. Zoomed about 0.175
addpath('../functions')
rvals=[0:0.001:1]; %r values used to perform the r-spectrum
mkdir("Simulation")
%% 
%Tvec=Tvec(14:end)
for run=[1]
    run
    rng(run*1234)  % Choose random number generator for reproducibility  
    %% Create WS Network
    Net=WattsStrogatz(N_nodes,5,0.6);
    Nei_list=table2array(Net.Edges);

    Weights_list=rand(size(Nei_list,1),1);
    Weights_list=-log(Weights_list)/lambda; % Generate exponential distribuion of weights 

    Nei_list=[Nei_list; [Nei_list(:,2), Nei_list(:,1)]]; % Organize the neighbors (is important because of the dynamic rules)
    Weights_list=[Weights_list; Weights_list];
    
    %% Compute the Neighbor List and the weights. 
    [Nei_list,order]=sortrows(Nei_list, [1 2]); % Sort neighbors and save original order indices
    Weights_list=Weights_list(order);

    %Find where the list of neighbors of i begins
    for i=1: N_nodes
       begin(i)=min(find(Nei_list(:,1)==i));
    end
    begin(N_nodes+1)=size(Weights_list,1)+1 ;
    Nei_list=Nei_list(:,2); 

    clear Net
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Dynamics
    for T=Tvec
        T
%         if exist(sprintf('GH_runsDS/ResultsWS%d_N%d_T%f.mat',run,N_nodes,T))==0
%             T
%             run
        clear CorrMat
        tic
        S=zeros(total,N_nodes);% Store the states for the Nnodes across time_steps (0=quiescent, 1=active, 2=refractory)
        S(1,randi(N_nodes,1,10))=1;
        for it=2:total

            S(it,:)=S(it-1,:);
            % Previously active states go to Refractory period
            S(it,S(it-1,:)==1)=2; 
            % Nodes in the refractory period go back to quiescent state with probability r2
            R=S(it-1,:)==2; 
            S(it,R&(rand(1,N_nodes)<p2))=0; 

            % Spontaneous activation of quiescent nodes with probability r1.
            Q=S(it-1,:)==0;
            S(it,Q&(rand(1,N_nodes)<p1))=1; 

            % Activation due to Neighbors
            for i=find(Q)
               allj=Nei_list(begin(i):begin(i+1)-1); % Vector of Neighbors of i
               Excited=S(it-1,allj)==1; % Active Neighbors at previous time
               Weight=Weights_list(begin(i):begin(i+1)-1); % Weight of Neighbors
               Contrib=dot(double(Excited),Weight);

               if Contrib>T % Activate node i
               S(it,i)=1;
               end
            end

        end % End iterations

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Delete the first 10% steps (let the system reach stationary state)
        S=S(fix(end/10)+1:end,:); 

        S=int8(S>0); 
        sp=sparse(double(S)); % Convert to sparse matrix format for efficiency
        Act=mean(S');
        Fa=mean(Act'); % Fraction of active nodes (the order parameter of the model)
        ACx=autocorr(Act);
        AC1=ACx(2); % First autocorrelation coefficient 
        CorrMat=corr(double(S)); %Pearson correlation matrix
        toc
        
         Cmt=triu(CorrMat);
         filename = sprintf('Simulation/CorrWS%d_N%d_T%f.mat',run,N_nodes,T);
         S_end=S(end,:);
         save(filename,'Cmt','Act','S_end','-v7.3')
         %         end %existence
    end %Tvalues
end % Move to another network
