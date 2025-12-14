% Code to generate data for Figs. 4, 5 and SM2 in Gabaldon et al, 
% "Heuristic Inferencence of the brain dynamical state from the r-spectrum of the correlation matrix"
% Code edited from  "Haimovici et al, Brain Organization into Resting State Networks Emerges at
% Criticality on a Model of the Human Connectome PRL 110, 178101 (2013)"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Header Information of the original code
%-------------------------------------------------------------------------%
% Simple Matlab code to implement the dynamics of the Greenberg-Hastings model
%
% J. M. Greenberg and S. P. Hastings, SIAM J. Appl. Math. 34, 515 (1978)
%
% running over the structure of the Human Brain Connectome described in:
%
% P. Hagmann, L. Cammoun, X. Gigandet, R. Meuli, C. J.
% Honey, V. J. Wedeen, and O. Sporns, PLoS Biol. 6, e159 (2008).
%
% The dynamics of the model as a function of the control parameter T was explored 
% by Haimovici et al, Brain Organization into Resting State Networks Emerges at
% Criticality on a Model of the Human Connectome PRL 110, 178101 (2013)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc
rvals=[0:0.001:1];  %r values used to perform the r-spectrum

%% Parameters for Fig. 4
%Tvec=[0.1:0.005:0.22]; % Vector of values of T: from 0.1 to 0.21
%total=1000000; % Number of iterations (time steps)
%Hemisphere=0 % Choose hemishpere 0 or 1. Full Brain: Hemisphere=2
%% Parameters for Fig. 5 and SM2
Tvec=[0.165:0.0005:0.21] % Vector of values of T. 
%Tvec=[0.165:0.0001:0.21] % Also use this with 5x resolution to estimate
%error bars for hemisphere 0
total=10000; % Number of iterations (time steps)
Hemisphere=0 %Hemisphere=1 in SM2
%% Connectome 

addpath('../functions')
fileaddr='../../HCP_MRI'
mkdir('Results')
DTIall=zeros(360);
files = dir([fileaddr '/*.mat']);

%%Compute the representative connectome
for R = 1:numel(files)
    fname = files(R).name;
    filepath = fullfile(fileaddr, fname);   
    load(filepath,'DTI');
    DTIall=DTIall+DTI;
end
DTIall=DTIall/R;

%Choose the hemisphere
if Hemisphere<2
    M=DTIall(1+180*Hemisphere:180+180*Hemisphere,1+180*Hemisphere:180+180*Hemisphere);  % the links to be used as the structure of the model
else
    M=DTIall;
end

%Normalize 
Norm=sum(M);
for i=1:size(M,2)
    for j=1:size(M,2)
        Mnew(i,j)=M(i,j)/Norm(j);
    end
end

M=Mnew';
%% SIMULATION 
%PARAMETERS

p1=0.001;   % Spontaneous activation probability (it set a very small rate of background activity)
p2=0.2;     % Refractory to Quiescent probability (it determine the lengh of time in which is active)
    
tindex=0;
for T=Tvec;     % The control parameter, is the activation threshold ( change for SUB-Critical, SUPER-Critical, or Critical dynamics)
    tindex=tindex+1;
    N=size(M,1);    % Number of nodes= 360
    S=zeros(total,N);  % Each Column in S is a node. Each row is a time step. States are labeled as 0: quiescent, 1: excited, 2: refractory.
    S(1,randi(N,1,10))=1;   % Initial excitation conditions at time 0, not really needed but here starts with 10 random active nodes
T
    for it=2:total
        % Single simulation step
        S(it,:)=S(it-1,:);
        S(it,S(it-1,:)==1)=2; % Previously active states go to Refractory period
        R=S(it-1,:)==2; 
        S(it,R&(rand(1,N)<p2))=0; % Nodes in the refractory period go back to quiescent state with probability r2
        Q=S(it-1,:)==0;
        S(it,Q&(rand(1,N)<p1))=1; % Spontaneous activation of quiescent nodes with probability r1.
        E=S(it-1,:)==1; 
        M1=M;
        M1(E,:)=0;
        M1(R,:)=0;
        M1(:,Q)=0;
        M1(:,R)=0; % Keep non zero values only at rows with quiescent nodes and columns with excited nodes (at t-1)
        W=sum(M1,2); % Sum contributions from excited neighbors(to determine if is going to fire)
        S(it,W> T)=1; % If sum is larger than threshold T, then fire.
    
    end % simulation for a given T ends

    % Delete the first 10% steps (let the system reach stationary state)
    S=S(fix(end/10):end,:);
    S=int8(S>0); 

    Act=mean(S');
    %Z=zscore(double(S));
    
    Fa(tindex)=mean(Act'); % Fraction of active nodes (the order parameter of the model)
    ACx=autocorr(Act);
    AC1(tindex)=ACx(2); % First autocorrelation coefficient 
    CorrMat=corr(double(S));  % Pearson correlation matrix
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % From the CorrMat, compute the diferent r-spectrum and r_c values
      [rc,results,MMat,NCL,AllSizesCell] = rc_function(CorrMat,rvals,1);
    
    % Store the rc values, and pairwise correlation mean (MMat) for diferent T
    Allrc_S2(tindex)=rc.S2;
    Allrc_sigm_S1(tindex)=rc.sigm_S1;
    Allrc_LG(tindex)=rc.LG;
    AllMMat(tindex)=MMat;


    lamb1(tindex)= eigs(CorrMat,1);

        
end % Move to another network

%Save values, for results visualization run the routine 'PlotHemispheres.m' 
save(sprintf('Results/ResultsHemi%d_Tvalues%d.mat',Hemisphere,numel(Tvec)),'Fa','AC1','lamb1','Allrc_S2','Allrc_sigm_S1','Allrc_LG','AllMMat','rvals','Tvec');
