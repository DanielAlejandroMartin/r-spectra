% Code to generate  Figs. SM9 in Gabaldon et al,
% "Heuristic Inferencence of the brain dynamical state from the r-spectrum of the correlation matrix"


clear all
close all
clc


%% Convert MMP (multimodal parcellation) parcellations to RSN (resting state network) parcellations

load('glasser_to_aal_mapping.mat','RSNname','aal_number_assignments')
RSNnames=unique(RSNname);
for i=1:116
    for j=1:8
        if isequal(RSNname{i},RSNnames{j})
            RSNnumAAL(i)=j;
        end
    end
end
RSN=RSNnumAAL(cell2mat(aal_number_assignments));  % RSN contains the RSN associated for each MMP ROI, in total 8 RSNs

for i=1:8
    RSNsize(i)=sum(RSN==i);
end

clear RSNname aal_number_assignments RSNnumAAL
%%

rvals=[0:0.001:1];% r values used to perform the r-spectrum
N_rois=360; % Number of MMP ROIs


addpath('../functions');

%% Simulation

files = dir('../../HCP_MRI/*.mat'); % set the correct path

aS1=zeros(numel(files),numel(rvals),N_rois); % Contain 1 if a ROI particpate in the giant cluster GC, for each r value and for each R subject. If not, contain a 0
aS2=zeros(numel(files),numel(rvals),N_rois); % Contain 1 if a ROI particpate in the second largest cluster SCL, for each r value and for each R subject. If not, contain a 0

for R = 1:numel(files)
    R
    tic
    
    fname = files(R).name;
    filepath = fullfile('../../HCP_MRI', fname); % Path for the proccesed BOLD signals (and DTI connectivity matrix) of the 996 subjects
    load(filepath,'ROI_ts');
    
    TS=ROI_ts{1}'; % Select the session 1,2,3,4 (similar results are obtained for each session)
    
    M=corr(TS); %Pearson correlation matrix
    for i=1:N_rois
        M(i,i)=0;
    end
    
    %Compute the r-spectrum to obtain the rc, and AllS1(AllS2), that contains wich rois participate to the GC (SCL).
    [rc,results,MMat,NCL,AllSizesCell,AllS1,AllS2] = rc_function(M,rvals,1);
    
    % rc for the r-spectrum
    Allrc_S2(R)=rc.S2;
    
    
    
    ridx=0
    for r=rvals
        ridx=ridx+1;
        % Compute which ROIS particiapte at Giant Cluster
        for v=AllS1{ridx}
            aS1(R,ridx,v)=1;
        end
        % Compute which ROIS particiapte at Second largest Cluster
        for v=AllS2{ridx}
            aS2(R,ridx,v)=1;
        end
    end % rvals
    toc
end % End 996 subjects


save('RSN.mat', '-V7.3')

clear rnorm

% Compute the probabilities

% Normalize r-axis for each subject (rnorm = r / r_c from S2)

for subj=1:numel(files)
    rnorm(subj,:)=rvals/Allrc_S2(subj);
end


% Window binning over normalized r-values

r_widx = 0;
dr=0.001;    % Width of the window
rwvals=0:dr:1.5;               % Common normalized r-axis
epsilon=0.00000001; % To avoid problems when r_window=0
Prob1=zeros(numel(rwvals),360);
Prob2=zeros(numel(rwvals),360);

for r_window=rwvals
    
    r_widx=r_widx + 1;
    
    % Find all (subject, r-index) pairs falling inside the current window
    [a,b] = find((rnorm < r_window + dr) .* (rnorm > r_window - epsilon));
    
    % Average S1 and S2 contributions across all matched subjects/points
    for i = 1:size(a,1)
        Prob1(r_widx,:)= Prob1(r_widx,:)+squeeze(aS1(a(i),b(i),:))'/size(a,1);
        Prob2(r_widx,:)= Prob2(r_widx,:)+squeeze(aS2(a(i),b(i),:))'/size(a,1);
    end
    
end % r_window loop


% Compute RSN-averaged probabilities across ROIs

RSNProb1 = zeros(8, numel(rwvals));
RSNProb2 = zeros(8, numel(rwvals));

for i = 1:8
    RSNProb1(i,:) = mean(Prob1(:, RSN == i), 2);
    RSNProb2(i,:) = mean(Prob2(:, RSN == i), 2);
end


% Plot RSN probabilities for GC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%dS1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure SM9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1900) %Panel A
grayColor = [.7 .7 .7];
%plot(rwvals,Prob1,'linewidth',0.25,'Color',grayColor) %Optional
%hold on
plot(rwvals, RSNProb1([1 4 5 6 7 8], :), 'linewidth', 4)
legend(RSNnames{[1 4 5 6 7 8]})
xlabel('r/rc')
ylabel('Prob(GC)')

%
