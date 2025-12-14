function [r_c_estimates,results,MMat,NCL,AllSizesCell,AllS1,AllS2] = rc_function(M,rvals,ComputeDistances)
%RC_FUNCTION    Compute the r-spectrum and rc values from the correlation matrix

%   This script analyses correlation networks by constructing graphs from a
%   correlation matrix. For each threshold value r, an edge is included
%   between two nodes if their correlation exceeds r. By sweeping across a
%   range of thresholds, we build a sequence of graphs and compute several
%   quantities of interest associated with percolation transitions — the
%   so-called r-spectrum.
%   In particular, we compute the r-spectra of four variables: the size of
%   the largest cluster, the size of the second-largest cluster, the
%   variability of the largest cluster, and the characteristic path length.
%   Each r-spectrum has an associated critical value r_c, defined as the
%   threshold at which the corresponding r-spectrum reaches its maximum.

%   Inputs:
%       M,
%           correlation matrix of the system (size: N_nodes x N_nodes)
%       rvals,
%           vector containing the threshold values to perform the r-spectra
%       ComputeDistances,
%           set to 1 to include the computation of the r-spectrum for the
%           characteristic path length (this step becomes slow for large
%           matrices). Set to 0 to skip it.

%   Outputs:
%
%   In r_c_estimates:
%       S2,    threshold r at which the size of the second largest cluster reaches its maximum value
%       sigm_S1,   threshold r at which the variance of the size of the largest cluster reaches its maximum value
%       LG,    threshold r at which the characteristic path lenght of the graph reaches its maximum value
%       Skew,  threshold r at which the skewness of the cluster sizes reaches its maximum value
%       cvdeg,  threshold r at which the coefficient of variation of the degree distribution reaches its maximum value
%
%   In results:
%
%       S1,       r-spectrum for the size of the largest cluster (i.e., S1
%       as a function of r values)
%       S2,       r-spectrum for the size of the second largest cluster
%       sigm_S1,      r-spectrum for the size of the variance of the size of the largest cluster
%       LG,       r-spectrum for the size of the characteristic path lenght of the graph
%       Skew,     r-spectrum for the size of the skewness of the cluster sizes
%       meandeg,  r-spectrum for the mean degree
%       sddeg,    r-spectrum for the standard deviation of the degree
%       cvdeg,    r-spectrum for the coefficient of variation of the degree
%   Others:
%
%       MMat,     Mean pairwise correlation value
%       NCL,      number of clusters in the system for each threshold r
%       AllSizesCell, Cell of all Cluster sizes for each threshold r
%       AllS1,  Cell of the elements of the largest cluster for each threshold r
%       AllS2,  Cell of the elements of the second largest cluster for each threshold r



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


r_c_estimates = struct();
results = struct();


LG=0;
N_rois=size(M,2); % Total nodes in network

for i=1:N_rois; M(i,i)=0; end % Set the Pearson self-correlations to zero

MMat=sum(sum(M))/(N_rois*N_rois-N_rois); % Mean pairwise correlation value
rindex=0;

for r=rvals
    rindex=rindex+1;
    
    W=(M>r); % Tresholded matrix
    U=adj2cluster(W);  % Cell array where each element contains node indices for a cluster
    
    % Sort the cells of clusters according to their size
    
    AllSizes=[];
    
    for i=1:size(U,2)
        AllSizes(i)=size(U{i},2);
    end
    
    [AllSizes,order_sizes]=sort(AllSizes,'descend');
    Isolated(rindex)=sum(AllSizes==1);
    NCL(rindex)=size(U,2); %-Isolated(rindex);
    % Extract the sizes of the largest and second largest clusters
    S1(rindex)=AllSizes(1);
    AllS1{rindex}=U{order_sizes(1)};
    AllS2{rindex}={};
    S2(rindex)=0;
    if size(AllSizes,2)>1
        S2(rindex)=AllSizes(2);
        AllS2{rindex}=U{order_sizes(2)};
    end
    
    Skew(rindex)=skewness(AllSizes(AllSizes>1));
    
    AllSizesCell{rindex}=AllSizes(AllSizes>1);
    
    
    
    %    [AllSizes,order_sizes]=sort(AllSizes,'descend');
    
    %% Degree
    deg=sum(W);
    meandeg(rindex)=mean(deg)+1; %avoid degree being zero
    sddeg(rindex)=std(deg);
    
    
    
    
    
    
    if ComputeDistances==1
        [LG(rindex),Deg,Nedges] = characteristic_path_length(W);        % Compute distance-associated magnitudes if desired (becomes slow for big graphs)
        
        
        
    end  %Compute Time-Consuming Routines
    
    
    
end %rvals

%For each r, compute the variance of S1 over a  10-point window about r
sigm_S1=zeros(size(rvals));
for rindex=1:size(rvals,2)
    mini=max(rindex-10,1);
    maxi=min(rindex+10,size(rvals,2));
    didx=randsample(mini:maxi,10);
    sigm_S1(rindex)=var(S1(didx)/N_rois);
end



%% Compute the rc values from the r-spectrum
cvdeg=sddeg./meandeg;


results.S1=S1;
results.S2=S2;
results.Skew=Skew;
results.meandeg=meandeg;
results.sddeg=sddeg;
results.cvdeg=cvdeg;
results.sigm_S1=sigm_S1;

%% Find the maximums of each meassure
[~,idx]=max(S2);
r_c_estimates.S2=rvals(idx);


[~,idx]=max(Skew);
r_c_estimates.Skew=rvals(idx);

[~,idx]=max(cvdeg);
r_c_estimates.cvdeg=rvals(idx);



[~,idx]=max(sigm_S1);
r_c_estimates.sigm_S1=rvals(idx);


if ComputeDistances>0
    
    results.LG=LG;
    
    
    [~,idx]=max(LG);
    r_c_estimates.LG=rvals(idx);
    
    
end


end



