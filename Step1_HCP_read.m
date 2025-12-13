% Code to generate data for Figs. 1, 2, 3 and SM1 in Gabaldon et al,
% "Heuristic Inferencence of the brain dynamical state from the r-spectrum of the correlation matrix"


clear all
close all

for ses= [1 5] % Number of Session 1-4. Sesion 5 is for all concatenated
    clearvars -except ses
    rvals=[0:0.001:1]; % r values used to perform the r-spectrum
    N_rois=360;
    fileaddr='../../HCP_MRI'; % Path for the proccesed BOLD signals (and DTI connectivity matrix) of the 996 subjects
    addpath('../functions');
    files = dir([fileaddr '/*.mat']);
    
    for R = 1:numel(files) %Subject
        fname = files(R).name;
        filepath = fullfile(fileaddr, fname);
        load(filepath);
        subjCode(R) = str2double(erase(fname, '.mat')); % Subject label
        R
        if ses<5
            TS=ROI_ts{ses}'; % Contains the proccesed BOLD signal. Each row in TS is a node. Each column is a time step.
        else
            TS=[];
            for i=1:4
                x=ROI_ts{i}';
                TS=[TS; x-mean(x)]; % Stack demeaned time series
            end
        end % if ends
        
        %% Compute the first autocorrelation coefficient
        Z=zscore(TS)';
        act=reshape(Z',[1 numel(Z)]);
        ac=autocorr(act);
        AC1(R)=ac(2);   % First autocorrelation coefficient
        
        
        CorrMat=corr(TS); % Pearson correlation matrix
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % From the CorrMat, compute the diferent r-spectrum and r_c values
        
        [rc,results,MMat,NCL,AllSizesCell] = rc_function(CorrMat,rvals,1);
        Allrc_S2(R)=rc.S2;
        Allrc_sigm_S1(R)=rc.sigm_S1;
        Allrc_LG(R)=rc.LG;
        Allrc_cvdeg(R)=rc.cvdeg;
        
        
        
        AllS1(R,:)=results.S1;
        AllS2(R,:)=results.S2;
        AllNCL(R,:)=NCL;
        
        
        if R<11 %Save for Fig1. and Fig2, 10 timeseries
            S1example(R,:)=results.S1;
            S2example(R,:)=results.S2;
            sigm_S1example(R,:)=results.sigm_S1;
            LGexample(R,:)=results.LG;
            cvdegexample(R,:)=results.cvdeg;
            degexample(R,:)=results.meandeg;
            sddegexample(R,:)=results.sddeg;
        end
        
        nbkt=false; % for pc estimation from the non backtracking matrix, set nbkt= true
        [pc, results] = pc_function(CorrMat, 'nonbacktracking',nbkt);
        
        
        Allpc_sp(R)=pc.spectral;
        Allpc_mom(R)=pc.degree_moments;
        
        
        
        
        
    end % Move to another subject R
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Extract subject age range
    datatable=readtable("HCP_YA_subjects_2025_10_15_11_40_46.csv");
    [tf, loc] = ismember(subjCode, datatable.Subject);
    Ages = datatable.Age(loc) ;
    %cell2mat(ages)==cell2mat(agedix(1))
    agedix=["22-25", "26-30", "31-35", "36+"];
    agedix(1)==Ages;
    
    % Relate the rc with age for each subject
    for ii=1:4
        subj=  agedix(ii)==Ages;
        Npat(ii) = sum(subj);
        Meanrc_S2(ii)=mean(Allrc_S2(subj));
        SDrc_S2(ii)=std(Allrc_S2(subj));
        Meanrc_sigm_S1(ii)=mean(Allrc_sigm_S1(subj));
        SDrc_sigm_S1(ii)=std(Allrc_sigm_S1(subj));
        Meanrc_LG(ii)=mean(Allrc_LG(subj));
        SDrc_LG(ii)=std(Allrc_LG(subj));
    end
    
    
    % Save values, for results visualization run the routine 'Step2_PlotHCP.m'
    
    save(sprintf('HCP_ses%d.mat',ses),'subjCode','Allrc_S2','Allrc_sigm_S1','Allrc_LG','Allrc_cvdeg','AC1','S1example','S2example','sigm_S1example','LGexample','cvdegexample','Npat','Meanrc_S2','SDrc_S2','AllS1','AllS2','AllNCL','Allpc_sp','Allpc_mom')
end %Ses