% Code to generate  Figs. 1, 2, 3  SM1, SM5, SM7 and SM8  in Gabaldon et al,
% "Heuristic Inferencence of the brain dynamical state from the r-spectrum of the correlation matrix"

clear all
close all
clc

%%
rvals=[0:0.001:1]; % r values used to perform the r-spectrum
N_Rois=360;
load("HCP_ses1.mat")
mkdir("FiguresData")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%dS1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(100) %Figure 1 Panel D

subplot(411)
plot(rvals,S1example(1,:))
xlabel("r_c"); ylabel("S1"); xlim([0 0.9])

subplot(412)
plot(rvals,S2example(1,:))
xlabel("r_c"); ylabel("S2"); xlim([0 0.9])

subplot(413)
plot(rvals,sigm_S1example(1,:))
xlabel("r_c"); ylabel("\sigma S1"); xlim([0 0.9])

subplot(414)
plot(rvals,LGexample(1,:))
xlabel("r_c")
ylabel("L(G)")
xlim([0 0.9])

%Export Data For Fig. 1D
TT=table(rvals',S1example(1,:)'/N_Rois);
writetable(TT,'FiguresData/f1_S1_ses1_1.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

TT=table(rvals',S2example(1,:)'/N_Rois);
writetable(TT,'FiguresData/f1_S2_ses1_1.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

TT=table(rvals',sigm_S1example(1,:)');
writetable(TT,'FiguresData/f1_sigm_S1_ses1_1.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

LGexample_clean = LGexample; LGexample_clean(isnan(LGexample_clean)) = 0;
TT = table(rvals', LGexample_clean');
writetable(TT,'FiguresData/f1_LG_ses1_1.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(200) %Figure 2 left Panels

subplot(411)
plot(rvals,S1example(:,:))
xlabel("r_c"); ylabel("S1"); xlim([0 0.9])

subplot(412)
plot(rvals,S2example(:,:))
xlabel("r_c"); ylabel("S2"); xlim([0 0.9])

subplot(413)
plot(rvals,sigm_S1example(:,:))
xlabel("r_c"); ylabel("\sigma S1"); xlim([0 0.9])


subplot(414)
plot(rvals,LGexample(:,:))
xlabel("r_c"); ylabel("L(G)"); xlim([0 0.9]);

%Export Data Fig 2A
TT=table(rvals',S1example'/N_Rois);
writetable(TT,'FiguresData/f2_S1_ses1_10.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

TT=table(rvals',S2example'/N_Rois);
writetable(TT,'FiguresData/f2_ses1_10.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

TT=table(rvals',sigm_S1example');
writetable(TT,'FiguresData/f2_sigm_S1_ses1_10.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')


LGexample_clean = LGexample;
LGexample_clean(isnan(LGexample_clean)) = 0;
TT = table(rvals', LGexample_clean');
writetable(TT,'FiguresData/f2_LG_ses1_10.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(201) %Figure 2 right Panel

plot(AC1(1:64),Allrc_S2(1:64),'.'); hold on
plot(AC1(1:64),Allrc_sigm_S1(1:64),'.')
plot(AC1(1:64),Allrc_LG(1:64),'.')

xlabel("AC(1)")
ylabel("r_c")
legend("S2" , "\sigma S1", "L(G)")

%Export Data

TT=table(AC1(1:64)',Allrc_S2(1:64)',Allrc_sigm_S1(1:64)',Allrc_LG(1:64)');
writetable(TT,'FiguresData/f2B_rc_ses1_64.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(300) %Figure 3 Left

plot(AC1,Allrc_S2,'.'); hold on
plot(AC1,Allrc_sigm_S1,'.') %Not plotted in the manuscript
plot(AC1,Allrc_LG,'.')
xlabel("AC(1)")
ylabel("r_c")
legend("S2" , "\sigma S1", "L(G)")

% Export Data

TT=table(AC1',Allrc_S2',Allrc_sigm_S1',Allrc_LG');
writetable(TT,'FiguresData/f3_rc_ses1_all.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(301) %Figure 3 Right
errorbar([(22+25)/2 (26+30)/2 (31+35)/2 36 ], Meanrc_S2, SDrc_S2./sqrt(Npat-1), 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);

%Export Data

TT=table([(22+25)/2 (26+30)/2 (31+35)/2 36 ]', Meanrc_S2', SDrc_S2'./sqrt(Npat'-1));
writetable(TT,'FiguresData/f3B_Ages_ses1.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')
xlabel("Age (yr.)")
ylabel("r_c (S2)")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure SM1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
rvals=[0:0.001:1];  % r values used to perform the r-spectrum
load("HCP_ses5.mat")
N_Rois=360;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1100) %Figure SM1 up-left Panels

subplot(411)
plot(rvals,S1example(:,:))
xlabel("r_c")
ylabel("S1")
xlim([0 0.9])

subplot(412)
plot(rvals,S2example(:,:))
xlabel("r_c")
ylabel("S2")
xlim([0 0.9])

subplot(413)
plot(rvals,sigm_S1example(:,:))
xlabel("r_c")
ylabel("\sigma S1")
xlim([0 0.9])

subplot(414)
plot(rvals,LGexample(:,:))
xlim([0 0.9])
xlabel("r_c")
ylabel("L(G)")


%Export Data

TT=table(rvals',S1example'/N_Rois);
writetable(TT,'FiguresData/fS1A_S1_all10.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

TT=table(rvals',S2example'/N_Rois);
writetable(TT,'FiguresData/fS1A_S2_all10.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

TT=table(rvals',sigm_S1example');
writetable(TT,'FiguresData/fS1A_sigm_S1_all10.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

LGexample_clean = LGexample;
LGexample_clean(isnan(LGexample_clean)) = 0;
TT = table(rvals', LGexample_clean');
writetable(TT,'FiguresData/fS1A_LG_all10.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1101) %Figure SM1 up-right Panel
plot(AC1(1:64),Allrc_S2(1:64),'.') ;hold on
plot(AC1(1:64),Allrc_sigm_S1(1:64),'.')
plot(AC1(1:64),Allrc_LG(1:64),'.')

xlabel("AC(1)")
ylabel("r_c")
legend("S2" , "\sigma S1", "L(G)","cvd")

%Export Data

TT=table(AC1(1:64)',Allrc_S2(1:64)',Allrc_sigm_S1(1:64)',Allrc_LG(1:64)');
writetable(TT,'FiguresData/fS1B_rc_all64.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1102) %Figure SM1 down-left panel

plot(AC1,Allrc_S2,'.'); hold on
plot(AC1,Allrc_sigm_S1,'.') %Not plotted in the manuscript
plot(AC1,Allrc_LG,'.')

xlabel("AC(1)")
ylabel("r_c")
legend("S2" , "\delta S1", "L(G)")

% Export Data

TT=table(AC1',Allrc_S2',Allrc_sigm_S1',Allrc_LG');
writetable(TT,'FiguresData/fS1CA_rcLG_all_all.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1103)  %Figure SM1 down-right panel
errorbar([(22+25)/2 (26+30)/2 (31+35)/2 36 ], Meanrc_S2, SDrc_S2./sqrt(Npat-1), 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel("Age (yr.)")
ylabel("r_c (S2)")

% Export Data

TT=table([(22+25)/2 (26+30)/2 (31+35)/2 36 ]', Meanrc_S2', SDrc_S2'./sqrt(Npat'-1));
writetable(TT,'FiguresData/fS1D_Ages_all.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure SM5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("HCP_ses1.mat")
figure(1500)
for i=1:100
    plot(meandeg(i,:),AllS1(i,:))
    %plot(AllS1(i,:),sddeg(i,:)./meandeg(i,:))
    hold on
end


figure(1501)
for i=1:size(sddeg,1)
    if i<20
        plot(rvals,sddeg(i,:)./meandeg(i,:))
        %plot(rvals,sddeg(i,:))
    end
    [a,b]=max(sddeg(i,:)./meandeg(i,:))
    Allrc_sddegNORM(i)=rvals(b);
    [a,b]=max(sddeg(i,:))
    Allrc_sddeg(i)=rvals(b);
    hold on
end

TT=table(rvals',sddeg(1:10,:)');
writetable(TT,'FiguresData/fS5A_sddeg.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

TT=table(rvals',(sddeg(1:10,:)./meandeg(1:10,:))');
writetable(TT,'FiguresData/fS5A_sddegnorm.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')


figure(182)
for i=1:50
    subplot(221)
    plot(rvals/Allrc_S2(i),AllS1(i,:))
    hold on
    subplot(222)
    plot(rvals/Allrc_sigm_S1(i),AllS1(i,:))
    hold on
    subplot(223)
    plot(rvals/Allrc_LG(i),AllS1(i,:))
    hold on
    subplot(224)
    plot(rvals/Allrc_sddegNORM(i),AllS1(i,:))
    hold on
end

figure(183)
plot(Allrc_sddeg,Allrc_LG, '.')
hold on
plot(Allrc_sddeg,Allrc_S2, '.')
plot(Allrc_sddeg,Allrc_sigm_S1, '.')


figure(184)
plot(Allrc_sddegNORM,Allrc_LG, '.')
hold on
plot(Allrc_sddegNORM,Allrc_S2, '.')
plot(Allrc_sddegNORM,Allrc_sigm_S1, '.')



TT=table(Allrc_sddeg(1:64)',Allrc_S2(1:64)',Allrc_sigm_S1(1:64)',Allrc_LG(1:64)');
writetable(TT,'FiguresData/fS5A_sddegCORREL.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')


TT=table(Allrc_sddegNORM(1:64)',Allrc_S2(1:64)',Allrc_sigm_S1(1:64)',Allrc_LG(1:64)');
writetable(TT,'FiguresData/fS5A_sddegCORREL_NORM.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure SM7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
load('HCP_ses1.mat',"Allrc_S2",'Allpc_sp','Allpc_mom')


figure(1700) 
plot(Allrc_S2,Allpc_sp,'.')
hold on
plot(Allrc_S2,Allpc_mom,'.')
xlabel("r_c [S_2]")
ylabel("p_c")

figure(1701) %Inset
loglog(Allrc_S2,Allpc_sp,'.')
hold on
loglog(Allrc_S2,Allpc_mom,'.')
grid on

% Export Data


TT=table(Allrc_S2',Allpc_sp',Allpc_mom');
writetable(TT,'FiguresData/fS7_pcrc.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure SM8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
rvals=[0:0.001:1];  % r values used to perform the r-spectrum
load("HCP_ses1.mat")
N_Rois=360;


for i=1:996 %Compute the relative r value (i.e., r/r_c for each subject)
    normrc(i,:)=rvals/Allrc_LG(i);
end

%Compute mean and SD S1 for different subjects with the same relative r value
dr=0.005
ridx=0
for r=10*dr:dr:1.5
    ridx=ridx+1;
    Idx= (normrc>r) .* (normrc<r+dr); %Group subjects values with the same relative r value
    
    mS1(ridx)=mean(AllS1(logical(Idx)));
    sS1(ridx)=std(AllS1(logical(Idx)));
end

figure(1800)
subplot(211)
for i=1:10
    plot(normrc(i,:),AllS1(i,:))
    hold on
end
xlabel("r/r_c  [L(G)]")
ylabel("S_1")



subplot(212)
plot(10*dr:dr:1.5,sS1)
grid on
xlabel("r/r_c  [L(G)]")
ylabel("s.d. (S_1"))

% Export Data

TT=table(normrc(1:10,:)',AllS1(1:10,:)'/N_Rois)
writetable(TT,'FiguresData/fS8A_S1_all.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')


TT=table([10*dr:dr:1.5]', mS1'/N_Rois,sS1'/N_Rois);
writetable(TT,'FiguresData/fS8B_s_all.dat','WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')
