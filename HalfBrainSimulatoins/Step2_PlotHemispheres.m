% Code to generate Figs.4, 5, and SM2 in Gabaldon et al, 
% "Heuristic Inferencence of the brain dynamical state from the r-spectrum of the correlation matrix"

clear all
close all
clc


%% Figure 4
load('Results/ResultsHemi2_Tvalues25.mat') % Hemisphere=2 --> Full Brain
figure(400)

subplot(311)
plot(Tvec,Fa,'-o')

xlabel('T')
ylabel('<F_a>')

subplot(312)
plot(Tvec,AC1,'-o')
hold on
subplot(313)
plot(Tvec,lamb1,'-o')
xlabel('T')
ylabel('AC(1)')




mkdir("Figures/Fig4/DataDS")

% Export data
X=cat(2,Tvec',Fa');
save('Figures/Fig4/DataDS/Fa_vs_T.dat', 'X', '-ascii', '-double', '-tabs')

X=cat(2,Tvec',AC1');
save('Figures/Fig4/DataDS/AC1_vs_T.dat', 'X', '-ascii', '-double', '-tabs')


X=cat(2,Tvec',lamb1');
save('Figures/Fig4/DataDS/l1_vs_T.dat', 'X', '-ascii', '-double', '-tabs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
load('Results/ResultsHemi0_Tvalues91.mat')


minl=1
lim=91;
figure(500)
subplot(311)
plot(Tvec,Allrc_S2,'-ro'); hold on
plot(Tvec,Allrc_LG,'-ko')
plot(Tvec,Allrc_sigm_S1,'-bo')
xlabel('T')
ylabel('r_c')

subplot(312)
plot(AC1,Allrc_S2,'r.'); hold on
plot(AC1,Allrc_LG,'k.')
plot(AC1,Allrc_sigm_S1,'b.')
xlabel('AC(1)')
ylabel('r_c')

subplot(313)
plot(AllMMat,Allrc_S2,'r.'); hold on
plot(AllMMat,Allrc_sigm_S1,'b.')
plot(AllMMat,Allrc_LG,'k.')
xlabel('<Pairwise Correlation>')
ylabel('r_c')





% Export data

mkdir("Figures/Fig5/DataDS")

X=cat(2,Tvec',Allrc_S2',Allrc_LG',Allrc_sigm_S1');
save('Figures/Fig5/DataDS/rc_vs_T.dat', 'X', '-ascii', '-double', '-tabs')

X=cat(2,AC1',Allrc_S2',Allrc_LG',Allrc_sigm_S1');
save('Figures/Fig5/DataDS/rc_vs_AC1.dat', 'X', '-ascii', '-double', '-tabs')

X=cat(2,AllMMat',Allrc_S2',Allrc_LG',Allrc_sigm_S1');
save('Figures/Fig5/DataDS/rc_vs_MMat.dat', 'X', '-ascii', '-double', '-tabs')

%% Add Error bars with 5x resolution
clear all
load('Results/ResultsHemi0_Tvalues451.mat')

minl=1
lim=451;

NTV=mean(reshape(Tvec(minl:lim-1), [25 18]))

NrcS2=mean(reshape(Allrc_S2(minl:lim-1), [25 18]))
NrcS2sd=std(reshape(Allrc_S2(minl:lim-1), [25 18]))
Nrcsigm_S1=mean(reshape(Allrc_sigm_S1(minl:lim-1), [25 18]))
Nrcsigm_S1sd=std(reshape(Allrc_sigm_S1(minl:lim-1), [25 18]))
NrcLG=mean(reshape(Allrc_LG(minl:lim-1), [25 18]))
NrcLGsd=std(reshape(Allrc_LG(minl:lim-1), [25 18]))




figure(500)
subplot(311)
errorbar(NTV,NrcS2,NrcS2sd)
errorbar(NTV,Nrcsigm_S1,Nrcsigm_S1sd)
errorbar(NTV,NrcLG,NrcLGsd)


X=cat(2,NTV',NrcS2',NrcS2sd');
save('Figures/Fig5/DataDS/rc_vs_T_eb1.dat', 'X', '-ascii', '-double', '-tabs')
X=cat(2,NTV',Nrcsigm_S1',Nrcsigm_S1sd');
save('Figures/Fig5/DataDS/rc_vs_T_eb2.dat', 'X', '-ascii', '-double', '-tabs')
X=cat(2,NTV',NrcLG',NrcLGsd');
save('Figures/Fig5/DataDS/rc_vs_T_eb3.dat', 'X', '-ascii', '-double', '-tabs')



%
clear NrcS2 NrcS2sd Nrcsigm_S1 Nrcsigm_S1sd NrcLG NrcLGsd
dAC=0.002
ridx=0
for AC=0.96:dAC:0.976
ridx=ridx+1
idx=(AC1>AC) .* (AC1<AC+dAC);
NAC1(ridx)=AC+dAC/2;
NrcS2(ridx)=mean(Allrc_S2(find(idx)));
NrcS2sd(ridx)=std(Allrc_S2(find(idx)));

Nrcsigm_S1(ridx)=mean(Allrc_sigm_S1(find(idx)));
Nrcsigm_S1sd(ridx)=std(Allrc_sigm_S1(find(idx)));

NrcLG(ridx)=mean(Allrc_LG(find(idx)));
NrcLGsd(ridx)=std(Allrc_LG(find(idx)));

end


subplot(312)
errorbar(NAC1,NrcS2,NrcS2sd)
errorbar(NAC1,Nrcsigm_S1,Nrcsigm_S1sd)
errorbar(NAC1,NrcLG,NrcLGsd)


X=cat(2,NAC1',NrcS2',NrcS2sd');
save('Figures/Fig5/DataDS/rc_vs_AC1eb1.dat', 'X', '-ascii', '-double', '-tabs')
X=cat(2,NAC1',Nrcsigm_S1',Nrcsigm_S1sd');
save('Figures/Fig5/DataDS/rc_vs_AC1eb2.dat', 'X', '-ascii', '-double', '-tabs')
X=cat(2,NAC1',NrcLG',NrcLGsd');
save('Figures/Fig5/DataDS/rc_vs_AC1eb3.dat', 'X', '-ascii', '-double', '-tabs')




clear NrcS2 NrcS2sd Nrcsigm_S1 Nrcsigm_S1sd NrcLG NrcLGsd
dMM=0.01
ridx=0
for MM=0.03:dMM:0.09
ridx=ridx+1
idx=(AllMMat>MM) .* (AllMMat<MM+dMM);
NMMAT(ridx)=MM+dMM/2;
NrcS2(ridx)=mean(Allrc_S2(find(idx)));
NrcS2sd(ridx)=std(Allrc_S2(find(idx)));

Nrcsigm_S1(ridx)=mean(Allrc_sigm_S1(find(idx)));
Nrcsigm_S1sd(ridx)=std(Allrc_sigm_S1(find(idx)));

NrcLG(ridx)=mean(Allrc_LG(find(idx)));
NrcLGsd(ridx)=std(Allrc_LG(find(idx)));

end

X=cat(2,NMMAT',NrcS2',NrcS2sd');
save('Figures/Fig5/DataDS/rc_vs_MMateb1.dat', 'X', '-ascii', '-double', '-tabs')
X=cat(2,NMMAT',Nrcsigm_S1',Nrcsigm_S1sd');
save('Figures/Fig5/DataDS/rc_vs_MMateb2.dat', 'X', '-ascii', '-double', '-tabs')
X=cat(2,NMMAT',NrcLG',NrcLGsd');
save('Figures/Fig5/DataDS/rc_vs_MMateb3.dat', 'X', '-ascii', '-double', '-tabs')

subplot(313)
errorbar(NMMAT,NrcS2,NrcS2sd)
errorbar(NMMAT,Nrcsigm_S1,Nrcsigm_S1sd)
errorbar(NMMAT,NrcLG,NrcLGsd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig SM 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
load('Results/ResultsHemi1_Tvalues91.mat')

minl=1;
lim=91;
figure(1200)
subplot(311)
plot(Tvec,Allrc_S2,'r.'); hold on
plot(Tvec,Allrc_LG,'k.')
plot(Tvec,Allrc_sigm_S1,'b.')

xlabel('T')
ylabel('r_c')

subplot(312)
plot(AC1,Allrc_S2,'r.'); hold on
plot(AC1,Allrc_LG,'k.')
plot(AC1,Allrc_sigm_S1,'b.')

xlabel('AC(1)')
ylabel('r_c')

subplot(313)
plot(AllMMat,Allrc_S2,'r.'); hold on
plot(AllMMat,Allrc_sigm_S1,'b.')
plot(AllMMat,Allrc_LG,'k.')
xlabel('<Pairwise Correlation>')
ylabel('r_c')





% Export Data
mkdir("Figures/FigSM2/DataDS")

X=cat(2,Tvec',Allrc_S2',Allrc_LG',Allrc_sigm_S1');
save('Figures/FigSM2/DataDS/rc_vs_T.dat', 'X', '-ascii', '-double', '-tabs')

X=cat(2,AC1',Allrc_S2',Allrc_LG',Allrc_sigm_S1');
save('Figures/FigSM2/DataDS/rc_vs_AC1.dat', 'X', '-ascii', '-double', '-tabs')

X=cat(2,AllMMat',Allrc_S2',Allrc_LG',Allrc_sigm_S1');
save('Figures/FigSM2/DataDS/rc_vs_MMat.dat', 'X', '-ascii', '-double', '-tabs')
