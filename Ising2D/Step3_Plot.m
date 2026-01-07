% Code to generate Fig. SM4 in Gabaldon et al, 
% "Data-driven inference of brain dynamical states from the r-spectrum of correlation matrices"

clear all
close all
clc
Tvec=unique(sort([1.2:0.2:2 2.05:0.05:2.5 2.6:0.2:3.4 ])); % Vector of values of T: from 1.2 to 3.4 . Zoomed about 2.25
rvals=[0:0.001:1]; %r values used to perform the r-spectrum
%% 
mkdir("FiguresData")
L=64;  % Side lenght of the LxL grid
for run=1:5
tindex=0
for T=Tvec
    tindex=tindex+1
    load(sprintf('Analysis/ResultsR%d_L%d_T%f.mat',run,L,T),'mag','AC1','rc_S2','rc_sigm_S1','rc_LG','rc_Skew')
    tic
    mag=mag(fix(end/10)+1:end);
    ACx=autocorr(mag);
    AC1=ACx(2);
    
    
    % figure(200+tindex) % If desired plot the magnetization for each T
    % plot(mag)
    % hold on

    AllAC1(run,tindex)=AC1;
    AllMag(run,tindex)=mean(abs(mag))/(L*L); % Order parameter of the system
    Allrc_S2(run,tindex)=rc_S2;
    Allrc_sigm_S1(run,tindex)=rc_sigm_S1;
    Allrc_LG(run,tindex)=rc_Skew;
    Allrc_Skew(run,tindex)=rc_Skew;
 
    

 
        
        
        
        
    end
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. SM4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_run=run;
figure(1400)

subplot(221)
hold on
errorbar(Tvec,mean(abs(AllMag)),std(abs(AllMag))/sqrt(tot_run))
xlabel('T')
ylabel('|Mag|')

subplot(222)
hold on
errorbar(Tvec,mean(AllAC1),std(AllAC1)/sqrt(tot_run))
ylabel('AC(1)')
xlabel('T')

subplot(223)
hold on
errorbar(Tvec,mean(Allrc_S2),std(Allrc_S2)/sqrt(tot_run))
errorbar(Tvec,mean(Allrc_sigm_S1),std(Allrc_sigm_S1)/sqrt(tot_run))
errorbar(Tvec,mean(Allrc_LG),std(Allrc_LG)/sqrt(tot_run))
ylabel('rc')
xlabel('T')
legend({'S2','\sigma S1','L(G)'})

subplot(224)
hold on
errorbar(mean(AllAC1),mean(Allrc_S2),std(Allrc_S2)/sqrt(tot_run),'-o')
ylabel('rc S2')
xlabel('AC(1)')


% Export data
%Panel A
X=cat(2,Tvec',mean(AllMag,1)',std(AllMag,1)'/sqrt(tot_run));
save('FiguresData/Fa_vs_T.dat', 'X', '-ascii', '-double', '-tabs')
%Panel B
X=cat(2,Tvec',mean(AllAC1,1)',std(AllAC1,1)'/sqrt(tot_run));
save('FiguresData/AC1_vs_T.dat', 'X', '-ascii', '-double', '-tabs')
%Panel C
X=cat(2,Tvec',mean(Allrc_S2,1)',std(Allrc_S2,1)'/sqrt(tot_run));
save('FiguresData/rc_S2_vs_T.dat', 'X', '-ascii', '-double', '-tabs')

X=cat(2,Tvec',mean(Allrc_sigm_S1,1)',std(Allrc_sigm_S1,1)'/sqrt(tot_run));
save('FiguresData/rc_sigmS1_vs_T.dat', 'X', '-ascii', '-double', '-tabs')

X=cat(2,Tvec',mean(Allrc_LG,1)',std(Allrc_LG,1)'/sqrt(tot_run));
save('FiguresData/rc_LG_vs_T.dat', 'X', '-ascii', '-double', '-tabs')
%Panel D

X=cat(2,mean(AllAC1(:,1:6))',mean(Allrc_S2(:,1:6))');
save('FiguresData/AC_vs_rcS2SUB.dat', 'X', '-ascii', '-double', '-tabs')

X=cat(2,mean(AllAC1(:,7:end))',mean(Allrc_S2(:,7:end))');
save('FiguresData/AC_vs_rcS2SUPER.dat', 'X', '-ascii', '-double', '-tabs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. SM6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except rvals L Allrc_S2 Allrc_sigm_S1 Allrc_LG  Allrc_Skew tot_run Tvec

figure(1604) %Panel E
errorbar(Tvec,mean(Allrc_Skew),std(Allrc_Skew,1)'/sqrt(tot_run))

xlabel("T")
ylabel("r_c [Skew]")
X=cat(2,Tvec',mean(Allrc_Skew,1)',std(Allrc_Skew,1)'/sqrt(tot_run));
save('FiguresData/rc_Skew_vs_T.dat', 'X', '-ascii', '-double', '-tabs')


figure(1605) %Panel F
plot(mean(Allrc_Skew),mean(Allrc_S2),'o')
xlabel("r_c [Skew]")
ylabel("r_c [S2]")

X=cat(2,mean(Allrc_Skew,1)',mean(Allrc_S2,1)');
save('FiguresData/rc_Skew_vs_rc_S2.dat', 'X', '-ascii', '-double', '-tabs')


clear Tvec
%Now compute Panels A-D 
Tvec= [1.4 2.3 3] %Run this code for T=1.4 (Fig SM6-Panel A),  T=2.3 (Fig SM6-Panel B) and
%T=3 (Fig SM6-Panel C) 

tindex=0;
    for T=Tvec
Cluslist1=[];
Cluslist2=[];
Cluslist3=[];

        tindex=tindex+1;
        
for run=[1:5] %1:tot_run


       load(sprintf('Analysis/ResultsR%d_L%d_T%f.mat',run,L,T),'AllSizesCell','Skew','rc_Skew')    
AllSkewness(run,:)=Skew;
  
        Allrc_Skew(run,tindex)=rc_Skew;

%Select the value of r that gives power-law like distribution of cluster
%sizes. Take a window of size drx to improve statistics.
drx=5
rccl=685 %T=2.3
if T==1.4
rccl=93 %T=1.4
elseif T==3
rccl=416
end

%Take the list of clusters at this r_c
for r=rccl-drx:rccl+drx
   Cluslist1=[Cluslist1 AllSizesCell{r}] ;
end

%Take the list of clusters for r>r_c (r_c + 20 units)
for r=rccl+20-drx:rccl+20+drx
   Cluslist2=[Cluslist2 AllSizesCell{r}] ;
end

%Take the list of clusters for r<r_c (r_c - 20 units)
for r=rccl-20-drx:rccl-20+drx
   Cluslist3=[Cluslist3 AllSizesCell{r}] ;
end

    end  %run




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ClusSizes = 2.^(1:0.5:14);
%Compute the histograms
figure(1)
A = histogram(Cluslist1,ClusSizes,'Normalization','probability'); 
figure(1600+tindex-1)
loglog(A.BinEdges(1:end-1)*sqrt(sqrt(2)),A.Values,'-o','LineWidth',2)
hold on
grid on

dist1=A.Values;

figure(1)
A = histogram(Cluslist2,ClusSizes,'Normalization','probability'); 
figure(1600+tindex-1)
loglog(A.BinEdges(1:end-1)*sqrt(sqrt(2)),A.Values,'-o','LineWidth',2)


dist2=A.Values;

figure(1)
A = histogram(Cluslist3,ClusSizes,'Normalization','probability'); 
figure(1600+tindex-1)
loglog(A.BinEdges(1:end-1)*sqrt(sqrt(2)),A.Values,'-o','LineWidth',2)

xlabel("s")
ylabel("P(s)")
legend("r_c","r>r_c","r<r_c")


ActualClusSiz=A.BinEdges(1:end-1)*sqrt(sqrt(2));
dist3=A.Values;

X=cat(2,ActualClusSiz',dist1',dist2',dist3');
save(sprintf('FiguresData/ClustSizeDistT%d.dat',T), 'X', '-ascii', '-double', '-tabs')

figure(1603) %Panel D
m=mean(AllSkewness);
idx=~isnan(m);

plot(rvals(idx),m(idx))
hold on
X=cat(2,rvals(idx)',m(idx)');
save(sprintf('FiguresData/SkewT%d.dat',T), 'X', '-ascii', '-double', '-tabs')

    end % Tvec

    











