setSynonymousCodonTable2('scFF.fa');

global cfG cfGT ctG cfD cfDT ctD cfE cfET ctE cfV cfVT ctV cfA cfAT ctA cfR cfRT ctR cfS cfST ctS cfK cfKT ctK cfN cfNT ctN cfM cfMT ctM cfI cfIT ctI cfT cfTT ctT cfW cfWT ctW cfZ cfZT ctZ cfC cfCT ctC cfY cfYT  ctY cfL cfLT ctL cfF cfFT ctF cfQ cfQT ctQ cfH cfHT ctH cfP cfPT ctP;

load('AveEntropy2f.mat');
load('AveEntropy3f.mat');
load('AveEntropy4f.mat');
load('AveEntropy6f.mat');

% fileID=fopen('scSUBxy.txt','a'); %%%%I:original;Ir:expected entropy;Ia:equal artificial;Iab: biased artificial.
% fprintf(fileID,'%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s;%s,%s,%s,%s\n','E','Er','Ea','Eab','H','Hr','Ha','Hab','Q','Qr','Qa','Qab','F','Fr','Fa','Fab','Y','Yr','Ya','Yab','C','Cr','Ca','Cab','N','Nr','Na','Nab','K','Kr','Ka','Kab','D','Dr','Da','Dab','I','Ir','Ia','Iab','P','Pr','Pa','Pab','T','Tr','Ta','Tab','A','Ar','Aa','Aab','V','Vr','Va','Vab','G','Gr','Ga','Gab','L','Lr','La','Lab','S','Sr','Sa','Sab','R','Rr','Ra','Rab');
% fclose(fileID);


% for tP=1:length(speciesName)
    
%fileName=sprintf('%s',[speciesName{tP},'.fa']);

pasteCodon=getCodonSequence('saccharomyces_cerevisiae.csv');
% load('pasteCodon42.mat');
% pasteCodon=pasteCodon42;


%%%%%%Ile begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector        
[NNi,~,Ip]=IleAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNi)||NNi>1000) 
Irp=AveEntropy3(NNi);
Ir(i,:)=-Irp/NNi;
I(i,:)=-log(Ip)/NNi;
Sia=replaceE(NNi,ctI);
[~,~,Iap]=IleAminoAcidH(Sia);
Ia(i,:)=-log(Iap)/NNi;
Siab=replaceT(NNi,ctI,cfIT);
[~,~,Iabp]=IleAminoAcidH(Siab);
Iab(i,:)=-log(Iabp)/NNi;
else 
I(i,:)=NaN;   
Ir(i,:)=NaN;
Ia(i,:)=NaN;
Iab(i,:)=NaN;
end 
end
%   


%%%%% Asp begins 
for i=1:length(pasteCodon)
[NNd,~,Dp]=AspAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNd)||NNd>1500) 
Drp=AveEntropy2(NNd);
Dr(i,:)=-Drp/NNd;
D(i,:)=-log(Dp)/NNd;
Sda=replaceE(NNd,ctD);
[~,~,Dap]=AspAminoAcidH(Sda);
Da(i,:)=-log(Dap)/NNd;
Sdab=replaceT(NNd,ctD,cfDT);
[~,~,Dabp]=AspAminoAcidH(Sdab);
Dab(i,:)=-log(Dabp)/NNd;
else 
D(i,:)=NaN;   
Dr(i,:)=NaN;
Da(i,:)=NaN;   
Dab(i,:)=NaN;
end 
end

%%%%%%%Glu begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector        
[NNe,~,Ep]=GluAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNe)||NNe>1500) 
Erp=AveEntropy2(NNe);
Er(i,:)=-Erp/NNe;
E(i,:)=-log(Ep)/NNe;
Sea=replaceE(NNe,ctE);
[~,~,Eap]=GluAminoAcidH(Sea);
Ea(i,:)=-log(Eap)/NNe;
Seab=replaceT(NNe,ctE,cfET);
[~,~,Eabp]=GluAminoAcidH(Seab);
Eab(i,:)=-log(Eabp)/NNe;
else 
E(i,:)=NaN;   
Er(i,:)=NaN;
Ea(i,:)=NaN;   
Eab(i,:)=NaN;
end 
end


%%% His begins
for i=1:length(pasteCodon)
[NNh,~,Hp]=HisAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNh)||NNh>1500)
Hrp=AveEntropy2(NNh);
Hr(i,:)=-Hrp/NNh;
H(i,:)=-log(Hp)/NNh;
Sha=replaceE(NNh,ctH);
[~,~,Hap]=HisAminoAcidH(Sha);
Ha(i,:)=-log(Hap)/NNh;
Shab=replaceT(NNh,ctH,cfHT);
[~,~,Habp]=HisAminoAcidH(Shab);
Hab(i,:)=-log(Habp)/NNh;
else 
H(i,:)=NaN;    
Hr(i,:)=NaN;
Ha(i,:)=NaN;    
Hab(i,:)=NaN;
end 
end

%%% Pro begins
for i=1:length(pasteCodon)
[NNp,~,Pp]=ProAminoAcidH(pasteCodon{1,i});
if ~(isnan(NNp)||NNp>700) 
Prp=AveEntropy4(NNp);
Pr(i,:)=-Prp/NNp;
P(i,:)=-log(Pp)/NNp;
Spa=replaceE(NNp,ctP);
[~,~,Pap]=ProAminoAcidH(Spa);
Pa(i,:)=-log(Pap)/NNp;
Spab=replaceT(NNp,ctP,cfPT);
[~,~,Pabp]=ProAminoAcidH(Spab);
Pab(i,:)=-log(Pabp)/NNp;
else 
P(i,:)=NaN;   
Pr(i,:)=NaN;
Pa(i,:)=NaN;   
Pab(i,:)=NaN;
end 
end


%%% Gly begins
for i=1:length(pasteCodon)  
[NNg,~,Gp]=GlyAminoAcidH(pasteCodon{1,i}); 
if ~(isnan(NNg)||NNg>700)
Grp=AveEntropy4(NNg);
Gr(i,:)=-Grp/NNg;
G(i,:)=-log(Gp)/NNg;
Sga=replaceE(NNg,ctG);
[~,~,Gap]=GlyAminoAcidH(Sga);
Ga(i,:)=-log(Gap)/NNg;
Sgab=replaceT(NNg,ctG,cfGT);
[~,~,Gabp]=GlyAminoAcidH(Sgab);
Gab(i,:)=-log(Gabp)/NNg;
else 
G(i,:)=NaN;    
Gr(i,:)=NaN;
Ga(i,:)=NaN;    
Gab(i,:)=NaN;
end 
end


%%% Val begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector        
[NNv,~,Vp]=ValAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNv)||NNv>700) 
Vrp=AveEntropy4(NNv);
Vr(i,:)=-Vrp/NNv;
V(i,:)=-log(Vp)/NNv;
Sva=replaceE(NNv,ctV);
[~,~,Vap]=ValAminoAcidH(Sva);
Va(i,:)=-log(Vap)/NNv;
Svab=replaceT(NNv,ctV,cfVT);
[~,~,Vabp]=ValAminoAcidH(Svab);
Vab(i,:)=-log(Vabp)/NNv;
else 
V(i,:)=NaN;   
Vr(i,:)=NaN;
Va(i,:)=NaN;   
Vab(i,:)=NaN;
end 
end


%%% Lys begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector        
[NNk,~,Kp]=LysAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNk)||NNk>1500) 
Krp=AveEntropy2(NNk);
Kr(i,:)=-Krp/NNk;
K(i,:)=-log(Kp)/NNk;
Ska=replaceE(NNk,ctK);
[~,~,Kap]=LysAminoAcidH(Ska);
Ka(i,:)=-log(Kap)/NNk;
Skab=replaceT(NNk,ctK,cfKT);
[~,~,Kabp]=LysAminoAcidH(Skab);
Kab(i,:)=-log(Kabp)/NNk;
else 
K(i,:)=NaN;   
Kr(i,:)=NaN;
Ka(i,:)=NaN;   
Kab(i,:)=NaN;
end 
end


%%% Asn begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector   
[NNn,~,Np]=AsnAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNn)||NNn>1500) 
Nrp=AveEntropy2(NNn);
Nr(i,:)=-Nrp/NNn;
N(i,:)=-log(Np)/NNn;
Sna=replaceE(NNn,ctN);
[~,~,Nap]=AsnAminoAcidH(Sna);
Na(i,:)=-log(Nap)/NNn;
Snab=replaceT(NNn,ctN,cfNT);
[~,~,Nabp]=AsnAminoAcidH(Snab);
Nab(i,:)=-log(Nabp)/NNn;
else 
N(i,:)=NaN;   
Nr(i,:)=NaN;
Na(i,:)=NaN;   
Nab(i,:)=NaN;
end 
end


%%%%% Thr begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector     
[NNt,~,Tp]=ThrAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNt)||NNt>700) 
Trp=AveEntropy4(NNt);
Tr(i,:)=-Trp/NNt;
T(i,:)=-log(Tp)/NNt;
Sta=replaceE(NNt,ctT);
[~,~,Tap]=ThrAminoAcidH(Sta);
Ta(i,:)=-log(Tap)/NNt;
Stab=replaceT(NNt,ctT,cfTT);
[~,~,Tabp]=ThrAminoAcidH(Stab);
Tab(i,:)=-log(Tabp)/NNt;
else 
T(i,:)=NaN;   
Tr(i,:)=NaN;
Ta(i,:)=NaN;   
Tab(i,:)=NaN;
end 
end


%%%% Cys begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector      
[NNc,~,Cp]=CysAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNc)||NNc>1500) 
Crp=AveEntropy2(NNc);
Cr(i,:)=-Crp/NNc;
C(i,:)=-log(Cp)/NNc;
Sca=replaceE(NNc,ctC);
[~,~,Cap]=CysAminoAcidH(Sca);
Ca(i,:)=-log(Cap)/NNc;
Scab=replaceT(NNc,ctC,cfCT);
[~,~,Cabp]=CysAminoAcidH(Scab);
Cab(i,:)=-log(Cabp)/NNc;
else 
C(i,:)=NaN;   
Cr(i,:)=NaN;
Ca(i,:)=NaN;   
Cab(i,:)=NaN;
end 
end

%%%%% Tyr beins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector   
[NNy,~,Yp]=TyrAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNy)||NNy>1500) 
Yrp=AveEntropy2(NNy);
Yr(i,:)=-Yrp/NNy;
Y(i,:)=-log(Yp)/NNy;
Sya=replaceE(NNy,ctY);
[~,~,Yap]=TyrAminoAcidH(Sya);
Ya(i,:)=-log(Yap)/NNy;
Syab=replaceT(NNy,ctY,cfYT);
[~,~,Yabp]=TyrAminoAcidH(Syab);
Yab(i,:)=-log(Yabp)/NNy;
else 
Y(i,:)=NaN;   
Yr(i,:)=NaN;
Ya(i,:)=NaN;   
Yab(i,:)=NaN;
end 
end


%%%%%Phe begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector    
[NNf,~,Fp]=PheAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNf)||NNf>1500) 
Frp=AveEntropy2(NNf);
Fr(i,:)=-Frp/NNf;
F(i,:)=-log(Fp)/NNf;
Sfa=replaceE(NNf,ctF);
[~,~,Fap]=PheAminoAcidH(Sfa);
Fa(i,:)=-log(Fap)/NNf;
Sfab=replaceT(NNf,ctF,cfFT);
[~,~,Fabp]=PheAminoAcidH(Sfab);
Fab(i,:)=-log(Fabp)/NNf;
else 
F(i,:)=NaN;   
Fr(i,:)=NaN;
Fa(i,:)=NaN;   
Fab(i,:)=NaN;
end 
end



%%%%%Gln begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector  
[NNq,~,Qp]=GlnAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNq)||NNq>1500) 
Qrp=AveEntropy2(NNq);
Qr(i,:)=-Qrp/NNq;
Q(i,:)=-log(Qp)/NNq;
Sqa=replaceE(NNq,ctQ);
[~,~,Qap]=GlnAminoAcidH(Sqa);
Qa(i,:)=-log(Qap)/NNq;
Sqab=replaceT(NNq,ctQ,cfQT);
[~,~,Qabp]=GlnAminoAcidH(Sqab);
Qab(i,:)=-log(Qabp)/NNq;
else 
Q(i,:)=NaN;   
Qr(i,:)=NaN;
Qa(i,:)=NaN;   
Qab(i,:)=NaN;
end 
end


%%%% Ala begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector   
[NNa,~,Ap]=AlaAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNa)||NNa>700) 
Arp=AveEntropy4(NNa);
Ar(i,:)=-Arp/NNa;
A(i,:)=-log(Ap)/NNa;
Saa=replaceE(NNa,ctA);
[~,~,Aap]=AlaAminoAcidH(Saa);
Aa(i,:)=-log(Aap)/NNa;
Saab=replaceT(NNa,ctA,cfAT);
[~,~,Aabp]=AlaAminoAcidH(Saab);
Aab(i,:)=-log(Aabp)/NNa;
else 
A(i,:)=NaN;   
Ar(i,:)=NaN;
Aa(i,:)=NaN;   
Aab(i,:)=NaN;
end 
end


%%%%%Arg begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector  
[NNr,~,Rp]=ArgAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNr)||NNr>400) 
Rrp=AveEntropy6f(NNr);
Rr(i,:)=-Rrp/NNr;
R(i,:)=-log(Rp)/NNr;
Sra=replaceE(NNr,ctR);
[~,~,Rap]=ArgAminoAcidH(Sra);
Ra(i,:)=-log(Rap)/NNr;
Srab=replaceT(NNr,ctR,cfRT);
[~,~,Rabp]=ArgAminoAcidH(Srab);
Rab(i,:)=-log(Rabp)/NNr;
else 
R(i,:)=NaN;   
Rr(i,:)=NaN;
Ra(i,:)=NaN;   
Rab(i,:)=NaN;
end 
end


%%%%Leu begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector      
[NNl,~,Lp]=LeuAminoAcidH(pasteCodon{1,i}); 
if ~(isnan(NNl)||NNl>400) 
Lrp=AveEntropy6f(NNl);
Lr(i,:)=-Lrp/NNl;
L(i,:)=-log(Lp)/NNl;
Sla=replaceE(NNl,ctL);
[~,~,Lap]=LeuAminoAcidH(Sla);
La(i,:)=-log(Lap)/NNl;
Slab=replaceT(NNl,ctL,cfLT);
[~,~,Labp]=LeuAminoAcidH(Slab);
Lab(i,:)=-log(Labp)/NNl;
else 
L(i,:)=NaN;   
Lr(i,:)=NaN;
La(i,:)=NaN;   
Lab(i,:)=NaN;
end 
end


%%%%%%Ser begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector   
[NNs,~,Sp]=SerAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNs)||NNs>400) 
Srp=AveEntropy6f(NNs);
Sr(i,:)=-Srp/NNs;
S(i,:)=-log(Sp)/NNs;
Ssa=replaceE(NNs,ctS);
[~,~,Sap]=SerAminoAcidH(Ssa);
Sa(i,:)=-log(Sap)/NNs;
Ssab=replaceT(NNs,ctS,cfST);
[~,~,Sabp]=SerAminoAcidH(Ssab);
Sab(i,:)=-log(Sabp)/NNs;
else 
S(i,:)=NaN;   
Sr(i,:)=NaN;
Sa(i,:)=NaN;   
Sab(i,:)=NaN;
end 
end

% dlmwrite('scSUBSnxy.txt',[E,Er,Ea,Eab,H,Hr,Ha,Hab,Q,Qr,Qa,Qab,F,Fr,Fa,Fab,Y,Yr,Ya,Yab,C,Cr,Ca,Cab,N,Nr,Na,Nab,K,Kr,Ka,Kab,D,Dr,Da,Dab,I,Ir,Ia,Iab,P,Pr,Pa,Pab,T,Tr,Ta,Tab,A,Ar,Aa,Aab,V,Vr,Va,Vab,G,Gr,Ga,Gab,L,Lr,La,Lab,S,Sr,Sa,Sab,R,Rr,Ra,Rab],'-append');
% dlmwrite('scSUBoriginal.txt',[E,H,Q,F,Y,C,N,K,D,I,P,T,A,V,G,L,S,R],'-append');
% dlmwrite('scSUBequalArti.txt',[Ea,Ha,Qa,Fa,Ya,Ca,Na,Ka,Da,Ia,Pa,Ta,Aa,Va,Ga,La,Sa,Ra],'-append');
% dlmwrite('scSUBbiasArti.txt',[Eab,Hab,Qab,Fab,Yab,Cab,Nab,Kab,Dab,Iab,Pab,Tab,Aab,Vab,Gab,Lab,Sab,Rab],'-append');
% dlmwrite('scSUBexpect.txt',[Er,Hr,Qr,Fr,Yr,Cr,Nr,Kr,Dr,Ir,Pr,Tr,Ar,Vr,Gr,Lr,Sr,Rr],'-append');
% 
% dlmwrite('scSUBoriginal.txt',[E,H,Q,F,Y,C,N,K,D,I,P,T,A,V,G,L,S,R],'delimiter',' ');
% dlmwrite('scSUBequalArti.txt',[Ea,Ha,Qa,Fa,Ya,Ca,Na,Ka,Da,Ia,Pa,Ta,Aa,Va,Ga,La,Sa,Ra],'delimiter',' ');
% dlmwrite('scSUBbiasArti.txt',[Eab,Hab,Qab,Fab,Yab,Cab,Nab,Kab,Dab,Iab,Pab,Tab,Aab,Vab,Gab,Lab,Sab,Rab],'delimiter',' ');
% dlmwrite('scSUBexpect.txt',[Er,Hr,Qr,Fr,Yr,Cr,Nr,Kr,Dr,Ir,Pr,Tr,Ar,Vr,Gr,Lr,Sr,Rr],'delimiter',' ');

%%%only used to generate plot
XX=dlmread('scSUBSnxy.txt',',',1,0);
text={'Glu','His','Gln','Phe','Tyr','Cys','Asn','Lys','Asp','Ile','Pro','Thr','Ala','Val','Gly','Leu','Ser','Arg'};

% for i=[2,3,4,5,6,7,8,9,11,12,13,14,16,17] %%%1,10,15,18 for main content
for i=[1,10,15,18]
    j=4*(i-1);
    figure
plot(XX(:,1+j),XX(:,2+j),'r.');
hold on;
plot(XX(:,4+j),XX(:,2+j),'g.');
plot(XX(:,3+j),XX(:,2+j),'y.');
x=[0,0.05,0.3];y=x;plot(x,y)
hold off
legend('Original','Biased Artificial','Equal Artificial','LINE: expected Sn=observed Sn');
xlim([0 1.4]);
xlabel('Empirical Sn','FontSize',14);
ylabel('Theorized Expected Sn','FontSize',14);
title(['Comparison between observed sequences and artificial sequences (',text{i},') in S.cerevisea'],'FontSize',12);
end
% savefig('A&Ofg.fig');