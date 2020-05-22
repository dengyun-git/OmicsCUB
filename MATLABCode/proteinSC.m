fileName1='scProtein2.txt';
fileid1=fopen(fileName1);
XX=textscan(fileid1,'%u %s %f %f','Delimiter','\t');
fclose(fileid1);

speciesName='saccharomyces_cerevisiae';
pro0=XX{1,3};
selID0=find(pro0>10); %%% include all the protein abundance regions
selID1=find(pro0>10000); %%% here select protein abundance regions
selID2=find(pro0<5);
prow=pro0(selID0);
proh=pro0(selID1);
prol=pro0(selID2);

geneWList=XX{1,2}(selID0); %%%% H: high gene, L: low gene, W: whole gene
geneHList=XX{1,2}(selID1);
geneLList=XX{1,2}(selID2);

fileName2=sprintf('%s',[speciesName,'Gene.fa']);%% only gene names in such file
geneNamesWhole=getHomologyList(fileName2); %% list all the gene names in the interesed species

for ctw=1:length(geneWList)
    idW(ctw)=find(contains(geneNamesWhole,geneWList{ctw}),1);
end

for cth=1:length(geneHList)
    idH(cth)=find(contains(geneNamesWhole,geneHList{cth}),1);
end

for ctl=1:length(geneLList)
    idL(ctl)=find(contains(geneNamesWhole,geneLList{ctl}),1);
end

fileName3='saccharomyces_cerevisiaeSUB.txt';
fileid3=fopen(fileName3);
fmt1=[repmat('%s ',1,36),'%s\n'];
fmt2=['%u ',repmat('%f ',1,35),'%d\n'];
temp=textscan(fileid3,fmt1,1,'Delimiter',',');
YY=textscan(fileid3,fmt2,'Delimiter',',');
fclose(fileid3);

AminoAcid={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};

for fk=10
%for fk=[2,3,4,6,7,8,9,10,13,14,15,16,17,18]
%for fk=[5,11,12,19] %% amino acid F P I R
    YYa=YY{1,fk};
    YYb=YY{1,fk+18};
%         figure %%%%protein abundance to Sn and length individually
%         yyaxis left
%         plot(prow,YYa(idW,:),'b.')
%         ylabel('Sn value','FontSize',14);
%     
%         yyaxis right
%         plot(prow, YYb(idW,:),'r.');
%         ylabel('subsequence length','FontSize',14);
%     
%         xlabel('protein abundance','FontSize',14);
%         
%         legend('Sn value','subsequence length');
%     %
%          title(['Protein abundance relationship analysis for amino acid ', AminoAcid{fk-1},' in S.cerevisiae']);
    %     %
    
    figure   %%%%% Sn to length in low high groups
    yyaxis left
    plot(YYa(idH,:),YYb(idH,:),'c.');
    ylabel('subsquence length','FontSize',14);

    yyaxis right
    plot(YYa(idL,:),YYb(idL,:),'r.');
    ylabel('subsquence length','FontSize',14);
    
    legend('gene group for high protein abundance','gene group for low protein abundance');
    xlabel('Sn value');
    title(['gene groups for lowest and highest protein abundance based Sn for amino acid: ',AminoAcid{fk-1}]);
    
%         figure
%         plot(YYa(idH,:),YYb(idH,:),'c.');
%         xlabel('Sn value');
%         ylabel('subsquence length');
%         title(['gene groups for highest protein abundance based Sn for amino acid: ',AminoAcid{fk-1}]);
%     
%     figure
%     plot(YYa(idL,:),YYb(idL,:),'r.');
%     xlabel('Sn value');
%     ylabel('subsquence length');
%     title(['gene groups for lowest protein abundance based Sn for amino acid: ',AminoAcid{fk-1}]);
%     
    %     figure
    %     scatter3(YYa(idW,:),YYb(idW,:),prow,'G.');
    %     xlabel('Sn value');
    %     ylabel('subsequence length');
    %     zlabel('protein abundance');
    %     title(['correlation analysis between protein and CUB based on Sn for amino acid: ',AminoAcid{fk-1}]);
    %
    
end

%cov11=corrcoef(xxx,E(disID));%%%%%correlationship Sn vs Protein
%         cov12=corrcoef(xxx,I(disID));
%         cov13=corrcoef(xxx,G(disID));
%         cov14=corrcoef(xxx,R(disID));

% Fit11=fitlm(xxx',(E(disID))','linear','RobustOpts','on');
%     Fit12=fitlm(xxx',(I(disID))','linear','RobustOpts','on');
%     Fit13=fitlm(xxx',(G(disID))','linear','RobustOpts','on');
%     Fit14=fitlm(xxx',(R(disID))','linear','RobustOpts','on');
% figure
%         plot(xxx,E(disID),'r*');
%         hold on
%         plot(xxx,Fit11.Coefficients.Estimate(2)*xxx+Fit11.Coefficients.Estimate(1),'r-');
%         plot(xxx,I(disID),'g*');
%         plot(xxx,Fit12.Coefficients.Estimate(2).*xxx+Fit12.Coefficients.Estimate(1),'g-')
%         plot(xxx,G(disID),'b*');
%         plot(xxx,Fit13.Coefficients.Estimate(2).*xxx+Fit13.Coefficients.Estimate(1),'b-')
%         plot(xxx,R(disID),'y*');
%         plot(xxx,Fit14.Coefficients.Estimate(2).*xxx+Fit14.Coefficients.Estimate(1),'y-')
%         hold off
% %     %
%         title(['Paralogous of ',homoList{homoc},': Protein Concentration vs Sn']);
%         xlabel('protein concentration');
%         ylabel('Sn');
%         legend('E','fitE','I','fitI','G','fitG','R','fitR');
%      figure
%     %     plot(Fit11.Residuals.Raw,'r*');
%     %     ylabel('Residual');
%     %     title(['Residuals for E of: ',stringName]);
%     %