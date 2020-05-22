%%%%compare species name list from Dominique whether the same as mine.
setSynonymousCodonTable(1,1);


fileName33='scProtein.txt';
fileid33=fopen(fileName33);
XXX=textscan(fileid33,'%u %s %f %f','Delimiter','\t');
fclose(fileid33);

homoList=XXX{1,2};

% homoList={'YGL103W','YOL127W','YLR075W','YBL092W','YOR063W','YLR340W','YOR369C',...
%     'YDR064W','YOL040C','YGL123W','YPL267W','YIL158W','YCL029C','YNL166C','YOR026W',...
%     'YCL014W','YLR175W','YCR002C','YLR215C','YKL022C','YCL040W','YFR053C','YGL059W',...
%     'YGL253W','YGR192C','YHR044C','YJL155C','YKL127W','YMR278W','YNL241C','YJL200C',...
%     'YGR204W','YAL012W','YIL116W','YMR108W','YLR451W','YDL182W','YNL277W','YDR007W',...
%     'YOR323C','YFR028C','YLR272C','YBR135W','YMR199W','YAL040C','YJL164C','YGL240W'...
%     'YML027W','YDR201W','YDR293C'};

% for homoc=[47,560,825,607,623,675,721]
% for homoc=[47,560,623,607]%%%% YAL005C,YGL053W,YFL062W,YDR342C; para1, para3, para4.fig
%%% PARA COUNT 8 8 8 15 12 25
for j=1:50
    geneName=homoList{j};
    
    speciesName='saccharomyces_cerevisiae';
    
    fileName1=sprintf('%s',[speciesName,'.fa']);%% find corresponding species file
    fileName2=sprintf('%s',[speciesName,'Gene.fa']);%% only gene names in such file
    geneNamesWhole=getHomologyList(fileName2); %% list all the gene names in the interesed species
    
    
    id=find(contains(geneNamesWhole,geneName));
    pasteCodonWhole=getCodonSequence(fileName1); %%%%% find interested homology
    pasteCodon=pasteCodonWhole{id};
    
    
    %%%%%%%%%%%%%Glu%%%%%%
    [NNe,~,Ep]=GluAminoAcidH(pasteCodon);
    E(j)=-log(Ep)/NNe;
    %             if NNe<=800
    %                 PlocationE(j) = entropy4Location(2,NNe,E);
    %             else
    %                 PlocationE(j)=NaN;
    %             end
    
    %%%%%%%%%%%%%%His%%%%%%
    %                 [NNh,~,Hp]=HisAminoAcidH(pasteCodon);
    %                 H(j)=-log(Hp)/NNh;
    %                 %                 if NNh<=800
    %                 %                     Plocation = entropy4Location(2,NNh,H);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %
    %
    %                 %%%%%%%%%%%%%Gln%%%%%%
    %                 [NNq,~,Qp]=GlnAminoAcidH(pasteCodon);
    %                 Q(j)=-log(Qp)/NNq;
    %                 %                 if NNq<=800
    %                 %                     Plocation = entropy4Location(2,NNq,Q);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %
    %                 %%%%%%%%%%%%%Lys%%%%%%
    %                 [NNk,~,Kp]=LysAminoAcidH(pasteCodon);
    %                 K(j)=-log(Kp)/NNk;
    %                 %                 if NNk<=800
    %                 %                     Plocation = entropy4Location(2,NNk,K);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %
    %
    %                 %%%%%%%%%%%%%%Asp%%%%%%
    %                 [NNd,~,Dp]=AspAminoAcidH(pasteCodon);
    %                 D(j)=-log(Dp)/NNd;
    %                 %                 if NNd<=800
    %                 %                     Plocation = entropy4Location(2,NNd,D);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %                 %
    %
    %                 %%%%%%%%%%%%%Phe%%%%%%
    %                 [NNf,~,Fp]=PheAminoAcidH(pasteCodon);
    %                 %                 F=-log(Fp)/NNf;
    %                 %                 if NNf<=800
    %                 %                     Plocation = entropy4Location(2,NNf,F);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %
    %                 %%%%%%%%%%%Tyr%%%%%%
    %                 [NNy,~,Yp]=TyrAminoAcidH(pasteCodon);
    %                 Y=-log(Yp)/NNy;
    %                 %                 if NNy<=800
    %                 %                     Plocation = entropy4Location(2,NNy,Y);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %
    %                 %%%%%%%%%%%%Cys%%%%%%
    %                 [NNc,~,Cp]=CysAminoAcidH(pasteCodon);
    %                 C=-log(Cp)/NNc;
    %                 %                 if NNc<=800
    %                 %                     Plocation = entropy4Location(2,NNc,C);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %
    %
    %                 %%%%%%%%%%Asn%%%%%%
    %                 [NNn,~,Np]=AsnAminoAcidH(pasteCodon);
    %                 N=-log(Np)/NNn;
    %                 %                 if NNn<=800
    %                 %                     Plocation = entropy4Location(2,NNn,N);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    
    
    %%%%%%%%%%Ile%%%%%%
    [NNi,~,Ip]=IleAminoAcidH(pasteCodon);
    I(j)=-log(Ip)/NNi;
    %             if NNi<=400
    %                 PlocationI(j) = entropy4Location(3,NNi,I);
    %             else
    %                 PlocationI(j)=NaN;
    %             end
    
    
    %%%%%%%%%%%Gly%%%%%%%
    [NNg,~,Gp]=GlyAminoAcidH(pasteCodon);
    G(j)=-log(Gp)/NNg;
    %             if NNg<=300
    %                 PlocationG(j) = entropy4Location(4,NNg,G);
    %             else
    %                 PlocationG(j)=NaN;
    %             end
    %
    
    %%%%%%%%%%%%%Pro%%%%%%%
    %                 [NNp,~,Pp]=ProAminoAcidH(pasteCodon);
    %                 P=-log(Pp)/NNp;
    %                 %                 if NNp<=300
    %                 %                     Plocation = entropy4Location(4,NNp,P);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %                 %
    %
    %                 %%%%%%%%%%%Val%%%%%%%
    %                 [NNv,~,Vp]=ValAminoAcidH(pasteCodon);
    %                 %                 V=-log(Vp)/NNv;
    %                 %                 if NNv<=300
    %                 %                     Plocation = entropy4Location(4,NNv,V);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %
    %
    %                 %%%%%%%%%%%%Thr%%%%%%%
    %                 [NNt,~,Tp]=ThrAminoAcidH(pasteCodon);
    %                 T=-log(Tp)/NNt;
    %                 %                 if NNt<=300
    %                 %                     Plocation = entropy4Location(4,NNt,T);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %                 %
    %
    %                 %%%%%%%%%%%%Ala%%%%%%%
    %                 [NNa,~,Ap]=AlaAminoAcidH(pasteCodon);
    %                 A=-log(Ap)/NNa;
    %                 %                 if NNa<=300
    %                 %                     Plocation = entropy4Location(4,NNa,A);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %
    %
    %%%%%%%%%%%%%%Arg%%%%%%%
    [NNr,~,Rp]=ArgAminoAcidH(pasteCodon);
    R(j)=-log(Rp)/NNr;
    %             if NNr<=400
    %                 PlocationR(j) = entropy4Location(6,NNr,R);
    %             else
    %                 PlocationR(j)=NaN;
    %             end
    
    
    
    %               %%%%%%%%%%%%Ser%%%%%%%
    %                 [NNs,~,Sp]=SerAminoAcidH(pasteCodon);
    %                 S=-log(Sp)/NNs;
    %                 %                 if NNs<=400
    %                 %                     Plocation = entropy4Location(6,NNs,S);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %
    %
    %                 %%%%%%%%%%%%Leu%%%%%%%
    %                 [NNl,~,Lp]=LeuAminoAcidH(pasteCodon);
    %                 L=-log(Lp)/NNl;
    %                 %                 if NNl<=400
    %                 %                     Plocation = entropy4Location(6,NNl,L);
    %                 %                 else
    %                 %                     Plocation=NaN;
    %                 %                 end
    %
    
    
    protein(j)=XXX{1,3}(j);
    
    
    clear pasteCodon
    
    %             catch
    %                 display([num2str(j),':error']);
    %             end
    
    
    
end





cov11=corrcoef(protein,E);%%%%%correlationship Sn vs Protein
cov12=corrcoef(protein,I);
cov13=corrcoef(protein,G);
cov14=corrcoef(protein,R);

%     cov21=corrcoef(paraVl(disID),Plocation(disID)); %%%%%%%%%%correlationship A vs Protein
%     cov22=corrcoef(paraVl(disID),Plocation(disID));
%     cov23=corrcoef(paraVl(disID),Plocation(disID));
%     cov24=corrcoef(paraVl(disID),Plocation(disID));
%
    Fit11=fitlm(protein',E','linear','RobustOpts','on');
    Fit12=fitlm(protein',I,'linear','RobustOpts','on');
    Fit13=fitlm(protein',G','linear','RobustOpts','on');
    Fit14=fitlm(protein',R','linear','RobustOpts','on');

%     %         figure
%     %     plot(Fit11.Residuals.Raw,'r*');
%     %     ylabel('Residual');
%     %     title(['Residuals for E of: ',stringName]);
%     %
%     %     figure
%     %     plot(Fit12.Residuals.Raw,'g*');
%     %     ylabel('Residual');
%     %     title(['Residuals for I of: ',stringName]);
%     %
%     %     figure
%     %     plot(Fit13.Residuals.Raw,'b*');
%     %     ylabel('Residual');
%     %     title(['Residuals for G of: ',stringName]);
%     %
%     %     figure
%     %     plot(Fit14.Residuals.Raw,'y*');
%     %     ylabel('Residual');
%     %     title(['Residuals for R of: ',stringName]);
%
% %     Fit21=fit((paraVl(disID))',(PlocationE(disID))',referenceType);
% %     Fit22=fit((paraVl(disID))',(PlocationI(disID))',referenceType);
% %     Fit23=fit((paraVl(disID))',(PlocationG(disID))',referenceType);
% %     Fit24=fit((paraVl(disID))',(PlocationR(disID))',referenceType);
% %
        figure
        plot(protein,E,'r*');
        hold on
        plot(protein,Fit11.Coefficients.Estimate(2)*protein+Fit11.Coefficients.Estimate(1),'r-');
        plot(protein,I,'g*');
        plot(protein,Fit12.Coefficients.Estimate(2).*protein+Fit12.Coefficients.Estimate(1),'g-')
        plot(protein,G,'b*');
        plot(protein,Fit13.Coefficients.Estimate(2).*protein+Fit13.Coefficients.Estimate(1),'b-')
        plot(protein,R,'y*');
        plot(protein,Fit14.Coefficients.Estimate(2).*protein+Fit14.Coefficients.Estimate(1),'y-')
        hold off
    %
%         title(['Paralogous of ',homoList{homoc},': Protein Concentration vs Sn']);
%         xlabel('protein concentration');
%         ylabel('Sn');
%         legend('E','fitE','I','fitI','G','fitG','R','fitR');
% 
%         figure
%         plot(paraVl(disID),PlocationE(disID),'r.');
%         hold on
%         plot(paraVl(disID),PlocationI(disID),'g.');
%         plot(paraVl(disID),PlocationG(disID),'b.');
%         plot(paraVl(disID),PlocationR(disID),'y.');
%         hold off
%         title(['Paralogous of ',homoList{homoc},': Protein Concentration vs Sn']);
%         xlabel('protein concentration');
%         ylabel('Avalue');
%         legend('E','I','G','R');
% 
