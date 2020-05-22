setSynonymousCodonTable(1,1);

% % GeneInterestList=getHomologyList('saccharomyces_arboricola_h_6Gene.fa');%%%%generate homoList which contain interested genes
% % indChosen=randi(3659,[500,1]);
% % homoListSel=GeneInterestList(indChosen,:);

homoList={'YGL103W','YOL127W','YLR075W','YBL092W','YOR063W','YLR340W','YOR369C',...
    'YDR064W','YOL040C','YGL123W','YPL267W','YIL158W','YCL029C','YNL166C','YOR026W',...
    'YCL014W','YLR175W','YCR002C','YLR215C','YKL022C','YCL040W','YFR053C','YGL059W',...
    'YGL253W','YGR192C','YHR044C','YJL155C','YKL127W','YMR278W','YNL241C','YJL200C',...
    'YGR204W','YAL012W','YIL116W','YMR108W','YLR451W','YDL182W','YNL277W','YDR007W',...
    'YOR323C','YFR028C','YLR272C','YBR135W','YMR199W','YAL040C','YJL164C','YGL240W'...
    'YML027W','YDR201W','YDR293C'};   %%homology sourse


% fileID=fopen('homoList5818Rand.txt');
% tt=textscan(fileID,'%s','Delimiter',',');
% homoList=tt{1,1};
% fclose(fileID);

% % % homoList={'YBR221C'};

% for homoc=1:length(homoList)

for homoc=1:1 %%homology counter
    
    stringName=homoList{homoc};
    
    [homoSpecies,homoGene,l] = getHomoInfor(stringName);%%%call function 'getHomoInfor' to retrieve homology information
    
    if l~=0
        
        
        for j=1:l
            try
                speciesName=homoSpecies{j};
                fileName1=sprintf('%s',[speciesName,'.fa']);%% find corresponding species file
                fileName2=sprintf('%s',[speciesName,'Gene.fa']);%% only gene names in such file
                geneNamesWhole=getHomologyList(fileName2); %% list all the gene names in the interesed species
                geneName=homoGene{j};
                
                id=find(contains(geneNamesWhole,geneName));
                pasteCodonWhole=getCodonSequence(fileName1); %%%%% find interested homology
                pasteCodon=pasteCodonWhole{id};
                
                %                 fileName3=[speciesName,'CAI.txt'];
                fileName3='PDI_substrates.csv';
                fileid3=fopen(fileName3);
                temp=textscan(fileid3,'%s %s %s %s',1,'Delimiter',',');
                %                 temp=textscan(fileid3,'%s %s',1,'Delimiter',',');
                %  XcaiW=textscan(fileid3,'%s %f','Delimiter',',');
                XX=textscan(fileid3,'%s %f %f %f','Delimiter',',');
                fclose(fileid3);
                
%                 Xcai(j)=XcaiW{1,2}(id);
                %%%%%%%%%%%%%Glu%%%%%%
                [NNe,~,Ep]=GluAminoAcidH(pasteCodon);
                E(j)=-log(Ep)/NNe;
                if NNe<=800
                    PlocationE(j) = entropy4Location(2,NNe,E);
                else
                    PlocationE(j)=NaN;
                end
                
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
                if NNi<=400
                    PlocationI(j) = entropy4Location(3,NNi,I);
                else
                    PlocationI(j)=NaN;
                end
                
                
                %%%%%%%%%%%Gly%%%%%%%
                [NNg,~,Gp]=GlyAminoAcidH(pasteCodon);
                G(j)=-log(Gp)/NNg;
                if NNg<=300
                    PlocationG(j) = entropy4Location(4,NNg,G);
                else
                    PlocationG(j)=NaN;
                end
                
                
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
                if NNr<=400
                    PlocationR(j) = entropy4Location(6,NNr,R);
                else
                    PlocationR(j)=NaN;
                end
                
                
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
            catch
                display('error');
            end
            
            clear pasteCodon
            
        end
    end
    figure
    plot(Xcai,E,'r*');
    hold on
    plot(Xcai,I,'g*');
    plot(Xcai,G,'b*');
    plot(Xcai,R,'y*');
    hold off
    
    figure
    plot(Xcai,PlocationE,'r*');
    hold on
    plot(Xcai,PlocationI,'g*');
    plot(Xcai,PlocationG,'b*');
    plot(Xcai,PlocationR,'y*');
    hold off
    
    clear E I G R PlocationE PlocationI PlocationG PlocationR;
    
end


