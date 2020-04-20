%%IMPORTANT don't forget to change 4 lines special for sc

%%%%%%% for Temperature raw data file: use original genome
%%%%%%% pi pmax calculation use codon usage table
%%%%%%%in the format file name: species+ForTb.txt
%%%%%%%inside each file: include degeneracy and geneID
global cfG ctG cfD ctD cfE ctE cfV ctV cfA ctA cfR ctR cfS ctS cfK ctK cfN ctN cfI ctI cfT ctT cfC ctC cfY ctY cfL ctL cfF ctF cfQ ctQ cfH ctH cfP ctP;

% speciesName={'saccharomyces_cerevisiae'};

% fileName0='protistNameList.csv';   %%% for protist list
fileName0='fungiNameListABC.csv';   %%% for fungi list
fileID0=fopen(fileName0,'r');
speciesNamep=textscan(fileID0,'%s','Delimiter','\n');
speciesName=speciesNamep{1,1};
fclose(fileID0);

fmt1=('%s %u %u %f\n');
fmt2=('%s %u %u %u %f\n');
fmt3=('%s %u %u %u %u %f\n');
fmt4=('%s %u %u %u %u %u %u %f\n');

fmt11=('%s %s %s %s\n');
fmt22=('%s %s %s %s %s\n');
fmt33=('%s %s %s %s %s %s\n');
fmt44=('%s %s %s %s %s %s %s %s\n');

fmList={fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt2,fmt3,fmt3,fmt3,fmt3,fmt3,fmt4,fmt4,fmt4};
fmList2={fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt22,fmt33,fmt33,fmt33,fmt33,fmt33,fmt44,fmt44,fmt44};

text={'GluE','HisH','GlnQ','PheF','TyrY','CysC','AsnN','LysK','AspD','IleI','ProP','ThrT','AlaA','ValV','GlyG','LeuL','SerS','ArgR'};

for tP=1:length(speciesName)
    %     setSynonymousCodonTable(tP,2,'synoTableProtist.txt');
    setSynonymousCodonTable(tP,2,'synoTableFungi.txt');
    
    species=speciesName{tP};
    %     fileName1=[speciesName{tP},'SUB.txt']; %%%% read Sn L files
    fileName2=[speciesName{tP},'ForFgTbp.txt']; %%%%write to this temperature file codon table replacement
    %     fileName2=[speciesName{tP},'ForPtTbp.txt'];
    %     fileName3=[speciesName{tP},'FFF.fa']; %%%%% read original genes files to get codon sequense
    %     fileName4=[speciesName{tP},'.cds.fa']; %%%%% for bacteria
    %     fileName4=[speciesName{tP},'.cds.pt.fa']; %%%%for protist
    fileName4=[speciesName{tP},'.fa']; %%%%for fungi
    
    cfWhole={cfE,cfH,cfQ,cfF,cfY,cfC,cfN,cfK,cfD,cfI,cfP,cfT,cfA,cfV,cfG,cfL,cfS,cfR};
    
    pasteCodon=getCodonSequence(fileName4);
    aaList={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
    synoL=[repmat(2,1,9),3,repmat(4,1,5),repmat(6,1,3)];
    funList={@GluAminoAcidH,@HisAminoAcidH,@GlnAminoAcidH,@PheAminoAcidH,@TyrAminoAcidH,@CysAminoAcidH,@AsnAminoAcidH,@LysAminoAcidH,@AspAminoAcidH,@IleAminoAcidH,@ProAminoAcidH,@ThrAminoAcidH,@AlaAminoAcidH,@ValAminoAcidH,@GlyAminoAcidH,@LeuAminoAcidH,@SerAminoAcidH,@ArgAminoAcidH};
    
    fileID=fopen(fileName2,'a');
    fmt1=[repmat('%s,',1,8),'%s\n'];
    fmt2='%s,%u,%d,%u,%d,%u,%u,%u,%u\n';
    
    fprintf(fileID,fmt1,'aa','L','Pi','DegPi','Pmax','DegPmax','Freq','GeneId','Codon1');
    
    %     XLLppre=dlmread(fileName1,',',1,0);
    
    for aa=1:18
        %%%%%%%%%%%%%
        cLeng=synoL(aa);
        %         XLLpre=unique(XLLppre(:,aa));   %%%%%% sepcial for sc
        % %         XLLpre=unique(XLLppre(:,18+aa)); %%%%% 461 species length column
        %XLL=XLLpre(~isnan(XLLpre));  %%%% L without overlap and delete NaN genes
        %%%%fungi
        
        %%%%%for bacteria get XXL
        fileName1=[speciesName{tP},text{aa},'ratioFg.txt']; %%%% for XXL
        %         fileName1=[speciesName{tP},text{aa},'ratioPt.txt']; %%%% for XXL
        fileID1 = fopen(fileName1);
        tt = textscan(fileID1,fmList2{aa},1,'delimiter',',');
        xx = textscan(fileID1,fmList{aa},'delimiter',',');
        fclose(fileID1);
        XLLppre=xx{1,end};
        XLLpre=unique(XLLppre);
        XLL=XLLpre(XLLpre~=0);
        %%%%%%%for bacteria XXL ends
        
        cfPart=cfWhole{aa}{:,:};
        
        for interestLid=1:length(XLL)
            
            interestL=XLL(interestLid);  %%% for each unique L
            
            %             Lid=getCertainLInd(aa,interestL,XLLppre);%%%% for fungi find all the certain length genes
            Lid=find(XLLppre==interestL);
            
            for i=1:length(Lid)
                
                codon=pasteCodon{1,Lid(i)};
                
                [~,Xg,Pi]=funList{aa}(codon);
                
                [Pmax,Xmax] = EforMore(cLeng,interestL,cfPart,2);
                
                %%%%%degeneracy too difficult need to think further
                DEG1=1;
                DEG2=1;
                Freq=1;
                
                fprintf(fileID,fmt2,aaList{aa},interestL,Pi,DEG1,Pmax,DEG2,Freq,Lid(i),Xg(1));
                
            end
            
        end
    end
    
    fclose(fileID);
end





% %%%%%% slope plots comparison
% %load('orderName461.mat');
% %tP=267;
% % fileName1='ForTempArti.txt';    %%sequences generation,Pi calculation and Pmax calculation all using equal frequencies;
% % fileName2='ForTempArtiB.txt';   %%sequences generation using codon table,Pi and Pmax calculation using equal frequencies;
% % fileName3='ForTemp.txt';
% %  fmt2=[repmat('%s ',1,8),'%s']; %%% write format the head of new file
% % fmt3=['%s %s %u ',repmat('%f ',1,5),'%f'];
% %
% % fileID1=fopen(fileName1,'r');
% % X0=textscan(fileID1,fmt2,1,'Delimiter',',');
% % X1=textscan(fileID1,fmt3,'Delimiter',',');
% % ID1=find(ismember(X1{1,1},speciesName{tP}));
% %
% % fileID2=fopen(fileName2,'r');
% % Y0=textscan(fileID2,fmt2,1,'Delimiter',',');
% % Y1=textscan(fileID2,fmt3,'Delimiter',',');
% % ID2=find(ismember(Y1{1,1},speciesName{tP}));
% %
% % fileID3=fopen(fileName3,'r');
% % Z0=textscan(fileID3,fmt2,1,'Delimiter',',');
% % Z1=textscan(fileID3,fmt3,'Delimiter',',');
% % ID3=find(ismember(Z1{1,1},speciesName{tP}));
% %
% % plot(X1{4}(ID1,:),'y.');
% % hold on
% % plot(Y1{4}(ID2,:),'c.');
% % plot(Z1{4}(ID3,:),'r.');
% %
% % hold off
% % title('slope of candida_tropicalis_mya_3404 comparison');
% % xlabel('18 aa in candida_tropicalis_mya_3404');
% % ylabel('slopes');
% % legend('equal artificial','codon table artificial','real sequences');
% % xticks(25:25:25*18);
% % xticklabels({'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'});
% % grid on