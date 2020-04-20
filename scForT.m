global cfG ctG cfD ctD cfE ctE cfV ctV cfA ctA cfR ctR cfS ctS cfK ctK cfN ctN cfI ctI cfT ctT cfC ctC cfY ctY cfL ctL cfF ctF cfQ ctQ cfH ctH cfP ctP;

%%%%%%%generate raw data for real sequence fit of Temprature calculatio
% load('orderName461.mat'); %%%% for fungi
% speciesName='saccharomyces_cerevisiae'; %%%%special for sc
% for tP=[14,116,307] %%%for in order 'Fusarium graminearum -  (fg)' 'Aspergillus niger -  (anr)' 'Schizzosaccharomyces pombe -  (sp)'

% fileName0='bacteriaNameListABC.csv';  %%%%for bacteria list
% fileName0='protistNameList.csv';   %%% for protist list
fileName0='fungiNameListABC.csv'; %%%%for fungi
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

setSynonymousCodonTable(1,1);
cfWhole={cfE,cfH,cfQ,cfF,cfY,cfC,cfN,cfK,cfD,cfI,cfP,cfT,cfA,cfV,cfG,cfL,cfS,cfR};


for tP=1:length(speciesName) %%% all species
    %     fileName1=[speciesName{tP},'SUB.txt']; %%%% for Fungi read Sn L files
        fileName2=[speciesName{tP},'ForFgTp.txt']; %%%%write to this temperature file
%     fileName2=[speciesName{tP},'ForPtTp.txt'];
        fileName3=[speciesName{tP},'.fa']; %%%%% for Fungi read original genes files to get codon sequense
%     fileName3=[speciesName{tP},'.cds.pt.fa']; %%%% for protist
    
    species=speciesName{tP};
    
    pasteCodon=getCodonSequence(fileName3);
    
    fileID=fopen(fileName2,'a');
    fmt1=[repmat('%s,',1,8),'%s\n'];
    fmt2='%s,%u,%d,%u,%d,%u,%u,%u,%u\n';
    
    fprintf(fileID,fmt1,'aa','L','Pi','DegPi','Pmax','DegPmax','Freq','GeneId','Codon1');
    
    aaList={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
    synoL=[repmat(2,1,9),3,repmat(4,1,5),repmat(6,1,3)];
    funList={@GluAminoAcidH,@HisAminoAcidH,@GlnAminoAcidH,@PheAminoAcidH,@TyrAminoAcidH,@CysAminoAcidH,@AsnAminoAcidH,@LysAminoAcidH,@AspAminoAcidH,@IleAminoAcidH,@ProAminoAcidH,@ThrAminoAcidH,@AlaAminoAcidH,@ValAminoAcidH,@GlyAminoAcidH,@LeuAminoAcidH,@SerAminoAcidH,@ArgAminoAcidH};
    
    for aa=1:18
        cLeng=synoL(aa);
        %%%%%%%%special for fungi
        %         % XLLppre=dlmread('scSublength.txt',',',1,0);%%%%%species for sc
        %         XLLppre=dlmread(fileName1,',',1,0);
        %         XLLpre=unique(XLLppre(:,18+aa)); %%%%% 461 species length column
        %         %         XLLpre=unique(XLLppre(:,aa));%%%%%%special for sc
        %         XLL=XLLpre(~isnan(XLLpre));
        %%%%%%%%for fungi ends
        
        %%%%%for Fungi get XXL
        fileName1=[speciesName{tP},text{aa},'ratioFg.txt']; %%%% for XXL
%                 fileName1=[speciesName{tP},text{aa},'ratioPt.txt']; %%%% for XXL
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
            
            interestL=XLL(interestLid);
            
            %             Lid=getCertainLInd(aa,interestL,XLLppre);%%%%for fungi find certain length genes
            
            Lid=find(XLLppre==interestL); %%%%%for fungi
            
            for i=1:length(Lid)
                
                codon=pasteCodon{1,Lid(i)};
                
                [~,Xg,Pi]=funList{aa}(codon);
                
                cfPart=zeros(1,cLeng);
                cfPart(1,:)=1/cLeng;
                [Pmax,Xmax] = EforMore(cLeng,interestL,cfPart,1);
                
                DEG1=getDegereracy(Xg,cLeng);
                DEG2=getDegereracy(Xmax,cLeng);
                Freq=1;
                %                 Freq=getSameSnCount(Lid,aa,i,XLLppre); %%%%%%% use previous SUB file, Sn not match Pi. So need rewriteFreq.m to calulate Freq
                
                fprintf(fileID,fmt2,aaList{aa},interestL,Pi,DEG1,Pmax,DEG2,Freq,Lid(i),Xg(1));
                
            end
            
        end
    end
    
    fclose(fileID);
end

