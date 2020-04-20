%%%%%%%generate raw data for artificial sequence fit of Temprature calculation
%%%%%%%in the format file name: species+ForTa.txt
%%%%%%%inside each file: include degeneracy and geneID
global cfG ctG cfD ctD cfE ctE cfV ctV cfA ctA cfR ctR cfS ctS cfK ctK cfN ctN cfI ctI cfT ctT cfC ctC cfY ctY cfL ctL cfF ctF cfQ ctQ cfH ctH cfP ctP;

setSynonymousCodonTable(1,1);

cfWhole={cfE,cfH,cfQ,cfF,cfY,cfC,cfN,cfK,cfD,cfI,cfP,cfT,cfA,cfV,cfG,cfL,cfS,cfR};

% speciesName={'saccharomyces_cerevisiae'};

%%%%%below for bacteria list
fileName0='fungiNameListABC.csv'; %%%%for fungi
% fileName0='protistNameList.csv';   %%% for protist list
fileID0=fopen(fileName0,'r');
speciesNamep=textscan(fileID0,'%s','Delimiter','\n');
speciesName=speciesNamep{1,1};
fclose(fileID0);
%%%%fungi list end

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


for tP=1:length(speciesName) %%% all species
    %     fileName1=[speciesName{tP},'SUB.txt']; %%%% read Sn L files
    fileName2=['/home/cug/yd46/MATLAB/CseqsTab/',speciesName{tP},'ForFgTap.txt']; %%%%%write to this file
    %         fileName2=[speciesName{tP},'ForPtTap.txt']; %%%%%write to this file
    %     fileName3=[speciesName{tP},'FFF.fa']; %%%%% read original genes files to get codon sequense
    
    fileID=fopen(fileName2,'a');
    fmt1=[repmat('%s,',1,8),'%s\n'];
    fmt2='%s,%u,%d,%u,%d,%u,%u,%u,%u\n';
    
    fprintf(fileID,fmt1,'aa','L','Pi','DegPi','Pmax','DegPmax','Freq','GeneId','Codon1');
    
    aaList={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
    synoL=[repmat(2,1,9),3,repmat(4,1,5),repmat(6,1,3)];
    
    
    for aa=1:18
        cLeng=synoL(aa);
        cfPart=cfWhole{aa}{:,:};%%%%% prepare P for multinomial generator
        
        %%%%%% length information for fungi
        % %         %
        % %         XLLppre=dlmread('scSublength.txt',',',1,0);%%%%%species for sc
        % %                          XLLpre=unique(XLLppre(:,aa));   %%%%%% sepcial for sc
        % % %         XLLppre=dlmread(fileName1,',',1,0);
        % % %         XLLpre=unique(XLLppre(:,18+aa)); %%%%% 461 species length column
        % % %         XLL=XLLpre(~isnan(XLLpre));  %%%% L without overlap and delete NaN genes
        % length information for fungi end
        %
        
        %%%length information for Fungi
        fileName1=[speciesName{tP},text{aa},'ratioFg.txt']; %%%% for fungi XXL
        %         fileName1=[speciesName{tP},text{aa},'ratioPt.txt']; %%%% for protist XXL
        fileID1 = fopen(fileName1);
        tt = textscan(fileID1,fmList2{aa},1,'delimiter',',');
        xx = textscan(fileID1,fmList{aa},'delimiter',',');
        fclose(fileID1);
        XLLppre=xx{1,end};
        XLLpre=unique(XLLppre);
        XLL=XLLpre(XLLpre~=0);
        %%%length information for Fungi end
        
        
        for interestLid=1:length(XLL)
            
            interestL=XLL(interestLid);  %%% for each unique L
            
            %             Lid=getCertainLInd(aa,interestL,XLLppre);%%%% for fungi find all the certain length genes
            Lid=find(XLLppre==interestL);
            
            Xg=mnrnd(interestL,cfPart,length(Lid)); %%%% generate equal random configurations for Lid counts
            
            %             Xg=mnrnd(interestL,cfWhole{aa}{1,:},length(Lid)); %%%% generate ctb random configurations for Lid counts
            Pi=mnpdf(Xg,cfPart); %%%% Lid counts Pi (with duplicated ones)
            %                         Pi=mnpdf(Xg,cfWhole{aa}{1,:}); %%%% Lid counts Pi (with duplicated ones)
            
            for i=1:length(Lid)
                
                [Pmax,Xmax] = EforMore(cLeng,interestL,cfPart,1);
                %                 [Pmax,Xmax] = EforMore(cLeng,interestL,cfT{aa}{1,:});%%%%use codon table frequency
                
                DEG1=getDegereracy(Xg(i,:),cLeng);
                DEG2=getDegereracy(Xmax,cLeng);
                
                %                 Freq=getSameSnCount(i,Pi);
                %                 Svalue=log(Pmax./(Pi.*DEG1));
                %
                Freq=1;
                
                fprintf(fileID,fmt2,aaList{aa},interestL,Pi(i),DEG1,Pmax,DEG2,Freq,Lid(i),Xg(i,1));
                
                
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