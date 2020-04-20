%%%%%generate raw data for artificial sequence fit of Temprature calculation
%%%%%%%in the format file name: species+ForFgTa.txt
%%%%%%%inside each file: include degeneracy and geneID
global cfG ctG cfD ctD cfE ctE cfV ctV cfA ctA cfR ctR cfS ctS cfK ctK cfN ctN cfI ctI cfT ctT cfC ctC cfY ctY cfL ctL cfF ctF cfQ ctQ cfH ctH cfP ctP;

% speciesName={'saccharomyces_cerevisiae'};

%%%%%%for fungi speices list
fileName0='protistNameList.csv';   %%% for protist list
% fileName0='fungiNameListABC.csv';   %%% for fungi list
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


for tP=1:length(speciesName) %%% all species
    %     setSynonymousCodonTable(tP,2,'synoTableFungi.txt'); %%%%for fungi
    setSynonymousCodonTable(tP,2,'synoTableProtist.txt'); %%%%for protist
    cfWhole={cfE,cfH,cfQ,cfF,cfY,cfC,cfN,cfK,cfD,cfI,cfP,cfT,cfA,cfV,cfG,cfL,cfS,cfR};
    
    
    %     fileName1=[speciesName{tP},'SUB.txt']; %read Sn L files,early method
    %     fileName2=[speciesName{tP},'ForBcTabp.txt']; %%%%write to this temperature file codon table replacement
    %     fileName2=[speciesName{tP},'ForPtTabp.txt'];
    %     fileName2=['/home/cug/yd46/MATLAB/CseqsTab/',speciesName{tP},'ForPtTabp.txt'];
    fileName2=[speciesName{tP},'ForPtTabp.txt'];
    %     fileName3=[speciesName{tP},'FFF.fa']; %%%%% read original genes files to get codon sequense
    
    %     pasteCodon=getCodonSequence(fileName4);
    aaList={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
    synoL=[repmat(2,1,9),3,repmat(4,1,5),repmat(6,1,3)];
    
    fileID=fopen(fileName2,'a');
    fmt1=[repmat('%s,',1,8),'%s\n'];
    fmt2='%s,%u,%d,%u,%d,%u,%u,%u,%u\n';
    
    fprintf(fileID,fmt1,'aa','L','Pi','DegPi','Pmax','DegPmax','Freq','GeneId','Codon1');
    
    %     XLLppre=dlmread(fileName1,',',1,0);
    
    for aa=1:18
        cLeng=synoL(aa);
        cfPart=cfWhole{aa}{:,:};%%%%% prepare P for multinomial generator
        %
        
        %%%% fungi length information
        % %         XLLpre=unique(XLLppre(:,aa));   %%%%%% sepcial for sc
        % %         %         XLLpre=unique(XLLppre(:,18+aa)); %%%%% 461 species length column
        % %         XLL=XLLpre(~isnan(XLLpre));  %%%% L without overlap and delete NaN genes
        %%%% fungi length information end
        
        %%%%%fungi length information
        %         fileName1=[speciesName{tP},text{aa},'ratioPt.txt']; %%%% for XXL
        fileName1=[speciesName{tP},text{aa},'ratioPt.txt'];
        fileID1 = fopen(fileName1);
        tt = textscan(fileID1,fmList2{aa},1,'delimiter',',');
        xx = textscan(fileID1,fmList{aa},'delimiter',',');
        fclose(fileID1);
        XLLppre=xx{1,end};
        XLLpre=unique(XLLppre);
        XLL=XLLpre(XLLpre~=0);
        %%%%fungi length infomration ends
        
        for interestLid=1:length(XLL)
            
            interestL=XLL(interestLid);  %%% for each unique L
            
            %             %%%%%following function need choose sc or 461
            %             Lid=getCertainLInd(aa,interestL,XLLppre);%%%%find all the certain length genes
            Lid=find(XLLppre==interestL);
            Xg=mnrnd(interestL,cfPart,length(Lid)); %%%% generate ctb random configurations for Lid counts
            Pi=mnpdf(Xg,cfPart); %%%% Lid counts Pi (with duplicated ones)
            
            for i=1:length(Lid)
                
                %                 [Pmax,Xmax] = EforMore(cLeng,interestL,Psyno,2);%%%%use codon table frequency
                [Pmax,~] = EforMore(cLeng,interestL,cfPart);
                %%%%%degeneracy too difficult need to think further
                DEG1=1;
                DEG2=1;
                Freq=1;
                
                fprintf(fileID,fmt2,aaList{aa},interestL,Pi(i),DEG1,Pmax,DEG2,Freq,Lid(i),Xg(i,1));
                
                
            end
            
        end
    end
    
    fclose(fileID);
end
