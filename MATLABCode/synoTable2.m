
synoN=[2,2,2,2,2,2,2,2,2,3,4,4,4,4,4,6,6,6];
% text={ctE,ctH,ctQ,ctF,ctY,ctC,ctN,ctK,ctD,ctI,ctP,ctT,ctA,ctV,ctG,ctL,ctS,ctR};
text={'GluE','HisH','GlnQ','PheF','TyrY','CysC','AsnN','LysK','AspD','IleI','ProP','ThrT','AlaA','ValV','GlyG','LeuL','SerS','ArgR'};

% speciesName={'saccharomyces_cerevisiae'};
fileName0='bacteriaNameListABC.csv';
% fileName0='protistNameList.csv';
% fileName0='fungiNameListABC.csv';
fileID0=fopen(fileName0,'r');
speciesNamep=textscan(fileID0,'%s','Delimiter','\n');
speciesName=speciesNamep{1,1};
fclose(fileID0);

% fileName1='synoTableFungi.txt'; %%%wrtie to this file: fungi 462
fileName1='synoTableBacteria.txt'; %%% write to this file: bacteria
% fileName1='synoTableProtist.txt'; %%% write to this file: protist
fileID1=fopen(fileName1,'a');

fprintf(fileID1,[repmat('%s,',1,18),'%s\n'],'speciesName','E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R');

fmt1=('%s %u %u %u\n');
fmt2=('%s %u %u %u %u\n');
fmt3=('%s %u %u %u %u %u\n');
fmt4=('%s %u %u %u %u %u %u %u\n');

fmt11=('%s %s %s %s\n');
fmt22=('%s %s %s %s %s\n');
fmt33=('%s %s %s %s %s %s\n');
fmt44=('%s %s %s %s %s %s %s %s\n');

fmList={fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt2,fmt3,fmt3,fmt3,fmt3,fmt3,fmt4,fmt4,fmt4};
fmList2={fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt22,fmt33,fmt33,fmt33,fmt33,fmt33,fmt44,fmt44,fmt44};

for tP=1:length(speciesName)
    for CTcount=1:18
        
        fileName2=[speciesName{tP},text{CTcount},'ratioBc.txt'];
        fileID2 = fopen(fileName2);
        tt = textscan(fileID2,fmList2{CTcount},1,'delimiter',',');
        xx = textscan(fileID2,fmList{CTcount},'delimiter',',');
        fclose(fileID2);
        
        for xxcount=1:synoN(CTcount)
            XX(:,xxcount)=xx{1,xxcount+1};
        end
        
        cfSum=sum(XX);
        cfWholep{CTcount}=cfSum./sum(cfSum);
        
        clear XX
        
    end
    
    cfWhole=cell2mat(cfWholep);
    
    
    fprintf(fileID1,['%s,',repmat('%.17f,',1,58),'%.17f\n'],speciesName{tP},cfWhole);
    
    
end


fclose(fileID1);