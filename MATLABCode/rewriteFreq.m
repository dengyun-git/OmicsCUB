%%%%%% count Sn occurrence, first need to produce *ForTbp.txt file, then
%%%%%% rewrite with right frequency

% load('orderName461.mat');
% speciesName={'saccharomyces_cerevisiae'};
% speciesName={'arthroderma_benhamiae_cbs_112371'};

fileName0='bacteriaNameList.csv';
fileID0=fopen(fileName0,'r');
speciesNamep=textscan(fileID0,'%s','Delimiter','\n');
speciesName=speciesNamep{1,1};
fclose(fileID0);

for tP=1:length(speciesName)
    
    fileName1=[speciesName{tP},'ForTbp.txt'];
    
    fileIDtemp=fopen(fileName1,'r');
    XX0=textscan(fileIDtemp,[repmat('%s ',1,8),'%s\n'],1,'Delimiter',',');
    XX=textscan(fileIDtemp,'%s %u %f %u %f %u %u %u %u\n','Delimiter',',');    %%%XX{1}:aa, XX{2}:sublength;XX{3}:Pi;XX{4}:DegPi;XX{5}:Pmax;XX{7}:frequency
    fclose(fileIDtemp);
    
    fmt1=[repmat('%s,',1,8),'%s\n'];
    fmt2='%s,%u,%d,%u,%d,%u,%u,%u,%u\n';
    
    aaList={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
    
    fileName2=[speciesName{tP},'ForBcTb.txt']; %%%%write to this temperature file codon table replacement
    fileID=fopen(fileName2,'a');
    fprintf(fileID,fmt1,'aa','L','Pi','DegPi','Pmax','DegPmax','Freq','GeneId','Codon1');
    
    for aa=1:18
        aaName=aaList{aa};
        aaID=find(ismember(XX{1},aaName));  %%%extract the same aa
        Lenw=XX{2}(aaID,1); %%%% lengw: length whole for one aa
        Pwp=XX{3}(aaID,1); %%%% Pwp: P whole for one aa
        DegPi1=XX{4}(aaID,1); %%%% DegPi1: Pi degeneracy whole for one aa
        Pmaxp=XX{5}(aaID,1);
        DegPmax1=XX{6}(aaID,1);
        GeneIDp=XX{8}(aaID,1);
        Codon1p=XX{9}(aaID,1);
        len=unique(Lenw);
        
        for interestLcount=1:length(len)
            interestL=len(interestLcount);
            interestID=find(Lenw==interestL);
            Pw=Pwp(interestID); %%%% Pw: P whole for a length within one aa
            DegPi2=DegPi1(interestID);
            Pmax=Pmaxp(interestID);
            DegPmax2=DegPmax1(interestID);
            GeneID=GeneIDp(interestID);
            Codon1=Codon1p(interestID);
            
            for i=1:length(Pw)
                %%%%Lenw: length whole
                Sn=Pw(i);
                Freq=sum(Pw==Sn);
                
                fprintf(fileID,fmt2,aaList{aa},interestL,Pw(i),DegPi2(i),Pmax(i),DegPmax2(i),Freq,GeneID(i),Codon1(i));
                
            end
        end
    end
    
    fclose(fileID);
    
end
