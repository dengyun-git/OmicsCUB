%%%%%%%%generate KL1 file for each species
load('orderName461.mat');
% speciesName={'wallemia_ichthyophaga_exf_994'};

aaList={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};

for tP=1:length(speciesName)
    
    fileName1=[speciesName{tP},'ForT.txt'];
    
    fileIDtemp=fopen(fileName1,'r');
    XX0=textscan(fileIDtemp,[repmat('%s ',1,7),'%s\n'],1,'Delimiter',',');
    XX=textscan(fileIDtemp,'%s %u %f %u %f %u %f %u\n','Delimiter',',');    %%%XX{1}:aa, XX{2}:sublength;XX{3}:Pi;XX{4}:DegPi;XX{5}:Pmax;XX{7}:frequency
    fclose(fileIDtemp);
            
    fileName2=[speciesName{tP},'ForTKL1.txt']; %%%%write to this temperature file codon table replacement
    fileID=fopen(fileName2,'a');
    fprintf(fileID,'%s,%s,%s,%s\n','aa','L','KL1','Nl');
    
    for aa=1:18
        aaName=aaList{aa};
        aaID=find(ismember(XX{1},aaName));  %%%extract the same aa
        Len1=XX{2}(aaID,1); %%%% len1: length whole for one aa
        Pi1=XX{3}(aaID,1); %%%% Pi1: P whole for one aa
        Freq1=XX{7}(aaID,1); %%%% Freq1: frequency for one aa
        lenUni=unique(Len1); %%%% find different lengths
        
        for interestLcount=1:length(lenUni)
            interestL=lenUni(interestLcount);
            interestID=find(Len1==interestL); %%%%interestID1 based on repeated Pi values
            Nl=length(interestID);

            Pi2=Pi1(interestID);    %%%%%Pi2 and Freq2 is for one lengths but repeated p values
            Freq2=Freq1(interestID);
            
            [Pi3,row1,col1]=unique(Pi2); %%%% Pi3 based on unique Pi values
            Freq3p=Freq2(row1);
            Freq3=Freq3p./sum(Freq3p);
            
            KL1=sum(Freq3.*log(Freq3./Pi3));
            
            fprintf(fileID,'%s,%u,%f,%u\n',aaList{aa},interestL,KL1,Nl);
                
        end
    end
    
    fclose(fileID);
    
end
