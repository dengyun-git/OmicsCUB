%%%%%%%%generate KL2 file for each species
load('orderName461.mat');
% speciesName={'saccharomyces_cerevisiae'};
aaList={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};


for tP=1:length(speciesName)
    
    fileName1=[speciesName{tP},'ForTab.txt'];
    
    fileIDtemp=fopen(fileName1,'r');
    XX0=textscan(fileIDtemp,[repmat('%s ',1,7),'%s\n'],1,'Delimiter',',');
    XX=textscan(fileIDtemp,'%s %u %f %u %f %u %f %u\n','Delimiter',',');    %%%XX{1}:aa, XX{2}:sublength;XX{3}:Pi;XX{4}:DegPi;XX{5}:Pmax;XX{7}:frequency
    fclose(fileIDtemp);
    
    fileName2='aForKL2.txt'; %%%%write to this temperature file codon table replacement
    fileID=fopen(fileName2,'a');
%     fprintf(fileID,'%s,%s,%s,%s,%s\n','speices','aa','KL2','Nltp','N'); %%% Nltp: how many types of l length; N: data points
    
    for aa=1:18
        aaName=aaList{aa};
        aaID=find(ismember(XX{1},aaName));  %%%extract the same aa
        Len1=XX{2}(aaID,1); %%%% len1: length whole for one aa
        Pi1=XX{3}(aaID,1); %%%% Pi1: P whole for one aa
        Freq1=XX{7}(aaID,1); %%%% Freq1: frequency for one aa
        
        Nltp=length(unique(Len1));
        N=length(Len1);

        [Pi2p,row1,col1]=unique(Pi1); %%%%%Pi2 based on unique Pi values
        Pi2=Pi2p./Nltp;
        Freq2p=Freq1(row1);
        Freq2=Freq2p./sum(Freq2p);
                
        KL2=sum(Freq2.*log(Freq2./Pi2));
        
        fprintf(fileID,'%s,%s,%f,%u,%u\n',speciesName{tP},aaList{aa},KL2,Nltp,N);
        
    end
end

fclose(fileID);

