load('orderName461.mat');
% speciesName={'saccharomyces_cerevisiae'};

for tP=1:length(speciesName)
    
    fileName1=[speciesName{tP},'CAI.txt']; %%%wrtie to this file
    fileID1=fopen(fileName1,'a');
    fprintf(fileID1,'%s,%s\n','GeneName','CAI');
    
    [RefCAI,ObsCAI] = getRefCAI(speciesName{tP});
    
    for i=1:length(ObsCAI(:,1))
    tempC=ObsCAI(i,:);
    CAI(i)=prod(RefCAI.^tempC);
    fprintf(fileID1,'%u,%f\n',i,CAI(i));
    end 
    
    fclose(fileID1);
end