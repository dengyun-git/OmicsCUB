%%%%% find all the wrong genes
fileName0='protistNameList.csv'; 
fileID0=fopen(fileName0,'r');
speciesListp=textscan(fileID0,'%s','Delimiter','\n');
speciesList=speciesListp{1,1};
fclose(fileID0);

filename1='wrongGeneRecordPt.csv';
fileID1=fopen(filename1,'a');
fprintf(fileID1,'%s,%s,%s\n','speciesName','geneNumber','WrongGene');


for tP=1:length(speciesList)
    
    speciesName=speciesList{tP};
    
    ctWrong=0;
    
    filename=[speciesName,'.cds.pt.fa'];
    fileID = fopen(filename);
    tt = textscan(fileID,'%s','delimiter',';');
    fclose(fileID);
    
    xxp = tt{1};
    xx= strrep(xxp,char(39),'');
    l = length(xx);
    datRaw = xx(2:2:l,:);
    
    for i=1:length(datRaw)
        
        dat=datRaw{i};
        
        if ~mod(length(dat),3)==0
            
            ctWrong=ctWrong+1;
        end
        
        clear dat
        
    end
    
    fprintf(fileID1,'%s,%u,%u\n',speciesName,length(datRaw),ctWrong);   
end

fclose(fileID1);

