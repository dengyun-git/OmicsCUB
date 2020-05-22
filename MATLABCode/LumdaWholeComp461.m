%%%%%%among whole genome: 461 lambda fit on S, comparison between real genome & artificial ones
%%%%%%pay attention to the order of speciesName, orderName461 are ordered
%%%%%%according to pylogenetic tree; speciesName are ordered by alphabat
%%%%%%SUBL461 are by alphabat

load('speciesName.mat');

fmt1=[repmat('%s,',1,72),'%s\n'];
fmt2=['%s,',repmat('%d,',1,71),'%d\n'];

fileID=fopen('LumdaWholeComp461.txt','a'); %%%%I:original;Ir:expected entropy;Ia:equal artificial;Iab: biased artificial.
fprintf(fileID,fmt1,'Species','Ealambda','Halambda','Qalambda','Falambda','Yalambda','Calambda','Nalambda','Kalambda','Dalambda','Ialambda','Palambda','Talambda','Aalambda','Valambda','Galambda','Lalambda','Salambda','Ralambda','EaStd','HaStd','QaStd','FaStd','YaStd','CaStd','NaStd','KaStd','DaStd','IaStd','PaStd','TaStd','AaStd','VaStd','GaStd','LaStd','SaStd','RaStd','Eablambda','Hablambda','Qablambda','Fablambda','Yablambda','Cablambda','Nablambda','Kablambda','Dablambda','Iablambda','Pablambda','Tablambda','Aablambda','Vablambda','Gablambda','Lablambda','Sablambda','Rablambda','EabStd','HabStd','QabStd','FabStd','YabStd','CabStd','NabStd','KabStd','DabStd','IabStd','PabStd','TabStd','AabStd','VabStd','GabStd','LabStd','SabStd','RabStd');


for tP=1:length(speciesName)
    seleLumda=zeros(1,18);
    fileName=[speciesName{tP},'SUB.txt'];
    XX=dlmread(fileName,',',1,18);
    
    text={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
    
    cf=cell(1,18);
    [cf{1},cf{2},cf{3},cf{4},cf{5},cf{6},cf{7},cf{8},cf{9},cf{10},cf{11},cf{12},cf{13},cf{14},cf{15},cf{16},cf{17},cf{18}]=CodonFrequency([speciesName{tP},'FFF.fa']);
    
    
    for aa=1:18
        xLength=XX(:,aa);
            for sampleSize=1:1000
                for sqID=1:length(xLength)
                    interestL=xLength(sqID);
                    if ~isnan(interestL)
                        DTa(sqID)=RandSqSn(interestL,cf,aa,2);
                        DTab(sqID)=RandSqSn(interestL,cf,aa,1);
                    else
                        DTa(sqID)=NaN;
                        DTab(sqID)=NaN;
                    end
                end
                seleLumdaA(sampleSize)=calLumda(DTa);
                seleLumdaAb(sampleSize)=calLumda(DTab);
            end
            
            meanA(aa)=mean(seleLumdaA); %%%%% mean of sample mean represent population u
            stdA(aa) = std(seleLumdaA); %%%%% Nomalised by N-1,unbiased estimator
            meanAb(aa)=mean(seleLumdaAb);
            stdAb(aa)=std(seleLumdaAb);
        
    end
    
    fprintf(fileID,fmt2,speciesName{tP},meanA,stdA,meanAb,stdAb);
end

fclose(fileID);


