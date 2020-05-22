%%%%% ARTIFICIAL certain sublengths Lumbda data sets; Currently probe into l=11,30,66

load('speciesName.mat');

interestL=30;

fmt1=[repmat('%s,',1,18),'%s\n'];
fmt2=['%s,',repmat('%d,',1,17),'%d\n'];

fileID=fopen('461LumdaArti30.txt','a'); %%%%I:original;Ir:expected entropy;Ia:equal artificial;Iab: biased artificial.
fprintf(fileID,'%s\n','Lambda sets for sublength: 30 for artificial sequences');
fprintf(fileID,fmt1,'Species','Ealambda','Halambda','Qalambda','Falambda','Yalambda','Calambda','Nalambda','Kalambda','Dalambda','Ialambda','Palambda','Talambda','Aalambda','Valambda','Galambda','Lalambda','Salambda','Ralambda');


for tP=1:length(speciesName)
    seleLumda=zeros(1,18);
    fileName=[speciesName{tP},'SUB.txt'];
    XLL=dlmread(fileName,',',1,18);
    
    for aaCount=1:18
        
        DT = getPiCertainL(aaCount,XLL,interestL);
        
        if (isempty(DT))
            seleLumda(aaCount)=NaN;
        else
            
            [yFit,xpFit]=histcounts(DT,'Normalization','cdf');
            
            xFit=xpFit(2:end);
            fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',0);
            referenceType = fittype('1-exp(a*xFit)','dependent','yFit','independent','xFit','coefficients','a','options',fo);
            yfFit=fit(xFit',yFit',referenceType);
            
            seleLumda(aaCount)=coeffvalues(yfFit);
        end
        
    end
    fprintf(fileID,fmt2,speciesName{tP},seleLumda);
end

fclose(fileID);
