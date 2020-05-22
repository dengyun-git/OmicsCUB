%%%%%%using S L datasets '*SUBarti.txt', generate Lambda table using artificial sequences 

load('speciesName.mat');

fmt1=[repmat('%s,',1,18),'%s\n'];
fmt2=['%s,',repmat('%d,',1,17),'%d\n'];

fileID=fopen('461ArtiLumda.txt','a'); %%%%I:original;Ir:expected entropy;Ia:equal artificial;Iab: biased artificial.
fprintf(fileID,fmt1,'Species','Elambda','Hlambda','Qlambda','Flambda','Ylambda','Clambda','Nlambda','Klambda','Dlambda','Ilambda','Plambda','Tlambda','Alambda','Vlambda','Glambda','Llambda','Slambda','Rlambda');


for tP=1:length(speciesName)
    seleLumda=zeros(1,18);
    fileName=[speciesName{tP},'SUBarti.txt'];
    XX=dlmread(fileName,',',1,0);
    
    for aaCount=1:18
        [yFit,xpFit]=histcounts(XX(:,aaCount),'Normalization','cdf');
%         [yFit,xpFit]=histcounts(XX(:,aaCount));

        xFit=xpFit(2:end);
        fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',0);
        referenceType = fittype('1-exp(a*xFit)','dependent','yFit','independent','xFit','coefficients','a','options',fo);
        yfFit=fit(xFit',yFit',referenceType);
        
        seleLumda(aaCount)=coeffvalues(yfFit);
        
%         plot(xFit,yFit,'r');
%         hold on
%         plot(xFit,yfFit,'g');
%         hold off
%         xlabel('S values');
%         ylabel('probability');
%         legend('original','fit');
%         title(['fit for species:',speciesNampe{tp}]);
    end
    fprintf(fileID,fmt2,speciesName{tP},seleLumda);
end

fclose(fileID);
