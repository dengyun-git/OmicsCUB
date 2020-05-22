% % % % for interestL=[11 30 66 70 90 120]  %%%can be used as main body
% % % % Sublumbda(interestL);
% % % % artiLumda(interestL);
% % % % clear
% % % % end


function [] = SubLlumbda (interestL)
%%%%%certain sublengths Lumbda data sets; Currently probe into l=11,30,66

load('speciesName.mat');

%%%interestL=50;

fmt1=[repmat('%s,',1,36),'%s\n'];
fmt2=['%s,',repmat('%d,',1,35),'%d\n'];

fileID=fopen(['461Lumda',num2str(interestL),'.txt'],'a'); %%%%I:original;Ir:expected entropy;Ia:equal artificial;Iab: biased artificial.
% % % % % % fprintf(fileID,'%s\n',['Lambda sets for sublength: ',num2str(interestL)]);
fprintf(fileID,fmt1,'Species','Elambda','Hlambda','Qlambda','Flambda','Ylambda','Clambda','Nlambda','Klambda','Dlambda','Ilambda','Plambda','Tlambda','Alambda','Vlambda','Glambda','Llambda','Slambda','Rlambda','Ecount','Hcount','Qcount','Fcount','Ycount','Ccount','Ncount','Kcount','Dcount','Icount','Pcount','Tcount','Acount','Vcount','Gcount','Lcount','Scount','Rcount');


for tP=1:length(speciesName)
    seleLumda=zeros(1,18);
    fileName=[speciesName{tP},'SUB.txt'];
    XX=dlmread(fileName,',',1,0);
    XLL=dlmread(fileName,',',1,18);
    
    for aaCount=1:18
        
        DT=getCertainL(aaCount,XX,XLL,interestL);
        
        if (isempty(DT))
            seleLumda(aaCount)=NaN;
        else
            
            [yFit,xpFit]=histcounts(DT,'Normalization','cdf');
            
            xFit=xpFit(2:end);
            fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',0);
            referenceType = fittype('1-exp(a*xFit)','dependent','yFit','independent','xFit','coefficients','a','options',fo);
            yfFit=fit(xFit',yFit',referenceType);
            
            PointAmount(aaCount)=length(DT);
            seleLumda(aaCount)=coeffvalues(yfFit);
        end
        
    end
    fprintf(fileID,fmt2,speciesName{tP},seleLumda,PointAmount);
end

fclose(fileID);

end