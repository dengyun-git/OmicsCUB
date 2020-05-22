function [] = artiLumda (interestL)

% load('speciesName.mat');

Ll=10000;

% fileID=fopen(['461Lumda',num2str(interestL),'.txt'],'a'); %%%%I:original;Ir:expected entropy;Ia:equal artificial;Iab: biased artificial.
% fprintf(fileID,'%s\n','Lambda values for sublength: 66, for artificial sequences');
% fileID=fopen('461LambdaAbove59.txt','a');

%%%%%%%%generate aritificial sequences according to codon usage table

% for SPcount=1:461

% i=1;
%     fileCUBName=[speciesName{SPcount},'FFF.fa'];
%     [cfE,cfH,cfQ,cfF,cfY,cfC,cfN,cfK,cfD,cfI,cfP,cfT,cfA,cfV,cfG,cfL,cfS,cfR]=CodonFrequency(fileCUBName);
%
% for aaCount = ['E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R']
%     switch aaCount
%         case 'E'
%             Pmax=EforT(2,interestL);
%             p=cell2mat(cfE);
%         case 'H'
%             Pmax=EforT(2,interestL);
%             p=cell2mat(cfH);
%         case 'Q'
%             Pmax=EforT(2,interestL);
%             p=cell2mat(cfQ);
%         case 'F'
%             Pmax=EforT(2,interestL);
%             p=cell2mat(cfF);
%         case 'Y'
%             Pmax=EforT(2,interestL);
%             p=cell2mat(cfY);
%         case 'C'
%             Pmax=EforT(2,interestL);
%             p=cell2mat(cfC);
%         case 'N'
%             Pmax=EforT(2,interestL);
%             p=cell2mat(cfN);
%         case 'K'
%             Pmax=EforT(2,interestL);
%             p=cell2mat(cfK);
%         case 'D'
%             Pmax=EforT(2,interestL);
%             p=cell2mat(cfD);
%         case 'I'
%             Pmax=EforT(3,interestL);
%             p=cell2mat(cfI);
%         case 'P'
%             Pmax=EforT(4,interestL);
%             p=cell2mat(cfP);
%         case 'T'
%             Pmax=EforT(4,interestL);
%             p=cell2mat(cfT);
%         case 'A'
%             Pmax=EforT(4,interestL);
%             p=cell2mat(cfA);
%         case 'V'
%             Pmax=EforT(4,interestL);
%             p=cell2mat(cfV);
%         case 'G'
%             Pmax=EforT(4,interestL);
%             p=cell2mat(cfG);
%         case 'L'
%             Pmax=EforT(6,interestL);
%             p=cell2mat(cfL);
%         case 'S'
%             Pmax=EforT(6,interestL);
%             p=cell2mat(cfS);
%         case 'R'
%             Pmax=EforT(6,interestL);
%             p=cell2mat(cfR);
%     end
%
%%%%%%%%%%generate artificial sequences according to equal probability
% for aaCount = [2 3 4 6]
for aaCount=6
    
    switch aaCount
        case 2
            Pmax=Efor(2,interestL);
            p=[0.5,0.5];
        case 3
            Pmax=Efor(3,interestL);
            p=[1/3,1/3,1/3];
        case 4
            Pmax=Efor(4,interestL);
            p=[1/4,1/4,1/4,1/4];
        case 6
            Pmax=Efor(6,interestL);
            p=[1/6,1/6,1/6,1/6,1/6,1/6];
    end
    
    
    R= mnrnd(interestL,p,Ll);
    DTp=mnpdf(R,p);
    DT=log(Pmax./DTp);
    
    %%%%count DT frequency
    DTocur = zeros(size(DT)); % vector for freqs
    % frequency for each value
    for ocur = 1:length(DT)
        DTocur(ocur) = sum(DT==DT(ocur));
    end
    
%     DTfinal=exp(DT.*DTocur);
    DTfinal=DT.*DTocur;
    
    [yFit,xpFit]=histcounts(DTfinal,'Normalization','cdf');
    xFit=xpFit(2:end);
    
     modelfun = @(a,xFit)1-exp(a*xFit);
%     modelfun=@(a,xFit)a*xFit;
    beta0=-1;
    
    opts=statset('nlinfit');  %%%%set opts to get a good residual plot representing a good fit
%     opts.RobustWgtFun='bisquare';
%     opts.Robust='on';
    opts.MaxIer=1000;
%     opts.TolFun=exp(-100);
    yfFit=fitnlm(xFit,yFit,modelfun,beta0,'Options',opts);
    
    yfFit.Coefficients.Estimate
    
    
    figure
    hold on
    plot(xFit,yFit,'r.');
    plot(xFit,yfFit.Fitted,'b.')
    legend('orginal','fit')
    hold off
    
    display(['adjusted R-squared is: ',num2str(yfFit.Rsquared.Adjusted),' corresponding to F test p-value: ',num2str(yfFit.Coefficients.pValue)]);
    
    figure
    plot(yfFit.Residuals.Raw,'*');
    title('Residual Plot');
    
    %     i=i+1;
    
end

%i=i-1;

% fprintf(fileID,['%s,',repmat('%d,',1,17),'%d\n'],'Aartificial',repmat(seleLumda(1),1,9),repmat(seleLumda(2),1,1),repmat(seleLumda(3),1,5),repmat(seleLumda(4),1,3));

% end

% fclose(fileID);

end
