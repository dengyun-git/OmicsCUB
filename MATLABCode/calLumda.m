%%%main body: lambda fit goodness
% aaN=11; %%%choose which amino acid to investigate
% load('orderName461.mat')
% species1=speciesName{100};
% SSS1=dlmread([species1,'SUB.txt'],',',1,0); 
% seleLumda=calLumda(aaN,SSS1(:,aaN))
%%%%% call the following function

function seleLumda = calLumda(aa,ySn)

text={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};

% [yFit,xpFit]=histcounts(ySn,'Normalization','cdf');
% xFit=xpFit(2:end);
[yFit,xpFit]=histcounts(ySn,'Normalization','probability');
xFit=xpFit(2:end);
phat=gamfit(xFit',yFit');
yfFit=gampdf(xFit',phat(1),phat(2));
resD=yFit'-yfFit;
plot(resD,'*');
% modelfun = @(a,x)1-exp(a*xFit);
% beta0=-0.5;
% opts=statset('nlinfit');  %%%%set opts to get a good residual plot representing a good fit
% opts.RobustWgtFun='bisquare';
% opts.Robust='on';
% opts.MaxIer=10000;
% opts.TolFun=exp(-100);
% yfFit=fitnlm(xFit',yFit',modelfun,beta0,'Options',opts);
% display(['amino acid: ',text{aa},' adjusted R-squared is: ',num2str(yfFit.Rsquared.Adjusted),' corresponding to F test p-value: ',num2str(yfFit.Coefficients.pValue)]);
% figure
% % histogram(yfFit.Residuals.Raw,'BinWidth',0.005);
% plot(yfFit.Residuals.Raw,'*');
xlabel('residuals');
ylabel('difference between original data and fit');
title(['amino acid: ',text{aa},' lambda fit residual plot']);
% seleLumda=yfFit.Coefficients.Estimate;

figure
hold on
plot(xFit,yFit,'r.');
plot(xFit,yfFit.Fitted,'b');
title(['lambda: (regression among whole sc genome)= ',num2str(seleLumda)]);
xlabel('S value');
ylabel('cumulative distribution');
legend('observed sequence','fit curve');
hold off
end