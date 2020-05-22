%%%% apply exhaustive data to test parallel lines
%%%% not prediction for artificial data or observed data
NN=25;
syno=3;

pref=getPref(NN,syno);
pmax=Efor(syno,NN);

S1pre=pmax./pref; %%%S1: original
S1=log(S1pre(S1pre>0));

[S2,uniID,~]=unique(S1); %%%%S2: no overlap

for i=1:length(uniID)
    Freq1(i)=sum(S1==S2(i));
end

Freq=log(Freq1/sum(Freq1));
aveFreq=sum(Freq1)/length(uniID);

lmF=fitlm(S2,Freq);

figure
plot(S2,Freq,'.')
hold on
plot(S2,lmF.Fitted,'-')


title(['line fit for length ',num2str(NN),' of syno ',num2str(syno)]);
xlabel('S=log(pmax/pi) value');
ylabel('frequency');

sprintf('%d %d %d %d %d',lmF.Coefficients{2,1},lmF.Coefficients{1,1},lmF.Rsquared.Adjusted,lmF.Rsquared.Adjusted,aveFreq)



