%%%%%%% Use Sn and expected Sn, calculate selection pressue, compare with
%%%%%%% artificial sequences
global cfG ctG cfD ctD cfE ctE cfV ctV cfA ctA cfR ctR cfS ctS cfK ctK cfN ctN cfI ctI cfT ctT cfC ctC cfY ctY cfL ctL cfF ctF cfQ ctQ cfH ctH cfP ctP;

load('AveEntropy2f.mat');
load('AveEntropy3f.mat');
load('AveEntropy4f.mat');
load('AveEntropy6f.mat');

SSS=dlmread('scSUBSnxy.txt',',',1,0);  %%%load Sn information
XXX=dlmread('scSublength.txt',',',1,0);%%%load Sub length information

setSynonymousCodonTable(359,2,'synoTableFungi.txt');
cf={cfE,cfH,cfQ,cfF,cfY,cfC,cfN,cfK,cfD,cfI,cfP,cfT,cfA,cfV,cfG,cfL,cfS,cfR};

%text={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
text={'Glu','His','Gln','Phe','Tyr','Cys','Asn','Lys','Asp','Ile','Pro','Thr','Ala','Val','Gly','Leu','Ser','Arg'};
textCt=1;

for aa=1:18
    aaCount=1+4*(aa-1);
    
    xLength=XXX(:,aa);
    ySn=SSS(:,aaCount);
    
    if aa<10
        aveEp=AveEntropy2;
    elseif aa==10
        aveEp=AveEntropy3;
    elseif aa>10 && aa<16
        aveEp=AveEntropy4;
    else
        aveEp=AveEntropy6f;
    end
    
    [Ll,aveReal,xorder]=getRealExpectedSn(xLength,ySn,aa);
    
    aveE=-aveEp(xorder)./xorder;
    
    %%%%%%%% generate artificial ones
    aveEa=getSubLSnMean(aa,xorder,cf,Ll,2);
    
    aveEab=getSubLSnMean(aa,xorder,cf,Ll,1);

    
    figure
    hold on
    plot(xorder,aveE,'b.');
    plot(xorder,aveReal,'r.');
    plot(xorder,aveEa,'y.');
    plot(xorder,aveEab,'g.');
    hold off
    
    xlabel('sublength');
    ylabel('Sn value');
    legend('Reference Expected Sn','Expected Sn Of Real Genome','Expected Sn of Equal Replaced Artificials','Expected Sn of Biased Weights Replaces Artificials');%%% corespond to each sublength
    title(['Expected Sn Comparison for Amino Acid ',text{textCt},' in S.erevisiae']);
    
    textCt=textCt+1;
    
    SelecPressure1(aa)=sum(abs(aveE-aveReal)/length(Ll));
    SelecPressure2(aa)=sum(abs(aveE-aveEa)/length(Ll));
    SelecPressure3(aa)=sum(abs(aveE-aveEab)/length(Ll));
    
end

%%%plot comparison M distance
figure
hold on
plot(SelecPressure1,'r*');
plot(SelecPressure2,'y*');
plot(SelecPressure3,'c*');
hold off
xlabel('amino acid');
ylabel('MD Values');
legend('Original Genome','Equally Replaced Artificial','Biased Weights Replaced Artificial');
title('Codon Usage Bias Comparison in S.cerevisiae')
set(gca,'xTick',1:1:18);
xtxt={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
set(gca,'xTickLabel',text);

