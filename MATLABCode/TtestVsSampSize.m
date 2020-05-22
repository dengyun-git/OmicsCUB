%%%test Temperature theory sample size influence on slope
% aaList={'E','I','G','R'};
% synoL=[2,3,4,6];

aaList={'E'};
synoL=2;

for aa=1:length(aaList)
    
    for sampSize=10
        
        cLeng=synoL(aa);
        Psyno=zeros(1,cLeng);%%%%% prepare P for multinomial generator
        Psyno(:,:)=1/cLeng;
        
        XLLppre=dlmread('scSublength.txt',',',1,0);%%%%%species for sc
        XLLpre=unique(XLLppre(:,aa)); %%%%% 461 species length column
        XLL=XLLpre(~isnan(XLLpre));
        
        interestL=3;
        
        Lid=getCertainLInd(aa,interestL,XLLppre);%%%%find certain length genes
        
        Xg=mnrnd(interestL,Psyno,sampSize*length(Lid)); %%%% generate random configurations for Lid counts
        Pi=mnpdf(Xg,Psyno);
        
%         Pmax = EforMore(cLeng,interestL,1,1);
        
        for i=1:sampSize*length(Lid)
            DEG1=getDegereracy(Xg(i,:),cLeng);
            Freq(i)=getSameSnCount(i,Pi);
            Svalue(i)=Pi(i).*DEG1;
        end
        
        freq=Freq/sum(Freq);
        
        
        lmF=fitlm(Svalue,freq);
        
        Slop=lmF.Coefficients{2,1}
        
    end
    
    
%     figure
%     hold on
%     plot((10:10:100),Slop,'.');
%     xlabel('sample size');
%     ylabel('slope');
%     title(['amino acid: ',aaList(aa),' in species sc']);
%     
%     hold off
    
end