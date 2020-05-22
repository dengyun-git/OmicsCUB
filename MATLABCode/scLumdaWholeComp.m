%%%%%%among whole genome: sc lambda fit on S, comparison between real genome & artificial ones

SSS=dlmread('scSUBSnxy.txt',',',1,0);  %%%load Sn information
XXX=dlmread('scSublength.txt',',',1,0);%%%load Sub length information

text={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};

cf=cell(1,18);
[cf{1},cf{2},cf{3},cf{4},cf{5},cf{6},cf{7},cf{8},cf{9},cf{10},cf{11},cf{12},cf{13},cf{14},cf{15},cf{16},cf{17},cf{18}]=CodonFrequency('scFF.fa');

fmt1=[repmat('%s,',1,90),'%s\n'];
fmt2=['%s,',repmat('%d,',1,89),'%d\n'];
fileID=fopen('LumdaWholeCompSc.txt','a'); %%%%I:original;Ir:expected entropy;Ia:equal artificial;Iab: biased artificial.
fprintf(fileID,fmt1,'Species','Elambda','Hlambda','Qlambda','Flambda','Ylambda','Clambda','Nlambda','Klambda','Dlambda','Ilambda','Plambda','Tlambda','Alambda','Vlambda','Glambda','Llambda','Slambda','Rlambda','Ealambda','Halambda','Qalambda','Falambda','Yalambda','Calambda','Nalambda','Kalambda','Dalambda','Ialambda','Palambda','Talambda','Aalambda','Valambda','Galambda','Lalambda','Salambda','Ralambda','EaStd','HaStd','QaStd','FaStd','YaStd','CaStd','NaStd','KaStd','DaStd','IaStd','PaStd','TaStd','AaStd','VaStd','GaStd','LaStd','SaStd','RaStd','Eablambda','Hablambda','Qablambda','Fablambda','Yablambda','Cablambda','Nablambda','Kablambda','Dablambda','Iablambda','Pablambda','Tablambda','Aablambda','Vablambda','Gablambda','Lablambda','Sablambda','Rablambda','EabStd','HabStd','QabStd','FabStd','YabStd','CabStd','NabStd','KabStd','DabStd','IabStd','PabStd','TabStd','AabStd','VabStd','GabStd','LabStd','SabStd','RabStd');


for aa=1:18
    aaCount=1+4*(aa-1);
    ySnp=SSS(:,aaCount);
    xLength=XXX(:,aa);
    ySn=ySnp.*xLength;
    
    seleLumda(aa)=calLumda(ySn);
    
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
    
    % %     savFile=['scLumdaComp',num2str(aa),'.mat'];  %%%%seleLumda: observed lambda,seleLumdaA: for equal replaced;seleLumdaAb: for codon table replaced.
    % %     save(savFile,'seleLumda','seleLumdaA','seleLumdaAb');
   
    fprintf(fileID,fmt2,'saccharomyces cerevisiae',seleLumdaSC,meanA,stdA,meanAb,stdAb); 
    
end
fclose(fileID);


%%%%%%%%%lambda distribution of saccharomyces group 
% % % load('speciesName.mat');
% % % nameAlphabet=speciesName;
% % % 
% % % text={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
% % % interestId=find(contains(nameAlphabet,'saccharomyces'));
% % % SSSp=dlmread('461Lumda.txt',',',1,1);
% % % SSS=SSSp(interestId,:);
% % % SSS(12,:)=seleLumdaSC; %%%%seleLumdaSC is 18 lambda for 18 aa in sc
% % % for aa=1:18  %%%%%plot each amino acid; save lambda and standard deviation
% % %     savFile=['scLumdaComp',num2str(aa),'.mat'];  %%%%seleLumda: observed lambda,seleLumdaA: for equal replaced;seleLumdaAb: for codon table replaced.
% % %     load(savFile);
% % %     figure
% % %     hold on  %%%sample size =1000
% % %     histogram(seleLumdaA,'Normalization','probability'); 
% % %     histogram(seleLumdaAb,'Normalization','probability');
% % %     histogram(SSS(:,aa),'BinWidth',0.05,'Normalization','probability');
% % %     plot(seleLumda,0,'r.'); %%%%% lambda for SC
% % %     xlabel('lambda value'); 
% % %     ylabel('probability');
% % %     legend('sc: equal replaced arti lambda distrituion','sc: codon table replaced arti lambda distribution','saccharomyces group lambda','lambda for sc')
% % %     title(['Lambda distribution (among whole genome): ',text{aa}]);
% % %     hold off
% % % end



%%%%%%%%%%%%%%%standard deviation analysis
% % Dt=dlmread('LumdaWholeCompSc.txt',',',1,1);
% % 
% % XXX=dlmread('scSublength.txt',',',1,0);
% % 
% % meanL=mean(XXX,'omitnan');
% % 
% % meanA=Dt(19:36);
% % meanAb=Dt(55:72);
% % stdA=Dt(37:54);
% % stdAb=Dt(73:90);
% % 
% % figure;plot(meanA,'r.');hold on;plot(meanAb,'b.');hold off;
% % title('mean');
% % figure;plot(stdA,'r.');hold on;plot(stdAb,'b.');hold off;
% % title('standard deviation');
% % 
% % figure;plot(meanL,stdA,'r.');hold on;plot(meanL,stdAb,'b.');hold off;
% % title('relationship between length and standard deviation');
