%%%%%%% for 461 Use Sn and expected Sn, calculate selection pressue, compare with
%%%%%%% artificial sequences
%%% orderName461.mat contain species list in phylogenetic order
%%% '461ExpectedComp.txt' have different list sort to FungiNameListABC.txt.
%%% about '_c...' 3 species
load('AveEntropy2f.mat');
load('AveEntropy3f.mat');
load('AveEntropy4f.mat');
load('AveEntropy6f.mat');
load('speciesName.mat');

fmt1=[repmat('%s,',1,18),'%s\n'];
fmt2=['%s,',repmat('%d,',1,17),'%d\n'];
fileID=fopen('461ExpectedComp.txt','a');
fprintf(fileID,fmt1,'speciesName','E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R');


for tP=1:length(speciesName)
    fileName=[speciesName{tP},'SUB.txt'];
    
    XXX=dlmread(fileName,',',1,0); %%%load sublength information
    %%%SSS=dlmread(fileName,',',1,19);  %%%load Sn information
    
    % fprintf(fileID,fmt1,'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R');
    
    for aa=1:18
        aaCount=aa+18;
        
        xLength=XXX(:,aaCount);
        ySn=XXX(:,aa);
        
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
% % % %         figure
% % % %         hold on
% % % %         plot(xorder,aveE,'r.')
% % % %         plot(xorder,aveReal,'g.')
% % % %         hold off
% % % %         xlabel('sublength');
% % % %         ylabel('Sn value');
% % % %         legend('reference expected Sn','real observed expected Sn');%%% corespond to each sublength
% % % %         title('expected Sn comparison for amino acid G in Species:sc')
% % % %         
        SelecPressure(aa)=sum(abs(aveE-aveReal))/length(Ll); %%%% sum difference among all sublengths
        
    end
    
    fprintf(fileID,fmt2,speciesName{tP},SelecPressure);
end

fclose(fileID);

%%%%%when have scExpectedComp.txt, make plot
% % % XX=dlmread('scExpectedComp.txt',',',1,0);
% % % figure
% % % hold on
% % % plot(XX(1,:),'r.');
% % % plot(XX(2,:),'g.');
% % % plot(XX(3,:),'b.');
% % % hold off
% % %
% % % xlabel('amino acid');
% % % ylabel('Selection Pressure');
% % % legend('observed sequence','equal replaced sequence','codon table replaced sequence');
% % % title('selection pressure comparison')

% % % set(gca,'xTick',1:1:18);
% % % xtxt={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
% % % set(gca,'xTickLabel',xtxt);


%%%%%%%%461 heatmap
heatmapdata=dlmread('461ExpectedComp.txt',',',1,1);
heatmapdata(221,11)=0;
heatmapdata(221,11)=mean(heatmapdata(221,:));%%%tackle with NaN
aaList={'E (2)','H (2)','Q (2)','F (2)','Y (2)','C (2)','N (2)','K (2)','D (2)','I (3)','P (4)','T (4)','A (4)','V (4)','G (4)','L (6)','S (6)','R (6)'};
HM=HeatMap(heatmapdata,'Colormap',[1 1 1;1 0.93,0.93;1 0.88 0.88;1 0.65 0.65;1 0.5 0.5; 1 0.2 0.2; 1 0 0],'ColumnLabels',aaList);
title='Codon Usage Bias in 462 species of Fungi kingdom (based on MD values)';
addTitle(HM,title);
ylabel='species';
addYLabel(HM,ylabel,'Fontsize',15);
%%%%%kmeans method to cluster species for each amino acid
% for hti=[1 10 15 18]
% idx = kmeans(heatmapdata(:,hti),14);
% speciesOrder=(1:461);
% figure
% plot(speciesOrder,idx,'r.');
% xlabel('461 species');
% ylabel()
% 
% end

%%%%kmeans to cluster amino acids
% idx = kmeans(heatmapdata',4);
% aaOrder=(1:18);
% figure
% plot(aaOrder,idx,'r*')
% 
% set(gca,'xTick',1:1:18);
% xtxt={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
% set(gca,'xTickLabel',xtxt);


%%%%%hyrachical cluster for 461 in 18 aa
heatmapdata(7,11)=0;%%%%tacle with NaN VALUE
heatmapdata(7,11)=mean(heatmapdata(:,11));
CGobj=clustergram(heatmapdata,'Colormap',[1 0.85 0.85;1 0.8,0.8;1 0.75 0.75;1 0.6 0.6;1 0.5 0.5;1 0 0],'ColumnLabels',aaList);
addXLabel(CGobj,'amino acid');
addYLabel(CGobj,'461 species--Phylogenetic tree order');
title='cluster tree based on SP3';
addTitle(CGobj,title);

heatmapdataPre(221,11)=0;
heatmapdataPre(221,11)=mean(heatmapdataPre(:,11));
CGobj2=clustergram(heatmapdataPre,'Colormap',[1 0.85 0.85;1 0.8,0.8;1 0.75 0.75;1 0.6 0.6;1 0.5 0.5;1 0 0],'ColumnLabels',aaList);
addXLabel(CGobj2,'amino acid');
addYLabel(CGobj2,'461 species--Alphabetic order');
title='cluster tree based on SP3';
addTitle(CGobj,title);


%%%%%%%%%%%% self organizing map for 461 species
% % % % fmt1=[repmat('%s ',1,18),'%s\n'];
% % % % fmt2=['%s ',repmat('%f ',1,17),'%f\n'];
% % % % 
% % % % fileID=fopen('461ExpectedComp.txt','r');
% % % % X0=textscan(fileID,fmt1,1,'Delimiter',',');
% % % % X1=textscan(fileID,fmt2,'Delimiter',',');
% % % % speciesName=X1{1};
% % % % s1=find(ismember(speciesName,'saccharomyces_arboricola_h_6'));
% % % % s2=find(ismember(speciesName,'saccharomyces_eubayanus'));
% % % % s3=find(ismember(speciesName,'saccharomyces_kudriavzevii_ifo_1802'));
% % % % s4=find(ismember(speciesName,'ashbya_gossypii'));
% % % % s5=find(ismember(speciesName,'yarrowia_lipolytica'));
% % % % s6=find(ismember(speciesName,'schizosaccharomyces_pombe'));
% % % % s7=find(ismember(speciesName,'schizosaccharomyces_japonicus'));
% % % % s8=find(ismember(speciesName,'schizosaccharomyces_octosporus'));
% % % % s9=find(ismember(speciesName,'aspergillus_clavatus'));
% % % % s10=find(ismember(speciesName,'aspergillus_flavus'));
% % % % s11=find(ismember(speciesName,'aspergillus_nidulans'));
% % % % s12=find(ismember(speciesName,'aspergillus_niger'));
% % % % s13=find(ismember(speciesName,'aspergillus_oryzae'));
% % % % s14=find(ismember(speciesName,'aspergillus_terreus'));
% % % % s15=find(ismember(speciesName,'aspergillus_fumigatus'));
% % % % s16=find(ismember(speciesName,'fusarium_fujikuroi'));
% % % % s17=find(ismember(speciesName,'fusarium_graminearum'));
% % % % s18=find(ismember(speciesName,'fusarium_oxysporum'));
% % % % s19=find(ismember(speciesName,'fusarium_verticillioides'));
% % % % s20=find(ismember(speciesName,'trichoderma_reesei'));
% % % % 
% % % % 
% % % % XX=cell2mat(X1(2:19));
% % % % [a,b]=find(isnan(XX));
% % % % XX(a,b)=mean(XX(~isnan(XX(:,b))));
% % % % 
% % % % % for j=1:18    %%%normalisation not necessary
% % % % %     for i=1:461
% % % % %     XXf(i,j)=XX(i,j)/sum(XX(:,j));
% % % % %     end
% % % % % end
% % % % 
% % % % XXf=XX([s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 s14 s15 s16 s17 s18 s19 s20],:);
% % % % net = selforgmap([6 6]);
% % % % net.trainParam.epochs=1000;
% % % % XXff=XXf';
% % % % XXfff=XXff(:,1:5);
% % % % net = train(net,XXfff);
% % % % view(net)
% % % % y = net(XXf);
% % % % classes = vec2ind(y');
% % % % 
% % % % % plotsomhits(net,XXf); %%%%plot sample hits
% % % % % plotsomplanes(net); %%%%plot samplehts from ith input to 2 dimensional map
% % % % % plotsomnc(net);     %%%%Plot self-organizing map neighbor connections
% % % % % plotsomnd(net);  %%%%%%%Plot self-organizing map neighbor distances
% % % % cbh=colorbar('v');;set(cbh,'XTickLabel',{'1','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'});

