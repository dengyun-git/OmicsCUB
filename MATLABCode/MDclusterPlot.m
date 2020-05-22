list={'saccharomyces_arboricola_h_6','saccharomyces_eubayanus','saccharomyces_cerevisiae_x_saccharomyces_kudriavzevii_vin7',...
    'saccharomyces_kudriavzevii_ifo_1802','saccharomyces_sp_boulardii_',...
    'fusarium_fujikuroi','fusarium_oxysporum','fusarium_verticillioides','fusarium_graminearum','fusarium_poae'};


fmt1=[repmat('%s ',1,18),'%s\n'];
List={'sah6','saeu','sc','saku','sp','ff','fo','fv','fg','fp'};
spCount=10;

fmt2=['%s ',repmat('%d ',1,17),'%d\n'];
fileID2=fopen('461ExpectedComp.txt');
tt2=textscan(fileID2,fmt1,1,'Delimiter',',');
XX2=textscan(fileID2,fmt2,'Delimiter',',');
fclose(fileID2);

filename1='Jeremie_Species_Context.csv';
fileID1=fopen(filename1);
tt1=textscan(fileID1,[repmat('%s ',1,10),'%s\n'],1,'Delimiter',',');
XX1=textscan(fileID1,'%s %s %f %f %f %s %s %s %s %s %s\n','Delimiter',',');
fclose(fileID1);

for i=1:spCount
ind(i)=find(ismember(XX1{1,1},list{i}));
end

xx1=XX1{1,3}(ind);
xx2p=dlmread('461ExpectedComp.txt',',',1,1);%%%%alphabetic order
xx2p(221,11)=0;
xx2p(221,11)=mean(xx2p(:,11));
xx2=xx2p(ind,:);

YY1=pdist(xx1);
ZZ1=linkage(YY1);

YY2=pdist(xx2);
ZZ2=linkage(YY2);

figure
[H1,~,xOrder1]=dendrogram(ZZ1);
set(H1,'LineWidth',4)
xtxt1=List(xOrder1);
set(gca,'xTick',1:1:spCount);
set(gca,'xTickLabel',xtxt1);
xlabel('species');
ylabel('distance');
title('cluster of species based on MD values','fontsize',14);

figure
dendrogram(ZZ2);
[H2,~,xOrder2]=dendrogram(ZZ2);
set(H2,'LineWidth',4)
xtxt2=List(xOrder2);
set(gca,'xTick',1:1:spCount);
set(gca,'xTickLabel',xtxt2);
xlabel('species');
ylabel('distance');
title('cluster of species based on NCBI taxonomy','fontsize',14);

