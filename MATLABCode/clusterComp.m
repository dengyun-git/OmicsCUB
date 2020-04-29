load('orderName461.mat')
species1=speciesName{309};
species2=speciesName{307};   %%%%%sah6:242;saeu:243;fg:14;ff:18;sp307;sj:309;

SSS1=dlmread([species1,'SUB.txt'],',',1,0); 
SSS2=dlmread([species2,'SUB.txt'],',',1,0);

xtxtp={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};

for aa=1:18
    XX1p(:,aa)=SSS1(:,aa);
    XX2p(:,aa)=SSS2(:,aa);
end

ZZ1=getDendroZZ(XX1p);
% figure
% [~,~,xOrder1]=dendrogram(ZZ1);
% xtxt1=xtxtp(xOrder1);
% set(gca,'xTick',1:1:18);
% set(gca,'xTickLabel',xtxt1);
% xlabel('amino acid');
% ylabel('distance');
% title(['cluster at amino acid level in: ',species1]);

ZZ2=getDendroZZ(XX2p);
% [~,~,xOrder2]=dendrogram(ZZ2);
% xtxt2=xtxtp(xOrder2);
% set(gca,'xTick',1:1:18);
% set(gca,'xTickLabel',xtxt2);
% xlabel('amino acid');
% ylabel('distance');
% title(['cluster at amino acid level in: ',species2]);

%%%%%%%compare similarity
kClust=6;
ClusterSn1=cluster(ZZ1,'maxclust',kClust);
ClusterSn2=cluster(ZZ2,'maxclust',kClust);
Cbase=(1:18); %%%%%%%original amino acids order
for kClustCount=1:kClust
GrpClust1{kClustCount}=Cbase(ClusterSn1==kClustCount); %%%%% put amino acids into group
GrpClust2{kClustCount}=Cbase(ClusterSn2==kClustCount);
end
[ClusterSam,BoundUpLow] = clustering_comparison(GrpClust1,GrpClust2,kClust)