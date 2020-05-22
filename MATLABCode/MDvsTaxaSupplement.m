
%%%%%now we try far distance species

spCount=10;

list=cell(spCount,1); %%% how many species into investigation
ind=zeros(spCount,1); %%%% id to select species information from file
% 
[ind1,ind2]=sort(XX1{1,3});  %%% taxonomy id = id1, index of species = id2

for i=1:5
list{i}=XX1{1,1}{ind2(i)}
ind(i)=ind2(i);
list{5+i}=XX1{1,1}{ind2(end-i+1)}
ind(5+i)=ind2(end-i+1);
end

%ind=[1,100,200,300,400,50,150,250,350,450];
%list={'ag','cj','pk','td','zr','lssc','pa','esc','psm','psn'};
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
dendrogram(ZZ1);
[~,~,xOrder1]=dendrogram(ZZ1);
xtxt1=list(xOrder1);
set(gca,'xTick',1:1:spCount);
set(gca,'xTickLabel',xtxt1);
xlabel('species');
ylabel('distance');
title('cluster of species based on MD values');

figure
dendrogram(ZZ2);
[~,~,xOrder2]=dendrogram(ZZ2);
xtxt2=list(xOrder2);
set(gca,'xTick',1:1:spCount);
set(gca,'xTickLabel',xtxt2);
xlabel('species');
ylabel('distance');
title('cluster of species based on NCBI taxonomy');

specialCount=1;

for clustCount=2:5 %%%comparison based on how many clusters
    
kClust=clustCount;
ClusterSn1=cluster(ZZ1,'maxclust',kClust);
ClusterSn2=cluster(ZZ2,'maxclust',kClust);

Cbase=(1:spCount); 
for kClustCount=1:kClust
GrpClust1{kClustCount}=Cbase(ClusterSn1==kClustCount); 
GrpClust2{kClustCount}=Cbase(ClusterSn2==kClustCount);
end
[ClusterSam,BoundUpLow] = clustering_comparison(GrpClust1,GrpClust2,kClust);


    KC(kClust)=kClust;
    CS(kClust)=ClusterSam;
    BDL(kClust)=BoundUpLow(1);
    BDU(kClust)=BoundUpLow(2);
 
 if ClusterSam<BoundUpLow(1) || ClusterSam>BoundUpLow(2)
   specialP(specialCount)=kClust;
   kClust
   CS(kClust)
   BDL(kClust)
   BDU(kClust)
   specialCount=specialCount+1;
 end
 
clear GrpClust1 GrpClust2
end

