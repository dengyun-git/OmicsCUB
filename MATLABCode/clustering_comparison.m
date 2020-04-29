function [ClusterSam,BoundUpLow] = clustering_comparison(ZZ1,ZZ2,kClust)%%%%ZZ1,ZZ2 arbitary index for individual clusters

M=zeros(kClust,kClust);
for M1count=1:kClust
    for M2count=1:kClust
        M(M1count,M2count)=length(find((ismember(ZZ1{M1count},ZZ2{M2count})))); %%%%evluate matrix: entry contain: counts of common elements in two branches in ZZ1 and ZZ2
    end
end

Mi=sum(M,2);
Mj=sum(M,1);
n=sum(M(:));

T=sum(M(:).^2)-n;

P=sum(Mi(:).^2)-n;

Q=sum(Mj(:).^2)-n;

ClusterSam= T/sqrt(P*Q);

PP=sum(Mi(:).*(Mi(:)-1).*(Mi(:)-2));

QQ=sum(Mj(:).*(Mj(:)-1).*(Mj(:)-2));

meanB=sqrt(P*Q)/(n*(n-1)); %%%% the mean and std: random unrelated cluster trees.

varB=2/(n*(n-1))+4*PP*QQ/(n*(n-1)*(n-2)*P*Q)+(P-2-4*PP/P)*(Q-2-4*QQ/Q)/(n*(n-1)*(n-2)*(n-3))-P*Q/(n.^2*(n-1).^2);

BoundUpLow=[meanB-2*sqrt(varB),meanB+2*sqrt(varB)];%%%% when Bk locate outside this area, means similaries is significance

end
