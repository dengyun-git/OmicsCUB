% function [Pmax,Xmax] = getPmax(cLeng,cfPart,Xpre,rmv)
function Pmax = getPmax(cLeng,cfPart,Xpre,rmv)

pcount=1;
switch cLeng
    case 2
        for i=0:rmv
            j=rmv-i;
            if j>=0
                mnvect=[i,j];    %% mnvect: all the possible rmv (remain values) partitions (order matters here)
                X{pcount}=Xpre+mnvect;
                p(pcount)=mnpdf(X{pcount},cfPart);
                pcount=pcount+1;
            end
        end
        
        
    case 3
        for i=0:rmv
            for j=0:rmv
                k=rmv-i-j;
                if k>=0
                    mnvect=[i,j,k];
                    X{pcount}=Xpre+mnvect;
                    p(pcount)=mnpdf(X{pcount},cfPart);
                    pcount=pcount+1;
                end
            end
        end
        
        
        
    case 4
        for i=0:rmv
            for j=0:rmv
                for k=0:rmv
                    r=rmv-i-j-k;
                    if r>=0
                        mnvect=[i,j,k,r];
                        X{pcount}=Xpre+mnvect;
                        p(pcount)=mnpdf(X{pcount},cfPart);
                        pcount=pcount+1;
                    end
                end
            end
        end
        
        
    case 6
        for i=0:rmv
            for j=0:rmv
                for k=0:rmv
                    for r=0:rmv
                        for s=0:rmv
                            t=rmv-i-j-k-r-s;
                            if t>=0
                                mnvect=[i,j,k,r,s,t];
                                X{pcount}=Xpre+mnvect;
                                p(pcount)=mnpdf(X{pcount},cfPart);
                                pcount=pcount+1;
                            end
                        end
                    end
                end
            end
        end
        
        
end
[Pmax,idMAX]=max(p);
% Xmax=X{idMAX};

end