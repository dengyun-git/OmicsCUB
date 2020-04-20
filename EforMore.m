% function [Pmax,X] = EforMore(cLeng,subSleng,cfPart,CaSe) %% pmax for synonymous codons equally likely show up, to assign occurrence as Domonique suggested
function [Pmax,X] = EforMore(cLeng,subSleng,cfPart)
% function Pmax = EforMore(cLeng,subSleng,cfPart,CaSe) %% pmax for synonymous codons equally likely show up, to assign occurrence as Domonique suggested

% switch CaSe
%     
%     case 1
%         
%         i=mod(subSleng,cLeng);  %%subSleng: the length of subsequence
%         
%         X=zeros(1,cLeng);
%         
%         for j=1:(i+1)
%             X(j)=ceil(subSleng/cLeng);
%         end
%         
%         
%         for k=j:cLeng
%             X(k)=floor(subSleng/cLeng);
%         end
%         
%         Pmax=mnpdf(X,cfPart);
        
%     case 2
        %%%%Pmax calculation using codon table frequencies
        P=cfPart;
        X=subSleng.*cfPart;
        Xtest=mod(X,1).*10; %%%%decimal part, *10, if not divisable, will not be 0
        if isempty(find(Xtest,1))
            Pmax=mnpdf(X,P); %% calculate maximum probability
        else
            Xpre=subSleng.*cfPart;
            rmv=subSleng-sum(floor(Xpre));  %%%remainer value
            %             [Pmax,X] = getPmax(cLeng,cfPart,floor(Xpre),rmv);
            Pmax=getPmax(cLeng,cfPart,floor(Xpre),rmv); %%% to save computping time, not care Xmax
        end
% end


end

%%%%%test pi > pmax among the whole file
%%%%fileName1='saccharomyces_cerevisiaeForTab.txt';

%     fileIDtemp=fopen(fileName1,'r');
%     XX0=textscan(fileIDtemp,[repmat('%s ',1,7),'%s\n'],1,'Delimiter',',');
%     XX=textscan(fileIDtemp,'%s %u %f %u %f %u %u %u\n','Delimiter',',');    %%%XX{1}:aa, XX{2}:sublength;XX{3}:Pi;XX{4}:DegPi;XX{5}:Pmax;XX{7}:frequency
%     fclose(fileIDtemp);
%
%     P=XX{3};
%     Pmax=XX{5};
%
%     length(find(P>Pmax))
