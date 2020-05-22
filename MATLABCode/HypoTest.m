%%% hypothesis test Sn significance
%%% each species
% fileName0='fungiNameListABC.csv'; %%%%for fungi
% fileID0=fopen(fileName0,'r');
% speciesNamep=textscan(fileID0,'%s','Delimiter','\n');
% speciesName=speciesNamep{1,1};
% fclose(fileID0);

speciesName={'saccharomyces_cerevisiae','saccharomyces_arboricola_h_6','saccharomyces_eubayanus','saccharomyces_kudriavzevii_ifo_1802',...
    'aspergillus_calidoustus','aspergillus_fumigatus','aspergillus_lentulus','aspergillus_niger',...
    'fusarium_fujikuroi','fusarium_graminearum','fusarium_oxysporum','fusarium_poae'};

aaList={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
synoL=[repmat(2,1,9),3,repmat(4,1,5),repmat(6,1,3)];
%%% easy reference setting file name
text={'GluEratio','HisHratio','GlnQratio','PheFratio','TyrYratio',...
    'CysCratio','AsnNratio','LysKratio','AspDratio',...
    'IleIratio','ProPratio','ThrTratio','AlaAratio','ValVratio',...
    'GlyGratio','LeuLratio','SerSratio','ArgRratio'};

%%% set format to read from file
fmt1=('%s %u %u %u\n');
fmt2=('%s %u %u %u %u\n');
fmt3=('%s %u %u %u %u %u\n');
fmt4=('%s %u %u %u %u %u %u %u\n');
fmt11=('%s %s %s %s\n');
fmt22=('%s %s %s %s %s\n');
fmt33=('%s %s %s %s %s %s\n');
fmt44=('%s %s %s %s %s %s %s %s\n');
fmList={fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt2,fmt3,fmt3,fmt3,fmt3,fmt3,fmt4,fmt4,fmt4};
fmList2={fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt22,fmt33,fmt33,fmt33,fmt33,fmt33,fmt44,fmt44,fmt44};

fileName2='HypoTest.txt';
fileID2=fopen(fileName2,'a'); %%% write hypo result into
fprintf(fileID2, '%s,%s,%s,%s,%s\n','speciesName{tP}','amino acid','reject H0/not the same distribution','accept H0/the same distribution','difference ratio');


for tP = 1: length(speciesName)
    
    %%% each amino acid
    for CTcount = 1:18
        a=0;
        b=0;
        fileName=[speciesName{tP},text{CTcount},'Fg.txt'];
        fileID = fopen(fileName,'r');
        tt = textscan(fileID,fmList2{CTcount},1,'delimiter',',');
        xx = textscan(fileID,fmList{CTcount},'delimiter',',');
        fclose(fileID);
        
        xxx = cell2mat(xx(:,2:end));
        
        for geneID= 1:size(xxx,1)
          
            [m,X,L] = getXL(xxx(geneID,:),synoL(CTcount));
            
            testVl = sum((X-L/m).^2./(L/m));
            
            p = chi2cdf(testVl,m-1);
            
            if p>0.95
                a = a+1; %%%reject, not the same
            else
                b = b+1; %%%accept, no difference
            end
            
        end
        fprintf(fileID2, '%s,%s,%u,%u,%d\n',speciesName{tP},aaList{CTcount},a,b,a/(a+b));
    end
end

fclose(fileID2);

%%% plot hypotest result in 'plotHopyTest.m'
