%%%%%%% calculate and write syno ratio for each gene within genome
ctE={'GAG','GAA'};  %%Glu
ctH={'CAT','CAC'};  %%His
ctQ={'CAG','CAA'};  %%Gln
ctF={'TTT','TTC'}; %%Phe
ctY={'TAT','TAC'};   %%Tyr
ctC={'TGT','TGC'};  %%Cys
ctN={'AAT','AAC'};   %%Asn
ctK={'AAG','AAA'};   %%Lys
ctD={'GAT','GAC'};  %%Asp
ctI={'ATA','ATT','ATC'}; %%Ile
ctP={'CCG','CCA','CCT','CCC'};  %%Pro
ctT={'ACG','ACA','ACT','ACC'};  %%Thr
ctA={'GCG','GCA','GCT','GCC'}; %%Ala
ctV={'GTG','GTA','GTT','GTC'}; %%Val
ctG={'GGG','GGA','GGT','GGC'}; %%Gly
ctL={'TTG','TTA','CTG','CTA','CTT','CTC'};  %%Leu
ctS={'AGT','AGC','TCG','TCA','TCT','TCC'};   %%Ser
ctR={'AGG','AGA','CGG','CGA','CGT','CGC'};   %%Arg

text={ctE,ctH,ctQ,ctF,ctY,ctC,ctN,ctK,ctD,ctI,ctP,ctT,ctA,ctV,ctG,ctL,ctS,ctR};

% speciesName={'pyrenochaeta_sp_ds3say3a'};
%
fileName0='bacteriaNameListABC.csv';
% fileName0='protistNameList.csv';
% fileName0='fungiNameListABC.csv'; %%% namelist ABC means alphabet order
fileID0=fopen(fileName0,'r');
speciesNamep=textscan(fileID0,'%s','Delimiter','\n');
speciesName=speciesNamep{1,1};
fclose(fileID0);

for tP=1:length(speciesName)
    %     fileName=[speciesName{tP},'.fa'];%%%fungi
    fileName=[speciesName{tP},'.cds.fa']; %%%%bacteria
    %     fileName=[speciesName{tP},'.cds.all.fa'];  %%%protist file name different
    geneName=getSequenceName(fileName);
    pasteCodon=getCodonSequence(fileName);
    
    fmt1=('%s,%u,%u,%u\n');
    fmt2=('%s,%u,%u,%u,%u\n');
    fmt3=('%s,%u,%u,%u,%u,%u\n');
    fmt4=('%s,%u,%u,%u,%u,%u,%u,%u\n');
    
    fileIDe=fopen([speciesName{tP},'GluEratioBc.txt'],'a');
    fprintf(fileIDe,'%s,%s,%s,%s\n','geneNeme','GAG','GAA','sublength');
    
    fileIDh=fopen([speciesName{tP},'HisHratioBc.txt'],'a');
    fprintf(fileIDh,'%s,%s,%s,%s\n','geneNeme','CAT','CAC','sublength');
    
    fileIDq=fopen([speciesName{tP},'GlnQratioBc.txt'],'a');
    fprintf(fileIDq,'%s,%s,%s,%s\n','geneNeme','CAG','CAA','sublength');
    
    fileIDf=fopen([speciesName{tP},'PheFratioBc.txt'],'a');
    fprintf(fileIDf,'%s,%s,%s,%s\n','geneNeme','TTT','TTC','sublength');
    
    fileIDy=fopen([speciesName{tP},'TyrYratioBc.txt'],'a');
    fprintf(fileIDy,'%s,%s,%s,%s\n','geneNeme','TAT','TAC','sublength');
    
    fileIDc=fopen([speciesName{tP},'CysCratioBc.txt'],'a');
    fprintf(fileIDc,'%s,%s,%s,%s\n','geneNeme','TGT','TGC','sublength');
    
    fileIDn=fopen([speciesName{tP},'AsnNratioBc.txt'],'a');
    fprintf(fileIDn,'%s,%s,%s,%s\n','geneNeme','AAT','AAC','sublength');
    
    fileIDk=fopen([speciesName{tP},'LysKratioBc.txt'],'a');
    fprintf(fileIDk,'%s,%s,%s,%s\n','geneNeme','AAG','AAA','sublength');
    
    fileIDd=fopen([speciesName{tP},'AspDratioBc.txt'],'a');
    fprintf(fileIDd,'%s,%s,%s,%s\n','geneNeme','GAT','GAC','sublength');
    
    fileIDi=fopen([speciesName{tP},'IleIratioBc.txt'],'a');
    fprintf(fileIDi,'%s,%s,%s,%s,%s\n','geneNeme','ATA','ATT','ATC','sublength');
    
    fileIDp=fopen([speciesName{tP},'ProPratioBc.txt'],'a');
    fprintf(fileIDp,'%s,%s,%s,%s,%s,%s\n','geneNeme','CCG','CCA','CCT','CCC','sublength');
    
    fileIDt=fopen([speciesName{tP},'ThrTratioBc.txt'],'a');
    fprintf(fileIDt,'%s,%s,%s,%s,%s,%s\n','geneNeme','ACG','ACA','ACT','ACC','sublength');
    
    fileIDa=fopen([speciesName{tP},'AlaAratioBc.txt'],'a');
    fprintf(fileIDa,'%s,%s,%s,%s,%s,%s\n','geneNeme','GCG','GCA','GCT','GCC','sublength');
    
    fileIDv=fopen([speciesName{tP},'ValVratioBc.txt'],'a');
    fprintf(fileIDv,'%s,%s,%s,%s,%s,%s\n','geneNeme','GTG','GTA','GTT','GTC','sublength');
    
    fileIDg=fopen([speciesName{tP},'GlyGratioBc.txt'],'a');
    fprintf(fileIDg,'%s,%s,%s,%s,%s,%s\n','geneNeme','GGG','GGA','GGT','GGC','sublength');
    
    fileIDl=fopen([speciesName{tP},'LeuLratioBc.txt'],'a');
    fprintf(fileIDl,'%s,%s,%s,%s,%s,%s,%s,%s\n','geneNeme','TTG','TTA','CTG','CTA','CTT','CTC','sublength');
    
    fileIDs=fopen([speciesName{tP},'SerSratioBc.txt'],'a');
    fprintf(fileIDs,'%s,%s,%s,%s,%s,%s,%s,%s\n','geneNeme','AGT','AGC','TCG','TCA','TCT','TCC','sublength');
    
    fileIDr=fopen([speciesName{tP},'ArgRratioBc.txt'],'a');
    fprintf(fileIDr,'%s,%s,%s,%s,%s,%s,%s,%s\n','geneNeme','AGG','AGA','CGG','CGA','CGT','CGC','sublength');
    
    fileList={fileIDe,fileIDh,fileIDq,fileIDf,fileIDy,fileIDc,fileIDn,fileIDk,fileIDd,fileIDi,fileIDp,fileIDt,fileIDa,fileIDv,fileIDg,fileIDl,fileIDs,fileIDr};
    fmList={fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt2,fmt3,fmt3,fmt3,fmt3,fmt3,fmt4,fmt4,fmt4};
    
    for CTcount=1:18
        synoEle=text{CTcount};
        
        for countCd=1:length(pasteCodon)
            
            for synoCount=1:length(synoEle)
                SYNO(synoCount)=length(find(ismember(pasteCodon{countCd},synoEle{synoCount})));
            end
            
            fprintf(fileList{CTcount},fmList{CTcount},geneName{countCd},SYNO,sum(SYNO));
            
        end
        
        fclose(fileList{CTcount});
        
        clear SYNO
        
    end
    
    
end