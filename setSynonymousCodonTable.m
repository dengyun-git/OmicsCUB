
%%%%%% initialise codon usage ratio by equal probability or codon usage table
function [] = setSynonymousCodonTable(speciesIndex,caseProb,kingdomName)
              %%kingdomName='synoTableFungi.txt'; sc index is 359

global cfG ctG cfD ctD cfE ctE cfV ctV cfA ctA cfR ctR cfS ctS cfK ctK cfN ctN cfI ctI cfT ctT cfC ctC cfY ctY cfL ctL cfF ctF cfQ ctQ cfH ctH cfP ctP;

switch caseProb
    case 1
        % prompt = 'Do you want equal probability for each synonymous codon? If so type y';
        % str=input(prompt,'s');
        %
        % if isempty(str)
        % use *FFF.fa file to calculate
        % species=speciesName{speciesIndex};
        % species='saccharomyces_cerevisiae';
        % GenomeFile=[species,'FFF.fa'];
        % [cfE,cfH,cfQ,cfF,cfY,cfC,cfN,cfK,cfD,cfI,cfP,cfT,cfA,cfV,cfG,cfL,cfS,cfR]=CodonFrequency(GenomeFile);
        % % else
        cfG={[1/4,1/4,1/4,1/4]};
        cfD={[1/2,1/2]};
        cfE={[1/2,1/2]};
        cfV={[1/4,1/4,1/4,1/4]};
        cfA={[1/4,1/4,1/4,1/4]};
        cfR={[1/6,1/6,1/6,1/6,1/6,1/6]};
        cfS={[1/6,1/6,1/6,1/6,1/6,1/6]};
        cfK={[1/2,1/2]};
        cfN={[1/2,1/2]};
        cfM={[1.00]};
        cfI={[1/3,1/3,1/3]};
        cfT={[1/4,1/4,1/4,1/4]};
        cfW={[1.00]};
        cfZ={[1/3,1/3,1/3]};
        cfC={[1/2,1/2]};
        cfY={[1/2,1/2]};
        cfL={[1/6,1/6,1/6,1/6,1/6,1/6]};
        cfF={[1/2,1/2]};
        cfQ={[1/2,1/2]};
        cfH={[1/2,1/2]};
        cfP={[1/4,1/4,1/4,1/4]};
        % end
        
    case 2
        [cfE,cfH,cfQ,cfF,cfY,cfC,cfN,cfK,cfD,cfI,cfP,cfT,cfA,cfV,cfG,cfL,cfS,cfR]=CodonFrequency2(speciesIndex,kingdomName);
end

ctG={'GGG','GGA','GGT','GGC'}; %%Gly
ctD={'GAT','GAC'};  %%Asp
ctE={'GAG','GAA'};  %%Glu
ctV={'GTG','GTA','GTT','GTC'}; %%Val
ctA={'GCG','GCA','GCT','GCC'}; %%Ala
ctR={'AGG','AGA','CGG','CGA','CGT','CGC'};   %%Arg
ctS={'AGT','AGC','TCG','TCA','TCT','TCC'};   %%Ser
ctK={'AAG','AAA'};   %%Lys
ctN={'AAT','AAC'};   %%Asn
% ctM={'ATG'};  %%Met
ctI={'ATA','ATT','ATC'}; %%Ile
ctT={'ACG','ACA','ACT','ACC'};  %%Thr
% ctW={'TGG'};   %%Trp
ctZ={'TGA','TAG','TAA'};  %%End
ctC={'TGT','TGC'};  %%Cys
ctY={'TAT','TAC'};   %%Tyr
ctL={'TTG','TTA','CTG','CTA','CTT','CTC'};  %%Leu
ctF={'TTT','TTC'}; %%Phe
ctQ={'CAG','CAA'};  %%Gln
ctH={'CAT','CAC'};  %%His
ctP={'CCG','CCA','CCT','CCC'};  %%Pro

end


%%%%%Glu ={'GAG','GAA'}; 
% His ={'CAT','CAC'};
% Gln ={'CAG','CAA'};
% Phe ={'TTT','TTC'};
% Tyr ={'TAT','TAC'}; 
% Cys ={'TGT','TGC'};
% Asn ={'AAT','AAC'}; 
% Lys ={'AAG','AAA'}; 
% Asp ={'GAT','GAC'};
% Ile ={'ATA','ATT','ATC'}; 
% Pro ={'CCG','CCA','CCT','CCC'}; 
% Thr ={'ACG','ACA','ACT','ACC'};  
% Ala ={'GCG','GCA','GCT','GCC'}; 
% Val ={'GTG','GTA','GTT','GTC'}; 
% Gly ={'GGG','GGA','GGT','GGC'};
% Leu ={'TTG','TTA','CTG','CTA','CTT','CTC'};  
% Ser ={'AGT','AGC','TCG','TCA','TCT','TCC'};   
% Arg ={'AGG','AGA','CGG','CGA','CGT','CGC'};   