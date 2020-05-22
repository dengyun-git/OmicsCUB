%%%%% on myrtle MdIDlist.csv have 344 species(except sc), which aims to match high/low
%%%%% 1450 species list. But *MdIDlist.txt in total have 345 species.

global cfG ctG cfD ctD cfE ctE cfV ctV cfA ctA cfR ctR cfS ctS cfK ctK cfN ctN cfI ctI cfT ctT cfC ctC cfY ctY cfL ctL cfF ctF cfQ ctQ cfH ctH cfP ctP;

fileID00=fopen('MdGlist.csv');
tt=textscan(fileID00,'%s','Delimiter',',');
homoList=tt{1,1};
fclose(fileID00);

%%%%if want to guarantee homologies are selected from the same specie list,
%%%% just change fileName0
fileName0='fungiNameListABC.csv';   %%% for fungi list
fileID0=fopen(fileName0,'r');
speciesNameWp=textscan(fileID0,'%s','Delimiter','\n');
speciesNameW=speciesNameWp{1,1};
fclose(fileID0);


for homoc=1:length(homoList) %%homology counter
    
    stringName=homoList{homoc};
    
    [homoSpecies,homoGene,l] = getHomoInfor(stringName);%%%call function 'getHomoInfor' to retrieve homology information: species+gene
    
    if l~=0
        for jcount=1:l
            speciesName=homoSpecies{jcount};
            if ismember(speciesName,speciesNameW)
                
                geneNamesWhole=getHomologyList([speciesName,'Gene.fa']); %% list all the gene names in the interesed species
                geneName=homoGene{jcount};
                
                id=find(contains(geneNamesWhole,geneName),1);
                if ~isempty(id) %%%%%%%%% some geneName not match geneName list
                    
                    fileID=fopen([speciesName,'MdIDlist.csv'],'a');
                    fprintf(fileID,'%u,',id);
                    fclose(fileID);
                    
                end
            end
        end
    end
    
end
