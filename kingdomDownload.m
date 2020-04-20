%%%%%%download species genome and list species
fileName1='bacteriaNameListW.csv';
fileID1=fopen(fileName1,'a');

fileName2='bacteriaWebList.csv';
fileID2=fopen(fileName2,'a');

fileName3='bacteriaPhylaNumber.csv';
fileID3=fopen(fileName3,'a');

fileName4='bacteriaNameList.csv';
fileID4=fopen(fileName1,'a');

%%%choose from list of wiki:
%%%%%%% Bacteria:
% % % % Actinobacteria % Aquificae % Armatimonadetes % Bacteroidetes % Caldiserica % Chlamydiae % Chlorobi
% % % % Chloroflexi % Chrysiogenetes % Coprothermobacterota % Cyanobacteria % Deferribacteres % Deinococcus-Thermus
% % % % Dictyoglomi % Elusimicrobia % Fibrobacteres % Firmicutes % Fusobacteria % Gemmatimonadetes % Lentisphaerae
% % % % Nitrospirae % Planctomycetes % Proteobacteria % Spirochaetes % Synergistetes % Tenericutes
% % % % Thermodesulfobacteria % Thermotogae % Verrucomicrobia

%%%%%%% Fungi: Fungi (only fungi is ok) or 'Microsporidia', 'Chytridiomycota', 'Blastocladiomycota', 'Neocallimastigomycota', 'Glomeromycota', 'Ascomycota', 'Basidiomycota'

%%%%%%% protist:
% % % % Rhodophyta %Glaucophyta %Stramenopiles %Alveolata %Rhizaria %Euglenozoa %Amoebozoa %Apusozoa %Choanozoa;

bacteriaTax={'Acidobacteria','Acidobacteria','Aquificae','Armatimonadetes','Bacteroidetes','Caldiserica','Chlamydiae','Chlorobi','Chloroflexi',...
    'Chrysiogenetes','Cyanobacteria','Deferribacteres','Deinococcus-Thermus',...
    'Dictyoglomi','Elusimicrobia','Fibrobacteres','Firmicutes','Fusobacteria','Gemmatimonadetes',...
    'Lentisphaerae','Nitrospirae','Planctomycetes','Proteobacteria','Spirochaetes','Synergistetes','Tenericutes','Thermodesulfobacteria',...
    'Thermotogae','Verrucomicrobia'};

for Btax=1:length(bacteriaTax)
    string=bacteriaTax{Btax};
    weblink=['http://rest.ensemblgenomes.org/info/genomes/taxonomy/',string,'?content-type=application/json'];
    option=weboptions('Timeout',120);
    strutPre=webread(weblink);
    [~,uniID]=unique({strutPre.species_taxonomy_id});%%%%eliminate same species different GCA
    strut=strutPre(uniID);
    
    fprintf(fileID3,'%s,%u\n',string,length(strut)); %%%% generally master spicies dispersion beneath each phyla
    
    for stCount=1:length(strut)
        dbnamePre=strut(stCount).dbname;
        dbnameAf=strrep(dbnamePre,'_core_41_94_1',''); %%%% stripped string should be modified according to current ftp address
        webFile=['ftp://ftp.ensemblgenomes.org/pub/bacteria/release-40/fasta/',dbnameAf,'/',strut(stCount).species,'/cds/*.cds.all.fa.gz'];
        
        %%%%special for fungi
        %if contains(strut(stCount).species,"_gca") %%unique Genome Collections Accession (GCA)after Genome Browser agreement, otherwise no '_gca'
        %    speciesSpecial=extractBefore(strut(stCount).species,"_gca");
        %    webFile=['ftp://ftp.ensemblgenomes.org/pub/fungi/release-36/fasta/',dbnameAf,'/',speciesSpecial,'/cds/*.cds.all.fa.gz'];
        %else
        %    webFile=['ftp://ftp.ensemblgenomes.org/pub/fungi/release-36/fasta/',strut(stCount).species,'/cds/*.cds.all.fa.gz'];
        %end
        
        %%%%special for protist
        %webFile=['ftp://ftp.ensemblgenomes.org/pub/protists/release-40/fasta/',dbnameAf,'/',strut(stCount).species,'/cds/*.cds.all.fa.gz'];
        %%%% different release can be adjusted from the above line 'release-*'
        
        if length(strut) < 25  %%%%% control how many species to be selected and avoid aggregation into one group
            fprintf(fileID1,'%s\n',[strut(stCount).species,'.',strut(stCount).assembly_name]);
            fprintf(fileID2,'%s\n',webFile);
            fprintf(fileID4,'%s\n',strut(stCount).species);
        elseif mod(stCount,floor(length(strut)/25))==0
            fprintf(fileID1,'%s\n',[strut(stCount).species,'.',strut(stCount).assembly_name]); %%%convenient for our own file referencing
            %%%%special for fungi%% fprintf(fileID1,'%s\n',speciesSpecial);
            fprintf(fileID2,'%s\n',webFile);
            fprintf(fileID4,'%s\n',strut(stCount).species);
            %%%%special for fungi%% fprintf(fileID1,'%s\n',speciesSpecial);
        end
        
    end
end

% fclose(fileID1);
% fclose(fileID2);
% fclose(fileID3);
% fclose(fileID4);



%%%%%%%%bacteria 442 namelist in taxonomy order bacteriaStill.csv
% % % filename1='bacteriaStill.csv';
% % % fileID1=fopen(filename1,'r');
% % % xx1=textscan(fileID1,'%s','Delimiter','\n');
% % % fclose(fileID1);
% % %
% % %
% % % xx2=sort(xx1{1,1});
% % %
% % % filename3='bacteriaNameListABC.csv';
% % % fileID3=fopen(filename3,'a');
% % %
% % % for i=1:length(xx2)
% % %
% % %     fprintf(fileID3,'%s\n',xx2{i});
% % %
% % % end
% % %
% % % fclose(fileID3);