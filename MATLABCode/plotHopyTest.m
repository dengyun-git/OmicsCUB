fileName='HypoTest.txt';
fileID=fopen(fileName,'r'); %%% write hypo result into
tt=textscan(fileID, '%s %s %s %s %s',1,'Delimiter',',');
xx=textscan(fileID, '%s %s %u %u %f','Delimiter',',');
fclose(fileID);

X=reshape(xx{1,5},18,12);

speciesName={'saccharomyces_cerevisiae','saccharomyces_arboricola_h_6','saccharomyces_eubayanus','saccharomyces_kudriavzevii_ifo_1802',...
    'aspergillus_calidoustus','aspergillus_fumigatus','aspergillus_lentulus','aspergillus_niger',...
    'fusarium_fujikuroi','fusarium_graminearum','fusarium_oxysporum','fusarium_poae'};
aaList={'E','H','Q','F','Y','C','N','K','D','I','P','T','A','V','G','L','S','R'};
speciesNameshort={'sc','sah6','saeu','saku','ac','af','al','an','ff','fg','fo','fp'};

h = heatmap(speciesNameshort,aaList,X);

h.Title = 'Hypothesis Test for Sn';
h.XLabel = 'Species';
h.YLabel = 'Amino Acids';