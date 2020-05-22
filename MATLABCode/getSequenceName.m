%%read from a file to creat input codon sequences

function SequenceName = getSequenceName(filename)

%%read the sequence from file 'inputSequence'

fileID = fopen(filename,'r');
%tt = textscan(fileID,'%s','delimiter',';');
tt = textscan(fileID,'%s','delimiter','\n'); %%special for Tobias toy genome
fclose(fileID);


%%seperate the sequence name & sequence itself
xx = tt{1};
l = length(xx);
SequenceName = xx(1:2:l,:);

end



%%%%%check with gene ratio file whether they match

% % % fileName='acremonium_chrysogenum_atcc_11550PheFratio.txt';
% % % 
% % % fileID = fopen(fileName);
% % % tt = textscan(fileID,'%s %s %s',1,'delimiter',',');
% % % xx = textscan(fileID,'%s %u %u','delimiter',',');
% % % fclose(fileID);
% % % 
% % % sum(xx{2})/(sum(xx{2})+sum(xx{3}))
