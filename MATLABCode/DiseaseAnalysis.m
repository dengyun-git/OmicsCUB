answer1 = questdlg('welcome to Entropy Based Analysis For Disease Related Species and Genes App!Have specific species or gene name?', ...
    'Entropy Based Diseased Related CUB App','yes','no','no');

switch answer1
    case 'yes'
        answer11 = questdlg('please choose species or gene under investigation',...
            'search method','species','gene','gene');
        if strcmp(answer11,'species')
            x=1;
        elseif strcmp(answer11,'gene')
            x=2;
        
        end
        
    case 'no'
        answer12 = questdlg('please select which part',...
            'disease part','head','lung','lung');
        if strcmp(answer12,'head')
            x=3;
        elseif strcmp(answer12,'lung')
            x=4;
        end
       
        
end

