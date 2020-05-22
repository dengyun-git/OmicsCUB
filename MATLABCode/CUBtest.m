% % % %%%%%check synoTable462v2 comply with gene ratio files
% % % 
% % % synoN=[2,2,2,2,2,2,2,2,2,3,4,4,4,4,4,6,6,6];
% % % 
% % % text={'GluE','HisH','GlnQ','PheF','TyrY','CysC','AsnN','LysK','AspD','IleI','ProP','ThrT','AlaA','ValV','GlyG','LeuL','SerS','ArgR'};
% % % 
% % % speciesName={'pichia_membranifaciens_nrrl_y_2026'};
% % % 
% % % fmt1=('%s %u %u\n');
% % % fmt2=('%s %u %u %u\n');
% % % fmt3=('%s %u %u %u %u\n');
% % % fmt4=('%s %u %u %u %u %u %u\n');
% % % 
% % % fmt11=('%s %s %s\n');
% % % fmt22=('%s %s %s %s\n');
% % % fmt33=('%s %s %s %s %s\n');
% % % fmt44=('%s %s %s %s %s %s %s\n');
% % % 
% % % fmList={fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt2,fmt3,fmt3,fmt3,fmt3,fmt3,fmt4,fmt4,fmt4};
% % % fmList2={fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt22,fmt33,fmt33,fmt33,fmt33,fmt33,fmt44,fmt44,fmt44};
% % % 
% % % for tP=1:length(speciesName)
% % %     
% % %     for CTcount=1:18
% % %         
% % %         fileName2=[speciesName{tP},text{CTcount},'ratio.txt'];
% % %         fileID2 = fopen(fileName2);
% % %         tt = textscan(fileID2,fmList2{CTcount},1,'delimiter',',');
% % %         xx = textscan(fileID2,fmList{CTcount},'delimiter',',');
% % %         fclose(fileID2);
% % %         
% % %         for xxcount=1:synoN(CTcount) 
% % %         XX(:,xxcount)=xx{1,xxcount+1};
% % %         end
% % %         
% % %         cfSum=sum(XX);
% % %         cfWholep{CTcount}=cfSum./sum(cfSum);
% % %         
% % %         clear XX
% % %         
% % %     end
% % %     
% % %     cfWhole=cell2mat(cfWholep);
% % % end
% % % 
% % % 
% % % 
% % % 
% % % %%%%%%check each gene has the right Pi
% % % 
% % % global cfE cfI cfG cfR ctE ctI ctG ctR
% % % 
% % % ctG={'GGG','GGA','GGT','GGC'}; %%Gly
% % % ctE={'GAG','GAA'};  %%Glu
% % % ctR={'AGG','AGA','CGG','CGA','CGT','CGC'};   %%Arg
% % % ctI={'ATA','ATT','ATC'}; %%Ile
% % % 
% % % 
% % % geneId=3537;
% % % 
% % % pasteCodon=getCodonSequence('pichia_membranifaciens_nrrl_y_2026.fa');
% % % 
% % % codon=pasteCodon{1,geneId};
% % % 
% % % cfE={[0.46609786900890915,0.53390213099109085]};
% % % 
% % % cfI={[0.21187676574141309,0.38494535735860341,0.40317787689998352]};
% % % 
% % % cfG={[0.14534999553445682,0.22919228628940841,0.33964234434146978,0.28581537383466499]};
% % % 
% % % cfR={[0.16663685228250408,0.44293952791052277,0.10995544879166561,0.07216784646443995,0.12123380440060311,0.08706652015026449]};
% % % 
% % % [~,~,Pie]=GluAminoAcidH(codon);
% % % 
% % % [~,~,Pii]=IleAminoAcidH(codon);
% % % 
% % % [~,~,Pig]=GlyAminoAcidH(codon);
% % % 
% % % [~,~,Pir]=ArgAminoAcidH(codon);
% % % 
% % % 
% % 
% % %%%%%% find syno ratio among certain sublengths
% % synoN=[2,2,2,2,2,2,2,2,2,3,4,4,4,4,4,6,6,6];
% % 
% % text={'GluE','HisH','GlnQ','PheF','TyrY','CysC','AsnN','LysK','AspD','IleI','ProP','ThrT','AlaA','ValV','GlyG','LeuL','SerS','ArgR'};
% % 
% % speciesName={'rosellinia_necatrix'};
% % 
% % fmt1=('%s %u %u\n');
% % fmt2=('%s %u %u %u\n');
% % fmt3=('%s %u %u %u %u\n');
% % fmt4=('%s %u %u %u %u %u %u\n');
% % 
% % fmt11=('%s %s %s\n');
% % fmt22=('%s %s %s %s\n');
% % fmt33=('%s %s %s %s %s\n');
% % fmt44=('%s %s %s %s %s %s %s\n');
% % 
% % fmList={fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt2,fmt3,fmt3,fmt3,fmt3,fmt3,fmt4,fmt4,fmt4};
% % fmList2={fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt22,fmt33,fmt33,fmt33,fmt33,fmt33,fmt44,fmt44,fmt44};
% % 
% % 
% % for tP=1:length(speciesName)
% %     
% %     for CTcount=1:1
% %         
% %         fileName2=[speciesName{tP},text{CTcount},'ratio.txt'];
% %         fileID2 = fopen(fileName2);
% %         tt = textscan(fileID2,fmList2{CTcount},1,'delimiter',',');
% %         xx = textscan(fileID2,fmList{CTcount},'delimiter',',');
% %         fclose(fileID2);
% %         
% %         lenW=xx{1,2}+xx{1,3};
% %         
% %         idL1=find(lenW>=5 & lenW<=10);
% %         idL2=find(lenW>10);
% %         
% %         cfSum11=xx{1,2}(idL1,:);
% %         cfSum12=xx{1,3}(idL1,:);
% %         cf5b10(1)=sum(cfSum11)./(sum(cfSum11)+sum(cfSum12)); %%%%%length between 5 to 10
% %         cf5b10(2)=sum(cfSum12)./(sum(cfSum11)+sum(cfSum12));
% %         
% %         cfSum21=xx{1,2}(idL2,:);
% %         cfSum22=xx{1,3}(idL2,:);
% %         cfa10(1)=sum(cfSum21)./(sum(cfSum21)+sum(cfSum22));  %%%% length above 10
% %         cfa10(2)=sum(cfSum22)./(sum(cfSum21)+sum(cfSum22));        
% %     end
% % end
