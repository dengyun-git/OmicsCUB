function [RefCAI,ObsCAI] = getRefCAI(speciesName)

fmt1=('%s %f %f\n');
fmt2=('%s %f %f %f\n');
fmt3=('%s %f %f %f %f\n');
fmt4=('%s %f %f %f %f %f %f\n');

fmt11=('%s %s %s\n');
fmt22=('%s %s %s %s\n');
fmt33=('%s %s %s %s %s\n');
fmt44=('%s %s %s %s %s %s %s\n');

fmList={fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt2,fmt3,fmt3,fmt3,fmt3,fmt3,fmt4,fmt4,fmt4};
fmList2={fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt11,fmt22,fmt33,fmt33,fmt33,fmt33,fmt33,fmt44,fmt44,fmt44};

synoN=[2,2,2,2,2,2,2,2,2,3,4,4,4,4,4,6,6,6];

text={'GluE','HisH','GlnQ','PheF','TyrY','CysC','AsnN','LysK','AspD','IleI','ProP','ThrT','AlaA','ValV','GlyG','LeuL','SerS','ArgR'};

for CTcount=1:18
    
    fileName2=[speciesName,text{CTcount},'ratio.txt'];
    fileID2 = fopen(fileName2);
    tt = textscan(fileID2,fmList2{CTcount},1,'delimiter',',');
    xx = textscan(fileID2,fmList{CTcount},'delimiter',',');
    fclose(fileID2);
    
    for xxcount=1:synoN(CTcount)
        XX(:,xxcount)=xx{1,xxcount+1};
    end
    
    
    RefCAIp{CTcount}=mean(XX)./max(mean(XX));
    ObsFrqp{CTcount}=XX;
    clear XX;
    
end
RefCAI=cell2mat(RefCAIp);
ObsFrq=cell2mat(ObsFrqp);


gSize=length(ObsFrq(:,1));

for rowi=1:gSize
    temp=ObsFrq(rowi,:);
    ObsCAI(rowi,:)=temp./sum(temp);
end

ObsCAI(isnan(ObsCAI)) = 0;

end