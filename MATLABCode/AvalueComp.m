%%%%%%% demonstrate parameter 'A' value advantagaes regarding 'different
%%%%%%% lengths' & 'different synonymous codons'

%%%0.8189 A1=entropy4Location(4,5,0.8189);
%%%0.8187 A2=entropy4Location(4,30,0.8187);
%%%0.8186 A3=entropy4Location(4,100,0.8186);


[edgesNew1,NHISTfinal1,CumSumNHIST1]=SnAinLen(2,30) ;
[edgesNew2,NHISTfinal2,CumSumNHIST2]=SnAinLen(3,30) ;
[edgesNew3,NHISTfinal3,CumSumNHIST3]=SnAinLen(4,30) ;
[edgesNew4,NHISTfinal4,CumSumNHIST4]=SnAinLen(6,30) ;


id1=find(edgesNew1>=0.2,1);
A1=CumSumNHIST1(id1);
id2=find(edgesNew2>=0.2,1);
A2=CumSumNHIST2(id2);
id3=find(edgesNew3>=0.2,1);
A3=CumSumNHIST2(id3);
id4=find(edgesNew4>=0.2,1);
A4=CumSumNHIST2(id4);
        
figure
subplot(4,1,1)
plot(edgesNew1,CumSumNHIST1,'r.');
hold on
plot(0.2,0,'g*');
ylabel('cummulative probability');
legend('syno=2',['Sn=0.2,','A=',num2str(A1)]);
title('Parameter A value comparison among lengths');
hold off

subplot(4,1,2)
plot(edgesNew2,CumSumNHIST2,'g.');
hold on
plot(0.2,0,'b*');
ylabel('cummulative probability');
legend('syno=3',['Sn=0.2,','A=',num2str(A2)]);
hold off

subplot(4,1,3)
plot(edgesNew3,CumSumNHIST3,'b.');
hold on
plot(0.2,0,'m*');
ylabel('cummulative probability');
legend('syno=4',['Sn=0.2,','A=',num2str(A3)]);

subplot(4,1,4)
plot(edgesNew4,CumSumNHIST4,'m.');
hold on
plot(0.2,0,'r*');
ylabel('cummulative probability');
legend('syno=6',['Sn=0.2,','A=',num2str(A3)]);

xlabel('exhaustive Sn values');
hold off
% 

figure
subplot(4,1,1)
plot(edgesNew1,NHISTfinal1,'r.');
hold on
plot(0.2,0,'g*');
ylabel('probability');
legend('syno=2',['Sn=0.2,','A=',num2str(A1)]);
title('Parameter A value comparison among lengths');
hold off

subplot(4,1,2)
plot(edgesNew2,NHISTfinal2,'g.');
hold on
plot(0.2,0,'b*');
ylabel('probability');
legend('syno=3',['Sn=0.2,','A=',num2str(A2)]);
hold off

subplot(4,1,3)
plot(edgesNew3,NHISTfinal3,'b.');
hold on
plot(0.2,0,'m*');
ylabel('probability');
legend('syno=4',['Sn=0.2,','A=',num2str(A3)]);

subplot(4,1,4)
plot(edgesNew4,NHISTfinal4,'m.');
hold on
plot(0.2,0,'r*');
ylabel('cummulative probability');
legend('syno=6',['Sn=0.2,','A=',num2str(A3)]);

xlabel('exhaustive Sn values');
hold off