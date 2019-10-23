clc;clear;close all;

N = 450;
zxy = textread('immTestData.txt','%f');
zxy = reshape(zxy,N,[]);
Zx = zxy(:,1:5);
Zy = zxy(:,6:end);

estx = textread('immEstDatax.txt','%f');
estx = reshape(estx,N,[]);

esty = textread('immEstDatay.txt','%f');
esty = reshape(esty,N,[]);

figure(1);
plot(Zx(:,5),Zy(:,5),'g*',estx,esty,'bo');
legend('观察样本','估计轨迹');

figure(2);
plot(esty);