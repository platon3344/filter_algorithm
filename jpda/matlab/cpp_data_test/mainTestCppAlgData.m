clc;clear;close all;
y1 = readFilterInputData("jpdaTestData.txt");
NOISE_sum_a = [];
NOISE_sum_b = [];
for i = 1 : length(y1)
    NOISE_sum_a = [NOISE_sum_a y1{i}(:,1)];
    NOISE_sum_b = [NOISE_sum_b y1{i}(:,2)];
    data_measurement(:,:,i) = y1{i}(:,end-1:end);
end
estData = readFilterInputData("jpdaEstData.txt");
for i = 1 : length(estData)
    x_Estimate(:,:,i) = estData{i};
end

%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%  画图  %%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%  1.滤波结果  %%%
%%%%%%%%%%%%%%%%%%%%
figure;
%目标a,b的观测位置
c = 2;
n = length(estData);
for i=1:c
    a=zeros(1,n);
    b=zeros(1,n);
    for j=1:n
        a(j)=data_measurement(i,1,j);
        b(j)=data_measurement(i,2,j);
    end
    if i==1
        plot(a(:),b(:),'bo')
    else
        plot(a(:),b(:),'ro')
    end
    hold on;
end
%目标a,b的杂波位置
for i=1:c
    if i==1
        plot(NOISE_sum_a(1,:),NOISE_sum_a(2,:),'b.');
    else
        plot(NOISE_sum_b(1,:),NOISE_sum_b(2,:),'r.');
    end
end
hold on;
%目标a,b的估计位置
for i=1:c
    a=zeros(1,n);
    b=zeros(1,n);
    for j=1:n
        a(j)=x_Estimate(1,i,j);
        b(j)=x_Estimate(2,i,j);
    end
    if i==1
        plot(a(:),b(:),'b*');
    else
        plot(a(:),b(:),'r*');
    end
    hold on;
end
xlabel('x/m'),ylabel('y/m');
legend('目标a的观测位置','目标b的观测位置','目标a的杂波','目标b的杂波','目标a的估计位置','目标b的估计位置');grid;

