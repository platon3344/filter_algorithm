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
%%%%%  ��ͼ  %%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%  1.�˲����  %%%
%%%%%%%%%%%%%%%%%%%%
figure;
%Ŀ��a,b�Ĺ۲�λ��
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
%Ŀ��a,b���Ӳ�λ��
for i=1:c
    if i==1
        plot(NOISE_sum_a(1,:),NOISE_sum_a(2,:),'b.');
    else
        plot(NOISE_sum_b(1,:),NOISE_sum_b(2,:),'r.');
    end
end
hold on;
%Ŀ��a,b�Ĺ���λ��
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
legend('Ŀ��a�Ĺ۲�λ��','Ŀ��b�Ĺ۲�λ��','Ŀ��a���Ӳ�','Ŀ��b���Ӳ�','Ŀ��a�Ĺ���λ��','Ŀ��b�Ĺ���λ��');grid;

