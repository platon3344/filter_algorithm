%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������  :����JPDA���ݹ����㷨ʵ�������������˶�Ŀ��ĵ㼣�뺽���Ĺ���
% �������  :
%           -target_position:  Ŀ��ĳ�ʼλ��
%           - n:  ��������
%           - T:  �������
%           -MC_number:�������
%           - c:  Ŀ�����
% �������  :
%           ��
% �ο�����  :
%           ����,�����ھ��ںϼ����о���Ӧ��,������ҵ��ѧ˶ʿѧλ����,2004��
% ����      ��
%           �ô���Ϊ���߱�ҵ������ݣ�����ѧ�������ĽǶȣ����ڹ��������ô���
%           �ô���Ǳ���ԭ�����޸���������һλ���ߵ�JPDA����
%           �ô��������ѧ�����������������κ�������ҵ��;�������Ծ�����
%           ��������øô�����в����ʵ���;���ô������߲��е��κ�����
%           ���������ߵ��Ͷ��ɹ���ת�������
% ��������  ��
%           wangzexun@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%  ��������  %%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
clc;close all ;clear;
n = 20;%�������ݸ���
T = 1;%����Ϊ1s
MC_number = 1;%montecarlo����
c = 2;%Ŀ�����

Pd=1;                                                                              %������
Pg=0.99;                                                                           %��ȷ��������������ڵø���
g_sigma=9.21;                                                                      %����
lambda=2;
gamma=lambda*10^(-6);                                                              %ÿһ����λ���(km^2)�ڲ���lambda���Ӳ�
Target_measurement=zeros(c,2,n);                                                   %Ŀ��۲⻥���洢����
target_delta=[100 100];                                                            %Ŀ���Ӧ�Ĺ۲��׼��
P=zeros(4,4,c);                                                                    %Э�������
P1=[target_delta(1)^2 0 0 0;0 0.01 0 0;0 0 target_delta(1)^2 0;0 0 0 0.01];        %��ʼЭ�������
P(:,:,1)=P1;
P(:,:,2)=P1;
phi = [1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];                                             %״̬ת�ƾ���
H = [1 0 0 0;0 0 1 0];                                                             %�۲����
R=[target_delta(1)^2 0;0 target_delta(1)^2];                                       %�۲�Э�������
Q=[4 0;0 4];                                                                       %ϵͳ��������Э����
G=[T^2/2 0;T 0;0 T^2/2;0 T];                                                       %������������
x_Estimate=zeros(4,c,n);                                                             %�洢Ŀ��ĸ�ʱ�̵��˲�ֵ
x_filter1=zeros(4,c,n,MC_number);                                                  %MC_number��Montle Carlo��������ȫ������洢����
data_measurement=zeros(c,2,n);                                                     %�۲�洢����
data_measurement1=zeros(c,4,n);                                                    %ʵ��λ������x,y����
target_position = [1500 300 500 400;
                            500 400 1500 300];                              %Ŀ�����ʼλ�ú��ٶ�(m,m/s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  ����Ŀ���ʵ��λ��  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_measurement1(:,:,1)=target_position;                                          %ʵ��λ�þ����ʼ��
for i=1:c
    for ii=2:n                                                                     %ʵ��λ��
        data_measurement1(i,:,ii)=(phi*data_measurement1(i,:,ii-1)')'+(G*sqrt(Q)*(randn(2,1)))';
    end
end

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%  ��������  %%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

for M=1:MC_number
    %%
    y1 = cell(1,n); %��ʵ����������
    NOISE_sum_a = [];%�Ӳ�����
    NOISE_sum_b = [];
    READ_WRITE_FLAG = 1;
    if READ_WRITE_FLAG == 0
        %%%%%%%%%%%%%%%%%%%%
        %%%  1.������������  %%%
        %%%%%%%%%%%%%%%%%%%%
        Noise=[];
        for i=1:n
            for j=1:c                                                                      %���������۲��λ��
                data_measurement(j,1,i)=data_measurement1(j,1,i)+rand(1)*target_delta(j);%ȡ����x
                data_measurement(j,2,i)=data_measurement1(j,3,i)+rand(1)*target_delta(j); %ȡ����y
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  2.�����Ӳ�,��ȷ����Ч�۲�  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        S=zeros(2,2,c);
        Z_predic=zeros(2,2);                                                               %�洢����Ŀ��Ĺ۲�Ԥ��ֵ,��ֻ����x,y����
        x_predic=zeros(4,2);                                                               %�洢����Ŀ���״̬Ԥ��ֵ,������x,y�����x,y�����ٶ�
        ellipse_Volume=zeros(1,2);
        fid = fopen("jpdaTestData.txt","w");
        %�����ǲ������ݵĲ��֣�
        for t = 1: n 
            Noise=[];
            for i=1:c %��ÿ��Ŀ�����������Ŀ��
                if t~=1
                    x_predic(:,i) = phi*data_measurement1(i,:,t-1)';                                       %��ǰһʱ�̵��˲�ֵ��Ԥ�⵱ǰ��ֵ(kalman�˲��ĵ�һ�����ʽ)
                else
                    x_predic(:,i)=target_position(i,:)';                                       %��һ�β�����������ʵλ�õ�Ԥ��ֵ
                end
                P_predic(:,:,i)=phi*P(:,:,i)*phi'+G*Q*G';                                                 %����x_predicЭ�������(kalman�˲��ĵڶ������ʽ)
                Z_predic(:,i)=H*x_predic(:,i);                                                 %��ȡԤ��ֵ��x,y���꣬����x,y�ٶ�
                R=[target_delta(i)^2 0; 0 target_delta(i)^2];
                S(:,:,i)=H*P_predic(:,:,i)*H'+R;                                                      %�����м����S
                ellipse_Volume(i)=pi*g_sigma*sqrt(det(S(:,:,i)));                              %������Բ�����ŵ����
                number_returns=floor(ellipse_Volume(i)*gamma+1);                               %��Բ�������ڵĴ���ز���
                side=sqrt((ellipse_Volume(i)*gamma+1)/gamma)/2;                                %����Բ�����ŵ�ЧΪ�����Σ�����������α߳��Ķ���֮һ
                Noise_x=x_predic(1,i)+side-2*rand(1,number_returns)*side;                      %��Ԥ��ֵ��Χ��������ز���ע�⣺��ĳһ��number_returnsС�ڵ���0ʱ�����������һ�μ��ɡ�
                Noise_y=x_predic(3,i)+side-2*rand(1,number_returns)*side;
                Noise=[Noise_x;Noise_y];
    %             NOISE=[NOISE Noise];
                if i==1
                    NOISE_sum_a=[NOISE_sum_a Noise];%ÿ������Ŀ����Ӳ�
                else
                    NOISE_sum_b=[NOISE_sum_b Noise];%ÿ������Ŀ����Ӳ�
                end
            end
            %��Ŀ��ͼ�Ŀ��ϲ�
            b = data_measurement(:,:,t);
            y1{t}=[Noise b'];                                                                    %����һ���Ӳ��ز�ʱ
            writeData = [Noise b'];
            fprintf(fid,"%d , " ,length(writeData));
            for jj = 1 : size(writeData,2)
                for kk = 1 : size(writeData,1)
                    fprintf(fid,"%f ,", writeData(kk,jj));
                end
            end
            fprintf(fid,"\n");        
        end
        fclose(fid);
    else
        data = textread("jpdaTestData.txt","%s","delimiter","\n");
        for i = 1 : length(data)
            lineData = split(data{i},',');
            measureData = [];
            for j = 1: length(lineData)
                measureData = [measureData str2num(lineData{j})];
            end
            targetNum = measureData(1);
            useData = measureData(2:end);
            y1{i} = reshape(useData,2,[]);
            NOISE_sum_a = [NOISE_sum_a y1{i}(:,1)];
            NOISE_sum_b = [NOISE_sum_b y1{i}(:,2)];
            data_measurement(:,:,i) = y1{i}(:,end-1:end);
        end
    end
    
    %% ����Ϊjpda����ִ�в��֡�
    for t=1:n
        for i=1:c %��ÿ��Ŀ�����������Ŀ��
            if t~=1
                x_predic(:,i) = phi*x_Estimate(:,i,t-1);                                       %��ǰһʱ�̵��˲�ֵ��Ԥ�⵱ǰ��ֵ(kalman�˲��ĵ�һ�����ʽ)
            else
                x_predic(:,i)=target_position(i,:)';                                       %��һ�β�����������ʵλ�õ�Ԥ��ֵ
            end
            P_predic(:,:,i)=phi*P(:,:,i)*phi'+G*Q*G';                                                 %����x_predicЭ�������(kalman�˲��ĵڶ������ʽ)
            Z_predic(:,i)=H*x_predic(:,i);                                                 %��ȡԤ��ֵ��x,y���꣬����x,y�ٶ�
            R=[target_delta(i)^2 0; 0 target_delta(i)^2];
            S(:,:,i)=H*P_predic(:,:,i)*H'+R;                                                      %�����м����S
            ellipse_Volume(i)=pi*g_sigma*sqrt(det(S(:,:,i)));                              %������Բ�����ŵ����
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  3.�����۲�ȷ�Ͼ���Q2  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m1=0;                                                                          %��¼��Ч�۲����
        y=[];
        inputData = y1{t};
        [n1,n2]=size(inputData);
        Q1=zeros(100,3);
        for j=1:n2 %n2Ϊ��ǰʱ��۲⵽Ŀ�����
            flag=0;
            for i=1:c
                d=inputData(:,j)-Z_predic(:,i);
                D=d'*inv(S(:,:,i))*d;
                if D<=g_sigma
                    flag=1;
                    Q1(m1+1,1)=1;
                    Q1(m1+1,i+1)=1;
                end
            end
            if flag==1
                y=[y inputData(:,j)];                                                      %������������е����лز�����y��
                m1=m1+1;                                                            %��¼��Ч�۲����
            end
        end
        Q2=Q1(1:m1,1:3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  4.������������A_matrix,����num��ʾ���������¼�����  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_matrix=zeros(m1,3,10000);
        A_matrix(:,1,1:10000)=1;%���������һ������Ϊ1
        if m1 > 0                                                                       %m1=0��ʾ����Ŀ�궼û�й۲�
            num=1;
            for i=1:m1
                if Q2(i,2)==1 %���Ƕ�1�ź���Ŀ�꣬������������
                    A_matrix(i,2,num)=1;
                    A_matrix(i,1,num)=0;
                    num=num+1;
                    for j=1:m1
                        if (i~=j)&(Q2(j,3)==1)
%                             A_matrix(i,2,num)=1;
%                             A_matrix(i,1,num)=0;
                            A_matrix(j,3,num)=1;
                            A_matrix(j,1,num)=0;
                            num=num+1;
                        end
                    end
                end
            end
            
            for i=1:m1
                if Q2(i,3)==1 %���Ƕ�3�ź���Ŀ�꣬������������
                    A_matrix(i,3,num)=1;
                    A_matrix(i,1,num)=0;
                    num=num+1;
                end
            end
            
%         else
%             flag=1;
        end
        A_matrix=A_matrix(:,:,1:num);                                                  %��ٷ���ֵĽ������A_matrix��
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  5.����������Pr,����False_num��ʾ������,mea_indicator��ʾ�۲�ָʾ��,target_indicator��ʾĿ��ָʾ��  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Pr=zeros(1,num);
        nn = [];
        for i=1:num
            targetNum = size(A_matrix(:,:,i),2);
            False_num=m1;
            N=1;
            for j=1:m1
                mea_indicator=sum(A_matrix(j,2:3,i));                                      %�ο�������ʽ4-48
                if mea_indicator==1
                    False_num=False_num-1;
                    for kk = 2 : targetNum
                        if A_matrix(j,kk,i)==1                                                  %����۲���Ŀ��1����
                            b=(y(:,j)-Z_predic(:,kk-1))'*inv(S(:,:,kk-1))*(y(:,j)-Z_predic(:,kk-1));
                            tmp1 = sqrt(det(2*pi*S(:,:,kk-1)));
                            tmp2 = exp(-1/2*b);
                            N=N/ tmp1* tmp2 ;                          %������̬�ֲ�����
                        end
                    end
                end
                nn = [nn N];
            end
            if Pd==1
                a=1;
            else
                a=1;
                for j=1:c
                    target_indicator=sum(A_matrix(:,j+1,i));                               %�ο�������ʽ4-49
                    a=a*Pd^target_indicator*(1-Pd)^(1-target_indicator);                   %���������
                end
            end
            V=ellipse_Volume(1)+ellipse_Volume(2);                                         %��ʾ������������
            
            a1=1;
            for j=1:False_num
                a1=a1*j;
            end
            Pr(i)=N*a*a1/(V^False_num);
        end
        Pr=Pr/sum(Pr);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  6.�����������U  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        U=zeros(m1+1,c);
        for i=1:c
            for j=1:m1
                for k=1:num
                    U(j,i)=U(j,i)+Pr(k)*A_matrix(j,i+1,k);
                end
            end
        end
        U(m1+1,:)=1-sum(U(1:m1,:));                                                      %��������Ŀ��T�����Ĺ������ʴ���U��m1+1,:),��һ��
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  7.Kalman�˲���ʼ  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:c                                                                          %����Э�������
            P_predic(:,:,i) = phi*P(:,:,i)*phi'+G*Q*G';
            K (:,:,i)= P_predic(:,:,i)*H'*inv(S(:,:,i));
            P(:,:,i)= P_predic(:,:,i)-(1-U(m1+1,i))*K(:,:,i)*S(:,:,i)*K(:,:,i)';
        end
        
        for i=1:c
            %7.1���¹���ֵ
            x_predic_tmp=0;                                                                   %������õ��м����
            for j=1:m1
                x_predic_tmp=x_predic_tmp+U(j,i)*(x_predic(:,i)+ K (:,:,i)*(y(:,j)- Z_predic(:,i)));
            end
            x_predic_tmp=x_predic_tmp + U(m1+1,i)*x_predic(:,i);
            x_Estimate(:,i,t)=x_predic_tmp;
            %7.2����Э�������
            tmpB = 0;
            for j=1:m1+1
                if j==m1+1
                    tmp=x_predic(:,i);
                else
                    tmp=x_predic(:,i)+ K (:,:,i)*(y(:,j)- Z_predic(:,i));
                end
                tmpB=tmpB+U(j,i)*(tmp*tmp'-x_predic_tmp*x_predic_tmp');
            end
            P(:,:,i)=P(:,:,i)+tmpB;
            x_filter1(:,i,t,M)=x_Estimate(:,i,t);
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%  ��ͼ  %%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
x_Estimate=sum(x_filter1,4)/MC_number;                                               %�˲�ֵ��ƽ��
%%%%%%%%%%%%%%%%%%%%
%%%  1.�˲����  %%%
%%%%%%%%%%%%%%%%%%%%
figure;
%Ŀ��a,b�Ĺ۲�λ��
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
        b(j)=x_Estimate(3,i,j);
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
%%
%%%%%%%%%%%%%%%%%%%%
%%%  2.�ٶ����  %%%
%%%%%%%%%%%%%%%%%%%%
figure;
a=0;
c1=zeros(c,n);
for j=1:n
    for i=1:MC_number                                                              %��С�������
        a=(x_filter1(1,1,j,i)-data_measurement1(1,1,j))^2+(x_filter1(3,1,j,i)-data_measurement1(1,3,j))^2;
        c1(1,j)=c1(1,j)+a;
    end
    c1(1,j)=sqrt(c1(1,j)/MC_number);
end
temp=c1(1,:);
a_extra=zeros(2,n);
b_extra=zeros(1,n);
c_extra=zeros(1,n);
a_extra(1,:)=temp;
a_extra(2,:)=1:1:n;
b_extra=a_extra(1,:);
[c_extra,pos]=sort(b_extra);                                                       %posΪ�������±�,cΪ��һ�е�������;
a_extra(2,:)=a_extra(2,pos);                                                       %�ڶ��а��յ�һ��������±��Ӧ
a_extra(1,:)=c_extra;                                                              %��һ�н�����¸���a �ĵ�һ��;
str1=num2str(a_extra(2,n));
str2=num2str(a_extra(1,n));
str=strcat('\itN=',str1,'\itError=',str2,'(m)');
text(a_extra(2,n),0.8*a_extra(1,n),str);
hold on;
plot([a_extra(2,n) a_extra(2,n)],[0 a_extra(1,n)],'r');
hold on;
plot(1:n,c1(1,:),'r:')
hold on;
a=0;
for j=1:n
    for i=1:MC_number                                                              %��С�������
        a=(x_filter1(1,2,j,i)-data_measurement1(2,1,j))^2+(x_filter1(3,2,j,i)-data_measurement1(2,3,j))^2;
        c1(2,j)=c1(2,j)+a;
    end
    c1(2,j)=sqrt(c1(2,j)/MC_number);
end
temp=c1(2,:);
a_extra=zeros(2,n);
b_extra=zeros(1,n);
c_extra=zeros(1,n);
a_extra(1,:)=temp;
a_extra(2,:)=1:1:n;
b_extra=a_extra(1,:);
[c_extra,pos]=sort(b_extra);                                                       %posΪ�������±�,cΪ��һ�е�������;
a_extra(2,:)=a_extra(2,pos);                                                       %�ڶ��а��յ�һ��������±��Ӧ
a_extra(1,:)=c_extra;                                                              %��һ�н�����¸���a �ĵ�һ��;
str1=num2str(a_extra(2,n));
str2=num2str(a_extra(1,n));
str=strcat('\itN=',str1,'\itError=',str2,'(m)');
text(a_extra(2,n),0.8*a_extra(1,n),str);
hold on;
plot([a_extra(2,n) a_extra(2,n)],[0 a_extra(1,n)],'b');
hold on;
plot(1:n,c1(2,:),'b:')
xlabel('times'),ylabel('����ֵ�����ֵ������/m');
legend('Ŀ��a��������ֵ','Ŀ��a�����','Ŀ��b��������ֵ','Ŀ��b�����');grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Revised on 26th June 2008 by wangzexun  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
extra11=zeros(1,n);
extra12=zeros(1,n);
extra13=zeros(1,n);
for j=1:n
    extra11(1,j)=sqrt(x_Estimate(1,1,j)-data_measurement1(1,1,j))^2+(x_Estimate(3,1,j)-data_measurement1(1,3,j))^2;
    extra12(1,j)=sqrt((data_measurement(1,1,j)-data_measurement1(1,1,j))^2+(data_measurement(1,2,j)-data_measurement1(1,3,j))^2);
    extra13(1,j)=extra12(1,j)/extra11(1,j);
end
plot(1:n,extra13(1,:),'k:');
xlabel('times'),ylabel('RMSE of a');
grid;

figure;
extra21=zeros(1,n);
extra22=zeros(1,n);
extra23=zeros(1,n);
for j=1:n
    extra21(1,j)=sqrt(x_Estimate(1,2,j)-data_measurement1(2,1,j))^2+(x_Estimate(3,2,j)-data_measurement1(2,3,j))^2;
    extra22(1,j)=sqrt((data_measurement(2,1,j)-data_measurement1(2,1,j))^2+(data_measurement(2,2,j)-data_measurement1(2,3,j))^2);
    extra23(1,j)=extra22(1,j)/extra21(1,j);
end
plot(1:n,extra23(1,:),'k:');
xlabel('times'),ylabel('RMSE of b');
grid;