clear all;clc;close all;
T =2;%ɨ������
M=5;%���ؿ���ģ����������ڲ��ԣ������˲�����ֵ������׼��
N = 900/T;%�ܵò��õ���
N1 = 400/T;%��һ��ת�䴦�������
N2 = 600/T;%��һ�����ٴ��������
N3 = 610/T;%�ڶ�ת�䴦�������
N4 = 660/T;%�ڶ����ٴ��������
Delta = 100;%����������׼��
Rx = zeros(N,1);
Ry = zeros(N,1);
Zx = zeros(N,M);
Zy = zeros(N,M);
%��y������ֱ���˶�
t = 2:T:400;
x0 = 2000 + 0*t';
y0 = 10000 - 15*t';
%��ת��
t= 400+T:T:600;
x1 = x0(N1) + 0.075*((t'-400).^2)/2;
y1 = y0(N1) - 15*(t'-400)+0.075*((t'-400).^2)/2;
%����
t = 600+T:T:610;
vx = 0.075*(600-400);
x2 = x1(N2-N1)+vx*(t'-600);
y2 = y1(N2-N1)+0*t';
%��ת��
t = 610+T : T:660;
x3=x2(N3-N2)+(vx*(t'-610)-0.3*((t'-610).^2)/2);
y3=y2(N3-N2)-0.3*((t'-610).^2)/2;
%����
t=660+T:T:900;
vy=-0.3*(660-610);
x4=x3(N4-N3)+0*t';
y4=y3(N4-N3)+vy*(t'-660);
%�켣�ϳ�
Rx = [x0;x1;x2;x3;x4];
Ry = [y0;y1;y2;y3;y4];
%ÿ���������������˲�����λ�õĳ�ʼ��
mtEstPx=zeros(M,N);
mtEstPy=zeros(M,N);
%�����۲����ݣ�Ҫ����M�Σ�������M�εĹ۲�����
nx = randn(N,M)*Delta;%������������
ny = randn(N,M)*Delta;
% Zx = Rx*ones(1,M)+nx;%��ʵ�Ĺ켣�ϵ�������
% Zy = Ry*ones(1,M)+ny;

%���ļ��ж�������
READ_WRITE_FLAG = 0;
if READ_WRITE_FLAG
    fid = fopen('immTestData.txt','w');
    zxr = reshape(Zx,1,M*N);
    zyr = reshape(Zy,1,M*N);
    for i = 1 : length(zxr)
        fprintf(fid,'%f\n ',zxr(i));
    end
    for i = 1 : length(zyr)
        fprintf(fid,'%f\n ',zyr(i));
    end
    fclose(fid);
else
    zxy = textread('immTestData.txt','%f');
    zxy = reshape(zxy,N,[]);
    Zx = zxy(:,1:5);
    Zy = zxy(:,6:end);
end

for m = 1:M
    %�˲���ʼ��
    mtEstPx(m,1)=Zx(1,m);%��ʼ����
    mtEstPy(m,1)=Zx(2,m);
    xn(1)=Zx(1,m);%�˲���ֵ
    xn(2)=Zx(2,m);
    yn(1)=Zy(1,m);
    yn(2)=Zy(2,m);
    %�ǻ���ģ�Ͳ���
    phi = [1 T 0 0 ;%�˶�ģ��
                0 1 0 0;
                0 0 1 T;
                0 0 0 1];
   h= [1 0 0 0;%��������
            0 0 1 0];
   g=[T/2 0;
            1 0 ;
            0 T;
            0 1];
%    q=0.01*[];
   R=[Delta.^2 0 ;%
            0 Delta.^2];
    vx=(Zx(2,m)-Zx(1,m))/2;
    vy=(Zy(2,m)-Zy(1,m))/2;
    %��ʼ״̬����
    x_est=[Zx(2,m);vx;Zy(2,m);vy];
    p_est=[Delta^2,Delta^2/T,0,0;
                    Delta^2/T,2*Delta^2/(T^2),0,0;
                    0,0,Delta^2,Delta^2/T;
                    0,0,Delta^2/T,2*Delta^2/(T^2)];
     mtEstPx(m,2)=x_est(1);
     mtEstPy(m,2)=x_est(3);
     %�˲���ʼ
     for r=3:N
        z=[Zx(r,m);Zy(r,m)];
        if r == 3%ǰ20��ֻ���ǻ���ģ���˲�
            [x_est,p_est] = kf(x_est,p_est,z,phi,h,R);
            xn(r)=x_est(1);%��¼�������˲�����
            yn(r)=x_est(3);
            mtEstPx(m,r)=x_est(1);%��¼��m�η����˲���������
            mtEstPy(m,r)=x_est(3);
        else
            if r==4
                x_est=[x_est;0;0];%��ά
                p_est=p_est;
                p_est(6,6)=0;%��ά
                for i =1:3
                    xn_est{i,1}=x_est;
                    pn_est{i,1}=p_est;
                end
                u=[0.8,0.1,0.1];%ģ�͸��ʳ�ʼ��
            end
            [x_est,p_est,xn_est,pn_est,u]=imm(xn_est,pn_est,T,z,Delta,u);
            xn(r)=x_est(1);
            yn(r)=x_est(3);
            mtEstPx(m,r)=x_est(1);
            mtEstPy(m,r)=x_est(3);
        end
     end%����һ���˲�
% figure(m);
% plot(Rx,Ry,'r',Zx,Zy,'g*',xn,yn,'b*');
% legend('��ʵ�켣','�۲�����','���ƹ켣');
end
err_x=zeros(N,1);
err_y=zeros(N,1);
delta_x = zeros(N,1);
delta_y=zeros(N,1);
for r = 1:N
    ex=sum(Rx(r)-mtEstPx(:,r));
    ey=sum(Ry(r)-mtEstPy(:,r));
    err_x(r)=ex/M;
    err_y(r)=ey/M;
    eqx=sum((Rx(r)-mtEstPx(:,r)).^2);
    eqy=sum((Ry(r)-mtEstPy(:,r)).^2);
    
    delta_x(r)=sqrt(abs(eqx/M-(err_x(r)^2)));
    delta_y(r)=sqrt(abs(eqy/M-(err_y(r)^2)));
end
figure(1);
plot(Rx,Ry,'r',Zx,Zy,'g*',xn,yn,'b+');
legend('��ʵ�켣','�۲�����','���ƹ켣');
figure(2);
subplot(211);
plot(err_x);
title('average of Error of x dirction ');
subplot(212);
plot(err_y);
title('average of Error of y dirction ');
figure(3);
subplot(211);
plot(delta_x);
title('standard deviation  of Error of x dirction ');
subplot(212);
plot(delta_y);
title('standard deviation  of Error of y dirction ');
