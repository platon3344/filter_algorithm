%x_est,p_est,�Ƿ��ص�m�η��棬��r����������˲����
%xn_est,pn_est,��¼ÿ��ģ�Ͷ�Ӧ��m�η��棬��r����������˲����
%T�״�����
%Z����ֵ
%Delta������������
%u��ʼģ��ѡ����ʾ���
function [x_est,p_est,xn_est,pn_est,u]=imm(xn_est,pn_est,T,Z,Delta,u)
p=[0.95 0.025 0.025;%����ģ��ת��������Ʒ�����ת�Ƹ��ʾ���
        0.025 0.95 0.025;
        0.025 0.025 0.95];
%����������ͬ��ģ�Ͳ�����ģ��1���ǻ���ģ�ͣ�ģ��2��3Ϊ����ģ�ͣ�Qֵ��ͬ��    
phi{1,1}=[1 T 0 0 ;
                0 1 0 0;
                0 0 1 T;
                0 0 0 1];
phi{1,1}(6,6)=0;%ģ��1,״̬ת�Ʒ��̡�

phi{2,1}=[1 T 0 0 T^2/2 0;
                0 1 0 0 T 0;
                0 0 1 T 0 T^2/2;
                0 0 0 1 0 T;
                0 0 0 0 1 0 ;
                0 0 0 0 0 1];%ģ��2,״̬ת�Ʒ��̡�
            
phi{3,1}=phi{2,1};%ģ��3,״̬ת�Ʒ��̡�

g{1,1}=[T/2,0;
                1,0;
                0 T/2;
                0 1];
g{1,1}(6,2)=0;%ģ��1

g{2,1}=[T^2/4,0;
                T/2,0;
                0,T^2/4;
                0 T/2;
                1 0 ;
                0 1];%ģ��2
            
g{3,1}=g{2,1};%ģ��3

q{1,1}=zeros(2);
q{2,1}=0.001*eye(2);
q{3,1}=0.014*eye(2);
% Q = g{j,1}*q{j,1}*g{j,1}';
H = [1 0 0 0 0 0;
        0 0 1 0 0 0 ];
R = eye(2)*Delta^2;%��������Э�������
c_mean = zeros(1,3);%��һ������
mu=zeros(3,3);%��ϸ��ʾ���
for i = 1:3
    c_mean=c_mean+p(i,:)*u(i);
end
for i = 1:3
    mu(i,:)=p(i,:)*u(i)./c_mean;
end
%��������
for j = 1:3
    x0{j,1}=zeros(6,1);
    p0{j,1}=zeros(6);
    for i = 1:3
        x0{j,1}=x0{j,1}+xn_est{i,1}*mu(i,j);
    end
    for i =1:3
        p0{j,1}=p0{j,1}+mu(i,j)*(pn_est{i,1}...
                    +(xn_est{i,1}-x0{j,1})*(xn_est{i,1}-x0{j,1})');
    end
end
%ģ�������˲�
a=zeros(1,3);
for j = 1:3
    x_pre{j,1}=phi{j,1}*x0{j,1};
    p_pre{j,1}=phi{j,1}*p0{j,1}*phi{j,1}'+g{j,1}*q{j,1}*g{j,1}';%  Q = g{j,1}*q{j,1}*g{j,1}';
    k{j,1}=p_pre{j,1}*H'*inv(H*p_pre{j,1}*H'+R);
    xn_est{j,1}=x_pre{j,1}+k{j,1}*(Z-H*x_pre{j,1});
    pn_est{j,1}=(eye(6)-k{j,1}*H)*p_pre{j,1};
end
%ģ�͸��ʸ���
for j = 1:3
    v{j,1}=Z-H*x_pre{j,1};%��Ϣ
    s{j,1}=H*p_pre{j,1}*H'+R;%�۲�Э�������
    n=length(s{j,1})/2;
    coeff =1/((2*pi)^n*sqrt(det(s{j,1}))); 
    a(1,j)=coeff *exp(-0.5*v{j,1}'*inv(s{j,1})*v{j,1});%�۲������ģ��j����Ȼ����
end
c=sum(a.*c_mean);%��һ������
u=a.*c_mean./c;%����ģ�͸���

%�������
xn=zeros(6,1);
pn=zeros(6);
for j = 1:3
    xn=xn+xn_est{j,1}.*u(j);
end
for j = 1:3
    pn=pn+u(j).*(pn_est{j,1}+(xn_est{j,1}-xn)*(xn_est{j,1}-xn)');
end
%�����˲����
x_est=xn;
p_est=pn;

