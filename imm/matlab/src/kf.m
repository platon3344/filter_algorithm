function [x_est,p_est] = kf(x_est,p_est,z,phi,h,R)

x_pre=phi*x_est;%Ԥ��,Warning :�õط�δ��ӹ�������Э����Q
p_pre=phi*p_est*phi';%Ԥ�����Э����,Warning :�õط�δ��ӹ�������Э����Q
k=p_pre*h'*inv(h*p_pre*h'+R);
x_est=x_pre+k*(z-h*x_pre);%�˲�
p_est=(eye(4)-k*h)*p_pre;%�˲�Э����
