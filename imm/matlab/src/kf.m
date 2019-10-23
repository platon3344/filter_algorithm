function [x_est,p_est] = kf(x_est,p_est,z,phi,h,R)

x_pre=phi*x_est;%预测,Warning :该地方未添加过程噪声协方差Q
p_pre=phi*p_est*phi';%预测误差协方差,Warning :该地方未添加过程噪声协方差Q
k=p_pre*h'*inv(h*p_pre*h'+R);
x_est=x_pre+k*(z-h*x_pre);%滤波
p_est=(eye(4)-k*h)*p_pre;%滤波协方差
