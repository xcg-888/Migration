function [Relative_permittivity1,Relative_permittivity2] = angle_hough_transform(B_scan_image_edge_dilate,B_scan_image_edge,Downsample_N,dt,dx)
beta = 1.3753;
k = 2; % 找出霍夫空间中投票最多的前k个角度
theta_range = -70:1:70;  % 角度范围缩小，避开横向的介质层杂波
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算霍夫变换，限制直线法向角度范围为theta_range
[H1, theta1, ~] = hough(B_scan_image_edge_dilate, 'Theta', theta_range);  
angle_H1 = var(H1,0,1);
[~,max_theta1] = mink(angle_H1,k);
k_line1 = -cot(theta1(max_theta1)/180*pi);
disp(['角度：' num2str(theta1(max_theta1))]);
disp(['斜率：' num2str(k_line1)]);
c = 3e8;
Relative_permittivity1 = mean((beta*abs(k_line1)*dt*Downsample_N*c/2/dx).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[H2, theta2, ~] = hough(B_scan_image_edge, 'Theta', theta_range);  
angle_H2 = var(H2,0,1);
[~,max_theta2] = mink(angle_H2,k);
k_line2 = -cot(theta2(max_theta2)/180*pi);
disp(['角度：' num2str(theta2(max_theta2))]);
disp(['斜率：' num2str(k_line2)]);
c = 3e8;
Relative_permittivity2 = mean((beta*abs(k_line2)*dt*Downsample_N*c/2/dx).^2);

end