function [Relative_permittivity,Relative_permittivity2,alpha,alpha2] = lines_hough_transform(data_out,B_scan_image_edge_dilate,B_scan_image_Mean_cancel,Downsample_N,dt,TrackInterval)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_N = 4;% 寻找最大峰值的个数
beta = 1/0.7773; % 1/α系数
eclipse_line = 9; % 假想εr
min_line_length = 50; % 直线最小长度
min_dx = 5; % 最小断点距离
max_theta = 75; % 霍夫变换最大角度

%% θ采样间距
% 1寻找连通域
stats = regionprops(data_out, 'Area', 'BoundingBox');
% 过滤掉面积小于 10 像素的噪点
valid_idx = [stats.Area] > 10; 
stats = stats(valid_idx);
% 计算每个类外接矩形的对角线长度
L_all = zeros(length(stats), 1);
for k = 1:length(stats)
    bb = stats(k).BoundingBox; % [x_min, y_min, width, height]
    L_all(k) = sqrt(bb(3)^2 + bb(4)^2);
end
% 提取最大对角线并计算自适应步长
L_max = max(L_all);
K = 2; % K倍过采样
delta_theta_adaptive = (1/K) * rad2deg(atan(2/L_max)); % 霍夫空间角度步进幅度，越小，检测幅度越精细
disp(['delta_theta_adaptive:',num2str(delta_theta_adaptive)]);
% delta_theta_adaptive = 0.1;
%% 直线霍夫变换——寻找渐近线（data_out）
theta_range = -max_theta:delta_theta_adaptive:max_theta;  % 角度范围缩小，避开横向的介质层杂波
% 计算霍夫变换，限制直线法向角度范围为theta_range
[H, theta, rho] = hough(data_out, 'Theta', theta_range);  

% 方法1：上下累加
% [rows , cols] = size(H); 
% H = [zeros(1,cols);H(1:rows-1,:)]+[H(2:rows,:);zeros(1,cols)]+H;

% 方法2：上下左右累加
% [rows, cols] = size(H);
% up    = [H(2:rows, :); zeros(1, cols)];
% down  = [zeros(1, cols); H(1:rows-1, :)];
% left  = [H(:, 2:cols), zeros(rows, 1)];
% right = [zeros(rows, 1), H(:, 1:cols-1)];
% H = H + up + down + left + right;

% 方法3：8领域累加
% kernel = [1, 1, 1; 
%           1, 1, 1; 
%           1, 1, 1];
% H = conv2(H, kernel, 'same');

% 方法4：上下两行相加，左右一列相加
% [rows , cols] = size(H); 
% H = [H(:, 2:cols), zeros(rows, 1)]+[zeros(rows, 1), H(:, 1:cols-1)]+[zeros(1,cols);H(1:rows-1,:)]+[H(2:rows,:);zeros(1,cols)]+[zeros(2,cols);H(1:rows-2,:)]+[H(3:rows,:);zeros(2,cols)]+H;

% 方法5：上下两行相加
[rows , cols] = size(H); 
H = [zeros(1,cols);H(1:rows-1,:)]+[H(2:rows,:);zeros(1,cols)]+[zeros(2,cols);H(1:rows-2,:)]+[H(3:rows,:);zeros(2,cols)]+H;

% 检测霍夫矩阵峰值
peaks = houghpeaks(H,max_N,'threshold',ceil(0.3*max(H(:))));   

% 展示霍夫空间并筛选的峰值
% figure;
% imshow(H,[],'XData',theta,'YData',rho,'InitialMagnification','fit');
% xlabel('\theta'), ylabel('\rho');
% axis on, axis normal, hold on;
% plot(theta(peaks(:,2)),rho(peaks(:,1)),'s','color','white');
% title('直线霍夫变换-参数空间（非膨胀）');

figure;% 1. [新增] 将整个 Figure 窗口背景设为白色，避免自带的灰色背景
set(gcf, 'Color', 'w'); % 显示图像
imshow(H, [], 'XData', theta, 'YData', rho, 'InitialMagnification', 'fit');
% 2. [关键] 反转颜色映射
% 默认 gray 是 0(黑)->255(白)。flipud 反转后变为 0(白)->255(黑)。
% 这样 H 矩阵中值越高（投票越多），显示越黑。
colormap(flipud(gray));
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
% 3. [修改] 修改标记点颜色
% 背景变白、线条变黑后，原来的白色('white')标记会看不见。
% 建议改成红色('r')以突出显示，或者蓝色('b')。
% 'LineWidth', 1.5 让标记框在论文里看更清楚。
plot(theta(peaks(:,2)), rho(peaks(:,1)), 's', 'color', 'k', 'LineWidth', 1.5);
% title('直线霍夫变换-参数空间（白底黑线）');
% 4. [可选] 优化字体，使其符合论文标准 (Times New Roman)
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');

k_line1 = -cot(theta(peaks(1:max_N,2))/180*pi);
disp('data_out');
disp(['角度：' num2str(theta(peaks(1:max_N,2)))]);
disp(['斜率：' num2str(k_line1)]);
c = 3e8;
% Relative_permittivity = mean((beta*abs(k_line1)*dt*Downsample_N*c/2/TrackInterval).^2);
Relative_permittivity = mean(beta*abs(k_line1)*dt*Downsample_N*c/2/TrackInterval).^2;
alpha = mean(abs(k_line1))*dt*Downsample_N*c/2/TrackInterval/sqrt(eclipse_line);

% 提取直线段（最小长度MinLength，最大间隙FillGap）
lines = houghlines(data_out, theta, rho, peaks,'FillGap',min_dx,'MinLength', min_line_length);
lines_show_y(B_scan_image_Mean_cancel,1-data_out,lines);


%% 直线霍夫变换——寻找渐近线（非膨胀）
theta_range = -max_theta:delta_theta_adaptive:max_theta;  % 角度范围缩小，避开横向的介质层杂波
% 计算霍夫变换，限制直线法向角度范围为theta_range
[H, theta, rho] = hough(B_scan_image_edge_dilate, 'Theta', theta_range);  
% [H, theta, rho] = hough(B_scan_image_edge_close);     % 计算霍夫变换

% 方法1：上下累加
[rows , cols] = size(H); 
H = [zeros(1,cols);H(1:rows-1,:)]+[H(2:rows,:);zeros(1,cols)]+H;

% 方法2：上下左右累加
% [rows, cols] = size(H);
% up    = [H(2:rows, :); zeros(1, cols)];
% down  = [zeros(1, cols); H(1:rows-1, :)];
% left  = [H(:, 2:cols), zeros(rows, 1)];
% right = [zeros(rows, 1), H(:, 1:cols-1)];
% H = H + up + down + left + right;

% 方法3：8领域累加
% kernel = [1, 1, 1; 
%           1, 1, 1; 
%           1, 1, 1];
% H = conv2(H, kernel, 'same');

% 方法4：上下两行相加，左右一列相加
% [rows , cols] = size(H); 
% H = [H(:, 2:cols), zeros(rows, 1)]+[zeros(rows, 1), H(:, 1:cols-1)]+[zeros(1,cols);H(1:rows-1,:)]+[H(2:rows,:);zeros(1,cols)]+[zeros(2,cols);H(1:rows-2,:)]+[H(3:rows,:);zeros(2,cols)]+H;

% 方法5：上下两行相加
% [rows , cols] = size(H); 
% H = [zeros(1,cols);H(1:rows-1,:)]+[H(2:rows,:);zeros(1,cols)]+[zeros(2,cols);H(1:rows-2,:)]+[H(3:rows,:);zeros(2,cols)]+H;

% 检测霍夫矩阵峰值
peaks = houghpeaks(H,max_N,'threshold',ceil(0.5*max(H(:))));   

% 展示霍夫空间并筛选的峰值
% figure;
% imshow(H,[],'XData',theta,'YData',rho,'InitialMagnification','fit');
% xlabel('\theta'), ylabel('\rho');
% axis on, axis normal, hold on;
% plot(theta(peaks(:,2)),rho(peaks(:,1)),'s','color','white');
% title('直线霍夫变换-参数空间（膨胀）');

figure;% 1. [新增] 将整个 Figure 窗口背景设为白色，避免自带的灰色背景
set(gcf, 'Color', 'w'); % 显示图像
imshow(H, [], 'XData', theta, 'YData', rho, 'InitialMagnification', 'fit');
% 2. [关键] 反转颜色映射
% 默认 gray 是 0(黑)->255(白)。flipud 反转后变为 0(白)->255(黑)。
% 这样 H 矩阵中值越高（投票越多），显示越黑。
colormap(flipud(gray));
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
% 3. [修改] 修改标记点颜色
% 背景变白、线条变黑后，原来的白色('white')标记会看不见。
% 建议改成红色('r')以突出显示，或者蓝色('b')。
% 'LineWidth', 1.5 让标记框在论文里看更清楚。
plot(theta(peaks(:,2)), rho(peaks(:,1)), 's', 'color', 'k', 'LineWidth', 1.5);
% title('直线霍夫变换-参数空间（白底黑线）');
% 4. [可选] 优化字体，使其符合论文标准 (Times New Roman)
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');

k_line2 = -cot(theta(peaks(1:max_N,2))/180*pi);
disp('非膨胀');
disp(['角度：' num2str(theta(peaks(1:max_N,2)))]);
disp(['斜率：' num2str(k_line2)]);
c = 3e8;
% Relative_permittivity2 = mean((beta*abs(k_line2)*dt*Downsample_N*c/2/TrackInterval).^2);
Relative_permittivity2 = mean(beta*abs(k_line2)*dt*Downsample_N*c/2/TrackInterval).^2;
alpha2 = mean(abs(k_line2))*dt*Downsample_N*c/2/TrackInterval/sqrt(eclipse_line);

% 提取直线段（最小长度MinLength，最大间隙FillGap）
lines = houghlines(B_scan_image_edge_dilate, theta, rho, peaks,'FillGap',min_dx,'MinLength', min_line_length);
lines_show_y(B_scan_image_Mean_cancel,1-B_scan_image_edge_dilate,lines);

end