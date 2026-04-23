% GPR_RTM_Bscan_Demo.m
% 探地雷达逆时偏移(RTM) - 多道扫描与干涉叠加演示
clear; clc; close all;

%% 1. 参数设置与模型构建
nx = 100; nz = 100;    
dx = 0.01; dz = 0.01;  
x = (0:nx-1)*dx; z = (0:nz-1)*dz;

ep = 4 * ones(nz, nx); 
ep(40:45, 45:55) = 16; % 地下目标

c = 0.3;               
v = c ./ sqrt(ep);     

dt = 0.02;             
nt = 450;              
t = (0:nt-1)*dt;

fc = 0.8;              
t0 = 1.5 / fc;
src = (1 - 2*(pi*fc*(t-t0)).^2) .* exp(-(pi*fc*(t-t0)).^2); 

% 【核心修改】：设置一条测线，天线从 x=20 移动到 x=80，共扫描 13 个点
sx_array = 20:5:80; 
num_traces = length(sx_array);
sz = 2; 

Image_Total = zeros(nz, nx); % 用于叠加所有单道RTM的最终图像
v_bg = c ./ sqrt(4 * ones(nz, nx)); % 背景速度
coeff_bg = (v_bg * dt / dx).^2;
coeff = (v * dt / dx).^2;

disp(['开始 B-scan 模拟，共需计算 ', num2str(num_traces), ' 道数据...']);

%% 开始多道循环扫描
for tr = 1:num_traces
    sx = sx_array(tr);
    fprintf('正在处理第 %d/%d 道 (位置: %.2fm)...\n', tr, num_traces, x(sx));
    
    % --- 正向传播 ---
    P_past = zeros(nz, nx); P_now = zeros(nz, nx); P_next = zeros(nz, nx);
    S_field = zeros(nz, nx, nt); 
    B_scan  = zeros(1, nt);      
    
    for i = 1:nt
        P_now(sz, sx) = P_now(sz, sx) + src(i);
        P_next(2:nz-1, 2:nx-1) = 2*P_now(2:nz-1, 2:nx-1) - P_past(2:nz-1, 2:nx-1) + ...
            coeff(2:nz-1, 2:nx-1) .* (P_now(3:nz, 2:nx-1) + P_now(1:nz-2, 2:nx-1) + ...
                                      P_now(2:nz-1, 3:nx) + P_now(2:nz-1, 1:nx-2) - ...
                                      4*P_now(2:nz-1, 2:nx-1));
        B_scan(i) = P_next(sz, sx);
        S_field(:,:,i) = P_now;
        P_past = P_now; P_now = P_next;
    end
    
    % --- 切除直达波 ---
    mute_steps = round(4.0 / dt); 
    B_scan_muted = B_scan;
    B_scan_muted(1:mute_steps) = 0; 
    
    % --- 逆时反传与成像 ---
    R_past = zeros(nz, nx); R_now = zeros(nz, nx); R_next = zeros(nz, nx);
    Image_Single = zeros(nz, nx);       
    
    for i = nt:-1:1
        R_now(sz, sx) = R_now(sz, sx) + B_scan_muted(i);
        R_next(2:nz-1, 2:nx-1) = 2*R_now(2:nz-1, 2:nx-1) - R_past(2:nz-1, 2:nx-1) + ...
            coeff_bg(2:nz-1, 2:nx-1) .* (R_now(3:nz, 2:nx-1) + R_now(1:nz-2, 2:nx-1) + ...
                                         R_now(2:nz-1, 3:nx) + R_now(2:nz-1, 1:nx-2) - ...
                                         4*R_now(2:nz-1, 2:nx-1));
        Image_Single = Image_Single + S_field(:,:,i) .* R_now;
        R_past = R_now; R_now = R_next;
    end
    
    % 将当前道的 RTM 结果叠加到总图像中
    Image_Total = Image_Total + Image_Single;
end

%% 压制低频伪影 (Laplacian 或 垂直差分)
% RTM 互相关通常伴随强烈的低频背景，这里取 z 方向的二阶导数进行图像锐化
Image_Final = del2(Image_Total);

%% 可视化最终成果
disp('计算完成，绘制最终高精度叠加图像！');
figure('Position', [150, 150, 800, 400], 'Name', '多道叠加 RTM 成像');

subplot(1,2,1);
imagesc(x, z, ep); colorbar;
title('真实相对介电常数模型');
xlabel('距离 (m)'); ylabel('深度 (m)'); 
axis image; set(gca,'FontSize',12);

subplot(1,2,2);
cmax = max(abs(Image_Final(:))) * 0.5; 
imagesc(x, z, Image_Final, [-cmax cmax]); 
colormap('gray'); colorbar;
title('多道叠加 RTM 成像结果');
xlabel('距离 (m)'); ylabel('深度 (m)'); 
axis image; set(gca,'FontSize',12);

%% 压制低频伪影 (Laplacian 锐化)
Image_Final = del2(Image_Total);

% ================= 新增：高端后处理管线 =================

% 1. 浅层强伪影切除 (Mute source footprint)
% 地表附近的强奇异性通常没有有效信息，直接将地表以下 8cm 置零
mute_depth_idx = round(0.08 / dz); 
Image_Final(1:mute_depth_idx, :) = 0;

% 2. 提取瞬时振幅包络 (Hilbert Transform)
% 消除子波的黑白相间震荡，提取纯粹的能量聚焦团
Image_Env = zeros(nz, nx);
for j = 1:nx
    % 对每一列（深度方向）做 Hilbert 变换求包络
    Image_Env(:, j) = abs(hilbert(Image_Final(:, j))); 
end

%% 可视化最终成果
disp('计算完成，绘制最终包络图像！');
figure('Position', [150, 150, 800, 400], 'Name', '多道叠加 RTM 包络成像');

subplot(1,2,1);
imagesc(x, z, ep); colorbar;
title('真实相对介电常数模型');
xlabel('距离 (m)'); ylabel('深度 (m)'); 
axis image; set(gca,'FontSize',12);

subplot(1,2,2);
% 包络是没有负数的，所以色标从 0 开始
cmax = max(Image_Env(:)) * 0.6; % 调节 0.6 这个系数可以控制图像对比度
imagesc(x, z, Image_Env, [0 cmax]); 
colormap('gray'); colorbar;
title('RTM 成像结果 (去浅层 + 能量包络)');
xlabel('距离 (m)'); ylabel('深度 (m)'); 
axis image; set(gca,'FontSize',12);