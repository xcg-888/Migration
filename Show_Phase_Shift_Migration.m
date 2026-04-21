% =========================================================================
% 相移偏移 (Phase-Shift Migration) 全流程步骤拆解与动态频域透视
% 特点：保留所有静态中间过程图 + 增加频域相位旋转的动态四面板演示
% =========================================================================
clear; clc; close all;

%% ================= 第0步：生成模拟测试数据 =================
disp('正在生成模拟探地雷达(GPR) B-scan 数据...');
% 1. 设定物理参数 (严格使用国际标准单位：米、秒、赫兹)
c = 3e8;                % 光速 (m/s)
epsilon_r = 9;          % 相对介电常数 (假设为湿土)
v = c / sqrt(epsilon_r);% 介质中的真实电磁波速 (1e8 m/s)
v_esm = v / 2;          % 爆炸反射面等效速度 (双程转单程)

% 2. 设定探测网格
Nx = 300;               % 空间采样点数 (道数)
dx = 0.02;              % 空间采样间隔 (2 cm)
x = (0:Nx-1) * dx;      % 空间坐标轴 (m)
Nt = 512;               % 时间采样点数
dt = 0.1e-9;            % 时间采样间隔 (0.1 ns)
t = (0:Nt-1) * dt;      % 时间坐标轴 (s)

% 3. 在地下放置一个“点状散射体”
target_x = Nx/2 * dx;   % 目标水平位置 (正中间, 3m处)
target_z = 1.5;         % 目标真实深度 (1.5m处)

% 4. 生成双曲线回波数据 (使用 Ricker 子波)
f0 = 800e6;             % 雷达中心频率 (800 MHz)
data_in = zeros(Nt, Nx);
for i = 1:Nx
    % 计算每一道天线到目标的真实几何距离
    dist = sqrt((x(i) - target_x)^2 + target_z^2);
    % 计算双程旅行时
    travel_time = 2 * dist / v; 
    % 叠加雷达子波
    tau = t - travel_time;
    pulse = (1 - 2*(pi*f0*tau).^2) .* exp(-(pi*f0*tau).^2);
    data_in(:, i) = pulse;
end

% [可视化] 原始 B-scan 图像
figure('Name', '步骤0：原始雷达数据 (时空域)', 'Position', [50, 500, 500, 400]);
imagesc(x, t*1e9, data_in);
colormap(gca, 'gray'); title('原始数据 D(x, t) - 呈现典型双曲线');
xlabel('水平位置 x (m)'); ylabel('双程时间 t (ns)');
set(gca, 'YDir', 'reverse'); % 雷达图习惯时间轴朝下

%% ================= 第1步：傅里叶变换到 f-k 域 =================
disp('步骤1：正在将数据转换至频率-波数域 (f-k 域)...');
% 1. 构建严格的离散频率轴 f 和角频率 w
df = 1 / (Nt * dt);
f = (0:Nt-1)' * df;
f(f >= 1/(2*dt)) = f(f >= 1/(2*dt)) - 1/dt; % 负频折叠
w = 2 * pi * f; 

% 2. 构建严格的离散水平波数轴 kx
dk = 1 / (Nx * dx);
k_x = (0:Nx-1) * dk;
k_x(k_x >= 1/(2*dx)) = k_x(k_x >= 1/(2*dx)) - 1/dx;
k_x = 2 * pi * k_x; 

% 3. 执行二维快速傅里叶变换
P_wk = fft2(data_in); 

% [可视化] 频率-波数域振幅谱
figure('Name', '步骤1：频率-波数域 (f-k 域)', 'Position', [550, 500, 500, 400]);
imagesc(fftshift(k_x), fftshift(f)/1e6, abs(fftshift(P_wk)));
colormap(gca, 'jet'); colorbar; title('二维傅里叶变换振幅谱 |P(\omega, k_x)|');
xlabel('水平波数 k_x (rad/m)'); ylabel('频率 f (MHz)');
ylim([-1500, 1500]); % 只看主要频带

%% ================= 第2步：计算波场外推参数与倏逝波处理 =================
disp('步骤2：正在构建深度网格与倏逝波滤波器...');
dz = v_esm * dt;  % 为了保持网格不畸变，深度步长由速度和时间步长决定
Nz = Nt;          % 深度层数与时间点数一致
[Kx_grid, W_grid] = meshgrid(k_x, w);

% 计算波束内部平方根项：(w/v)^2 - kx^2
inner_term = (W_grid ./ v_esm).^2 - Kx_grid.^2;
% 构建倏逝波遮罩 (Evanescent Mask)
evanescent_mask = (inner_term < 0);
% 计算垂直波数 kz (倏逝波区域强制置零)
kz = zeros(size(inner_term));
kz(~evanescent_mask) = sqrt(inner_term(~evanescent_mask));

% [可视化] 倏逝波与传播波的边界
figure('Name', '步骤2：波场状态分区', 'Position', [1050, 500, 500, 400]);
imagesc(fftshift(k_x), fftshift(f)/1e6, fftshift(double(~evanescent_mask)));
colormap(gca, 'bone'); title('波场状态：白色为传播波，黑色为倏逝波(将滤除)');
xlabel('水平波数 k_x (rad/m)'); ylabel('频率 f (MHz)');
ylim([-1500, 1500]);

%% ================= 第3步：逐层向下延拓与 t=0 成像 (静态运算) =================
disp('步骤3：正在执行逐层波场向下延拓与成像 (静态运算)...');
migrated_img_static = zeros(Nz, Nx);
P_current = P_wk; 
K = exp(1j * sign(W_grid) .* kz * dz);
K(evanescent_mask) = 0; % 滤除倏逝波

h = waitbar(0, '正在生成静态偏移结果...');
for iz = 1:Nz
    P_current = P_current .* K;           % 向下延拓一层
    P_kx_t0 = sum(P_current, 1);          % t=0 成像条件 (沿频率轴求和)
    migrated_img_static(iz, :) = real(ifft(P_kx_t0)); % 转换回空间域
    if mod(iz, 50) == 0, waitbar(iz/Nz, h); end
end
close(h);
migrated_img_static = migrated_img_static / max(abs(migrated_img_static(:)));

% [可视化] 最终的偏移聚焦图像
z_axis = (0:Nz-1) * dz; 
figure('Name', '步骤3：最终偏移成像 (空深域)', 'Position', [50, 50, 500, 400]);
imagesc(x, z_axis, migrated_img_static);
colormap(gca, 'gray'); title('静态最终相移偏移结果 M(x, z)');
xlabel('水平位置 x (m)'); ylabel('真实深度 z (m)');
set(gca, 'YDir', 'reverse');
hold on; plot(target_x, target_z, 'ro', 'MarkerSize', 10, 'LineWidth', 2); hold off;
legend('偏移图像', '物理坐标');

%% ================= 第4步：核心动态展示 (四域联动与相位演化) =================
disp('步骤4：准备进入动态重播，请观察弹出的四面板大图...');
fig_anim = figure('Name', '步骤4：频域相移原理全景透视 (动态)', 'Position', [550, 50, 1000, 700]);

migrated_img_dyn = zeros(Nz, Nx);
P_current_dyn = P_wk;  % 重新初始化波场
render_step = 4;       % 渲染步长，保证动画流畅度

% 预先计算用于 f-k 域画图的网格中心化坐标
k_x_shift = fftshift(k_x);
f_ax_shift = fftshift(f)/1e6;

for iz = 1:Nz
    % [核心公式] 频域向下推演
    P_current_dyn = P_current_dyn .* K;
    
    % [成像提取] 提取 t=0 作为当前深度的成像行
    P_kx_t0_dyn = sum(P_current_dyn, 1); 
    migrated_img_dyn(iz, :) = real(ifft(P_kx_t0_dyn)); 
    
    % --- 动态绘图刷新 ---
    if mod(iz, render_step) == 0 || iz == Nz
        
        % 将当前下推后的频域波场，完整变换回时空域用于展示
        current_wavefield_xt = real(ifft2(P_current_dyn));
        
        % [左上] 时空波场 (x-t) - 观察双曲线如何上移并坍缩
        subplot(2, 2, 1);
        imagesc(x, t*1e9, current_wavefield_xt); colormap(gca, 'gray'); caxis([-0.3 0.3]);
        title(sprintf('虚拟时空波场 (天线等效下推深度 z = %.2f m)', z_axis(iz)));
        xlabel('x (m)'); ylabel('t (ns)'); yline(0, 'r', 'LineWidth', 2, 'Label', 't=0 成像线');
        
        % [右上] 偏移成像结果 (x-z) - 观察图像如何逐层生成
        subplot(2, 2, 2);
        imagesc(x, z_axis, migrated_img_dyn); colormap(gca, 'gray'); caxis([-max(migrated_img_dyn(:))*0.3, max(migrated_img_dyn(:))*0.3]);
        title('最终成像 M(x, z) 逐层构建');
        xlabel('x (m)'); ylabel('z (m)');
        if z_axis(iz) >= target_z
            hold on; plot(target_x, target_z, 'ro', 'MarkerSize', 8); hold off;
        end
        
        % [左下] f-k 域振幅谱 - 观察能量守恒
        subplot(2, 2, 3);
        imagesc(k_x_shift, f_ax_shift, abs(fftshift(P_current_dyn))); colormap(gca, 'jet');
        title('f-k域振幅谱 |P| (保持不变，说明能量守恒)');
        xlabel('k_x (rad/m)'); ylabel('f (MHz)'); ylim([-1500, 1500]);
        
        % [右下] f-k 域相位谱 - 核心！观察相位旋转与频散
        subplot(2, 2, 4);
        imagesc(k_x_shift, f_ax_shift, angle(fftshift(P_current_dyn))); colormap(gca, 'hsv');
        title('f-k域相位谱 \angleP (不同频率的条纹随深度非线性加密)');
        xlabel('k_x (rad/m)'); ylabel('f (MHz)'); ylim([-1500, 1500]);
        
        drawnow; 
    end
end
disp('=== 全部演示完成 ===');