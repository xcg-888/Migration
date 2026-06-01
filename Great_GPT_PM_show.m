clear; clc; close all;

%% ===================== 1. 参数设置 =====================
c = 0.299792458;        % 光速，单位 m/ns
epsilon_r = 9;          % 相对介电常数
v = c / sqrt(epsilon_r);
vm = v / 2;             % 爆炸反射面模型半速度

x0_true = 0;            % 真实目标横向位置 m，即真实双曲线顶点
z0_true = 0.45;         % 真实目标深度 m

Nx = 81;                
x = linspace(-0.8, 0.8, Nx);    
dx = x(2) - x(1);

dt = 0.05;              
t = -25:dt:60;          
Nt = numel(t);

f0 = 0.8;               % Ricker 子波主频，GHz，也即 1/ns

%% ===================== 2. 构造点目标双曲线数据 =====================
data = zeros(Nt, Nx);

for ix = 1:Nx
    R = sqrt((x(ix) - x0_true)^2 + z0_true^2);
    tau = R / vm;       % 爆炸反射面模型单程时间，等效原始双程走时

    tt = t(:) - tau;
    a = pi * f0 * tt;
    wavelet = (1 - 2*a.^2) .* exp(-a.^2);

    data(:, ix) = wavelet;
end

% 加一点背景偏置和弱噪声，便于观察聚焦
data = data + 0.03 * randn(size(data)) + 0.05;

%% ===================== 3. 频率轴与频域数据 =====================
if mod(Nt, 2) == 0
    f = [0:Nt/2-1, -Nt/2:-1] / (Nt * dt);
else
    f = [0:(Nt-1)/2, -(Nt-1)/2:-1] / (Nt * dt);
end

f = f(:);               % Nt x 1

DATA_F = fft(data, [], 1);

%% ===================== 4. 设置两个演示场景 =====================
z_scan = linspace(0.1, 0.8, 120);   % 假设深度扫描范围
[~, t0_idx] = min(abs(t));          % 找到最接近 t = 0 的采样点

% 两个横向成像位置：
% 第一个是顶点位置，第二个是非顶点位置
x_focus_list = [x0_true, 0.30];

case_name = {
    '阶段 1：顶点位置正确补偿';
    '阶段 2：非顶点位置错误补偿'
};

case_short = {
    '顶点位置';
    '非顶点位置'
};

Ncase = numel(x_focus_list);
Nz = numel(z_scan);

image_amp_all = zeros(Ncase, Nz);
stack_record = cell(Ncase, Nz);
comp_record  = cell(Ncase, Nz);

%% ===================== 5. 分别计算顶点与非顶点补偿结果 =====================
for icase = 1:Ncase

    x_focus = x_focus_list(icase);

    for iz = 1:Nz

        z = z_scan(iz);

        % 当前假设成像点为 (x_focus, z)
        % 若 x_focus = x0_true，则对应顶点正确补偿
        % 若 x_focus ~= x0_true，则对应非顶点错误补偿
        tau_pred = sqrt((x - x_focus).^2 + z^2) / vm;   % 1 x Nx

        % 频域相位补偿
        phase = exp(1i * 2*pi * f * tau_pred);          % Nt x Nx

        % 回到时间域
        comp = real(ifft(DATA_F .* phase, [], 1));      % Nt x Nx

        % 所有道叠加
        stack_trace = sum(comp, 2) * dx;                % Nt x 1

        % t = 0 处成像幅值
        image_amp_all(icase, iz) = abs(stack_trace(t0_idx));

        stack_record{icase, iz} = stack_trace;
        comp_record{icase, iz} = comp;
    end
end

%% ===================== 6. 制作融合演示视频：先顶点，后非顶点 =====================
video_name = 'migration_apex_then_non_apex_demo.mp4';

vwriter = VideoWriter(video_name, 'MPEG-4');
vwriter.FrameRate = 12;
open(vwriter);

figure('Color', 'w', 'Position', [100 100 1350 820]);

% 真实双曲线
tau_true = sqrt((x - x0_true).^2 + z0_true^2) / vm;

for icase = 1:Ncase

    x_focus = x_focus_list(icase);

    for iz = 1:Nz

        z = z_scan(iz);

        tau_pred = sqrt((x - x_focus).^2 + z^2) / vm;

        comp = comp_record{icase, iz};
        stack_trace = stack_record{icase, iz};

        % 当前横向成像位置在真实双曲线上的时间
        tau_focus_on_true = sqrt((x_focus - x0_true)^2 + z0_true^2) / vm;

        clf;

        %% ---------- 图 1：原始 B-scan 与当前假设走时曲线 ----------
        subplot(2,2,1);
        imagesc(x, t, data);
        set(gca, 'YDir', 'reverse');
        hold on;

        h1 = plot(x, tau_true, 'c--', 'LineWidth', 1.5);       % 真实双曲线
        h2 = plot(x, tau_pred, 'w', 'LineWidth', 2.2);         % 当前假设曲线

        h3 = plot(x0_true, z0_true/vm, 'ro', ...
            'MarkerSize', 8, 'LineWidth', 2);                 % 真实顶点

        h4 = plot(x_focus, tau_focus_on_true, 'bo', ...
            'MarkerSize', 8, 'LineWidth', 2);                 % 当前横向位置在真实双曲线上的点

        h5 = plot(x_focus, z/vm, 'yo', ...
            'MarkerSize', 7, 'LineWidth', 2);                 % 当前假设成像点对应的曲线顶点

        xline(x0_true, 'r--', 'LineWidth', 1);
        xline(x_focus, 'b--', 'LineWidth', 1);

        xlabel('x / m');
        ylabel('t / ns');

        title([case_name{icase}, ...
               '，x_f = ', num2str(x_focus, '%.2f'), ...
               ' m，z = ', num2str(z, '%.3f'), ' m']);

        legend([h1 h2 h3 h4 h5], ...
            {'真实双曲线', '当前假设走时曲线', '真实顶点', ...
             '当前横向位置对应的真实回波点', '当前假设成像点'}, ...
            'Location', 'southoutside');

        colormap(gray);
        colorbar;
        ylim([0 30]);

        %% ---------- 图 2：相位补偿后的各道信号 ----------
        subplot(2,2,2);
        hold on;

        step = 6;
        scale = 0.12;

        for ix = 1:step:Nx
            plot(t, comp(:, ix) * scale + x(ix), 'k');
        end

        xline(0, 'r--', 'LineWidth', 1.5);

        xlabel('t / ns');
        ylabel('trace position x / m');

        if icase == 1
            title('顶点补偿：正确深度处，各道波包会压到 t = 0');
        else
            title('非顶点补偿：预测走时不匹配，波包难以同时压到 t = 0');
        end

        xlim([-12 12]);
        grid on;

        %% ---------- 图 3：叠加后的波形 ----------
        subplot(2,2,3);
        plot(t, stack_trace, 'k', 'LineWidth', 1.5);
        hold on;

        plot(t(t0_idx), stack_trace(t0_idx), 'ro', ...
            'MarkerSize', 8, 'LineWidth', 2);

        xline(0, 'r--', 'LineWidth', 1.5);

        xlabel('t / ns');
        ylabel('stack amplitude');

        title([case_short{icase}, ...
               '叠加波形：I(z)=|P(t=0)| = ', ...
               num2str(image_amp_all(icase, iz), '%.4f')]);

        xlim([-12 12]);
        grid on;

        %% ---------- 图 4：顶点与非顶点成像幅值对比 ----------
        subplot(2,2,4);
        plot(z_scan, image_amp_all(1,:), 'k', 'LineWidth', 2);
        hold on;
        plot(z_scan, image_amp_all(2,:), 'b--', 'LineWidth', 2);

        xline(z0_true, 'r--', '真实深度', 'LineWidth', 1.5);

        plot(z, image_amp_all(icase, iz), 'ro', ...
            'MarkerSize', 8, 'LineWidth', 2);

        xlabel('假设深度 z / m');
        ylabel('|P(t=0)|');

        title('深度扫描成像结果对比');
        legend('顶点位置补偿', '非顶点位置补偿', 'Location', 'best');
        grid on;

        if icase == 1
            sgtitle('阶段 1：顶点位置正确补偿 → 同相叠加 → 能量聚焦');
        else
            sgtitle('阶段 2：非顶点位置错误补偿 → 走时不匹配 → 聚焦变弱');
        end

        frame = getframe(gcf);
        writeVideo(vwriter, frame);
    end
end

close(vwriter);

disp(['视频已生成：', video_name]);

%% ===================== 7. 静态对比结果图 =====================
figure('Color', 'w', 'Position', [200 150 1250 450]);

%% ---------- 静态图 1：原始 B-scan ----------
subplot(1,2,1);
imagesc(x, t, data);
set(gca, 'YDir', 'reverse');
hold on;

x_non_apex = x_focus_list(2);
tau_non_apex = sqrt((x_non_apex - x0_true)^2 + z0_true^2) / vm;

plot(x, tau_true, 'w', 'LineWidth', 2);
plot(x0_true, z0_true/vm, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot(x_non_apex, tau_non_apex, 'bo', 'MarkerSize', 8, 'LineWidth', 2);

xline(x0_true, 'r--', 'LineWidth', 1);
xline(x_non_apex, 'b--', 'LineWidth', 1);

xlabel('x / m');
ylabel('t / ns');
title('模拟点目标 B-scan：顶点位置与非顶点位置');
legend('真实双曲线', '真实顶点', '非顶点位置', 'Location', 'southoutside');

colormap(gray);
colorbar;
ylim([0 30]);

%% ---------- 静态图 2：顶点与非顶点成像幅值曲线 ----------
subplot(1,2,2);
plot(z_scan, image_amp_all(1,:), 'k', 'LineWidth', 2);
hold on;
plot(z_scan, image_amp_all(2,:), 'b--', 'LineWidth', 2);

xline(z0_true, 'r--', '真实深度', 'LineWidth', 1.5);

% 标出两个曲线最大值
[~, idx_max_apex] = max(image_amp_all(1,:));
[~, idx_max_non_apex] = max(image_amp_all(2,:));

plot(z_scan(idx_max_apex), image_amp_all(1, idx_max_apex), ...
    'ko', 'MarkerSize', 8, 'LineWidth', 2);

plot(z_scan(idx_max_non_apex), image_amp_all(2, idx_max_non_apex), ...
    'bo', 'MarkerSize', 8, 'LineWidth', 2);

xlabel('假设深度 z / m');
ylabel('|P(t=0)|');
title('顶点补偿与非顶点补偿的成像幅值对比');

legend('顶点位置补偿', '非顶点位置补偿', ...
       '真实深度', '顶点曲线峰值', '非顶点曲线峰值', ...
       'Location', 'best');

grid on;