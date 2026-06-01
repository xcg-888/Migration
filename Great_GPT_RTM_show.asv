clear; clc; close all;

%% ===================== 1. 基本参数 =====================
c = 0.299792458;        % 光速，m/ns
epsilon_r = 9;          % 相对介电常数
v = c / sqrt(epsilon_r);
vm = v / 2;             % 爆炸反射面模型速度，m/ns

Nx = 121;               % A-Scan 道数
Nt = 320;               % 时间采样点数

x = linspace(-0.8, 0.8, Nx);     % 测线位置，m
dx = x(2) - x(1);

dt = 0.05;              % 时间采样间隔，ns
t = (0:Nt-1) * dt;      % 时间轴，ns

z = vm * t;             % 等效深度轴，m
dz = z(2) - z(1);

%% ===================== 2. 目标参数 =====================
x0 = 0;                 % 目标横向位置，m
z0 = 0.45;              % 目标深度，m

%% ===================== 3. 构造理想双曲线图像 =====================
data = zeros(Nt, Nx);

sigma_t = 0.08;         % 双曲线线宽，ns，越小越细

for ix = 1:Nx
    r = sqrt((x(ix) - x0)^2 + z0^2);
    ti = r / vm;        % 爆炸反射面模型下的单程等效时间
    data(:, ix) = exp(-((t - ti).^2) / (2 * sigma_t^2));
end

data = data / max(data(:));

% 原图对比度增强
dataShow = data .^ 0.45;
dataShow = min(dataShow * 1.25, 1);

%% ===================== 4. 动画参数 =====================
ampThreshold = 0.10;    % 信号阈值，越小参与传播的点越多
frontWidth = 1;         % 波前厚度，单位：像素
pauseTime = 0.08;       % 动画速度，越大越慢

fieldGamma = 0.45;      % 传播场对比度增强参数，越小越黑越明显
fieldGain = 1.8;        % 传播场显示增益，越大波前越明显

saveGif = false;
gifName = 'reverse_time_circle_migration.gif';

%% ===================== 5. 显示窗口 =====================
figure('Color', 'w', 'Position', [80, 100, 1250, 520]);

% ===================== 左图：原始双曲线 =====================
ax1 = subplot(1, 2, 1);

hLeft = imagesc(x, t, dataShow);
set(ax1, 'YDir', 'reverse');
set(ax1, 'Color', 'w');
colormap(ax1, flipud(gray));    % 0 = 白色，1 = 黑色
caxis(ax1, [0 1]);

xlabel('x / m');
ylabel('t / ns');
title('原始理想双曲线图像');
box on;
hold on;

% 动态横线：表示当前正在逆时传播的那一行
hLine = plot([x(1), x(end)], [t(end), t(end)], ...
    'r-', 'LineWidth', 1.6);

% ===================== 右图：逆时传播场 =====================
ax2 = subplot(1, 2, 2);

hImg = imagesc(x, z, zeros(Nt, Nx));
set(ax2, 'YDir', 'reverse');
set(ax2, 'Color', 'w');
colormap(ax2, flipud(gray));    % 0 = 白色，1 = 黑色
caxis(ax2, [0 1]);

xlabel('x / m');
ylabel('z / m');
title('逆时半圆传播场');
box on;
hold on;

% 小红圆点表示真实目标位置
plot(x0, z0, 'ro', ...
    'MarkerSize', 5, ...
    'LineWidth', 1.2, ...
    'MarkerFaceColor', 'none');

%% ===================== 6. 逆时半圆传播动画 =====================
for frame = 1:Nt
    
    field = zeros(Nt, Nx);
    
    % 当前正在从原图中取出的那一行
    currentRow = Nt - frame + 1;
    set(hLine, 'YData', [t(currentRow), t(currentRow)]);
    
    % 当前帧中，所有已经发射过的行都继续向下传播
    for emitStep = 1:frame
        
        % 从原图最后一行开始，逐行向上逆时取出
        srcRow = Nt - emitStep + 1;
        
        % 该行信号已经传播的时间步数
        ageStep = frame - emitStep;
        
        % 当前半圆传播半径
        radius = vm * ageStep * dt;
        
        % 当前被发射的顶部信号
        srcSignal = data(srcRow, :);
        
        % 只让较强信号参与传播
        srcCols = find(srcSignal > ampThreshold);
        
        for k = 1:numel(srcCols)
            
            srcCol = srcCols(k);
            amp = srcSignal(srcCol);
            xSource = x(srcCol);
            
            % 半圆传播，只向下传播
            horizontalDist = x - xSource;
            inside = abs(horizontalDist) <= radius;
            
            if ~any(inside)
                continue;
            end
            
            xCols = find(inside);
            zFront = sqrt(radius^2 - horizontalDist(inside).^2);
            
            rowFront = round(zFront / dz) + 1;
            
            valid = rowFront >= 1 & rowFront <= Nt;
            rowFront = rowFront(valid);
            xCols = xCols(valid);
            
            % 给波前一点厚度，避免太细看不清
            for w = -frontWidth:frontWidth
                
                rr = rowFront + w;
                valid2 = rr >= 1 & rr <= Nt;
                
                % 波前厚度权重
                weight = exp(-(w^2) / 2);
                
                ind = sub2ind(size(field), rr(valid2), xCols(valid2));
                field(ind) = field(ind) + amp * weight;
            end
        end
    end
    
    %% ===================== 7. 显示增强 =====================
    if max(field(:)) > 0
        showField = field / max(field(:));
    else
        showField = field;
    end
    
    % 增强对比度：让弱波前也更明显
    showField = showField .^ fieldGamma;
    showField = min(showField * fieldGain, 1);
    
    % 直接显示 showField：
    % 因为使用 flipud(gray)，所以 0 是白色，1 是黑色
    set(hImg, 'CData', showField);
    
    subplot(ax2);
    title(['逆时半圆传播场  |  当前逆时行数：', num2str(frame), '/', num2str(Nt)]);
    
    drawnow;
    pause(pauseTime);
    
    %% ===================== 8. 保存 GIF，可选 =====================
    if saveGif
        frameImg = getframe(gcf);
        [A, map] = rgb2ind(frame2im(frameImg), 256);
        
        if frame == 1
            imwrite(A, map, gifName, 'gif', ...
                'LoopCount', inf, ...
                'DelayTime', pauseTime);
        else
            imwrite(A, map, gifName, 'gif', ...
                'WriteMode', 'append', ...
                'DelayTime', pauseTime);
        end
    end
end