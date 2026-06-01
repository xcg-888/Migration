clear; clc; close all;

%% ===================== 1. 基本参数 =====================
c = 0.299792458;        % 光速，m/ns
epsilon_r = 9;          % 相对介电常数
v = c / sqrt(epsilon_r);
vm = v / 2;             % 爆炸反射面模型速度，m/ns

Nx = 121;               % A-Scan 道数
Nt = 320;               % 时间采样点数

x = linspace(-0.8, 0.8, Nx);     % 测线位置，m
dt = 0.05;              % 时间采样间隔，ns
t = (0:Nt-1) * dt;      % 时间轴，ns

z = vm * t;             % 等效深度轴，m
dz = z(2) - z(1);

%% ===================== 2. 目标参数 =====================
x0 = 0;                 % 目标横向位置，m
z0 = 0.45;              % 目标深度，m

%% ===================== 3. 构造理想双曲线图像 =====================
data = zeros(Nt, Nx);
sigma_t = 0.08;         % 双曲线线宽，ns

ti_all = zeros(1, Nx);  % 每一道对应的双曲线峰值时刻

for ix = 1:Nx
    r = sqrt((x(ix) - x0)^2 + z0^2);
    ti = r / vm;        % 爆炸反射面模型下的单程等效时间
    ti_all(ix) = ti;
    data(:, ix) = exp(-((t - ti).^2) / (2 * sigma_t^2));
end

data = data / max(data(:));

% 原图显示增强
dataShow = data .^ 0.45;
dataShow = min(dataShow * 1.25, 1);

%% ===================== 4. 自动选取 7 个对称代表点 =====================
% 顶点道
[~, apexIdx] = min(abs(x - x0));

% 只考虑图像中“可见”的双曲线部分
visibleIdx = find(ti_all <= t(end));

% 左右两翼可见部分
leftVisible  = visibleIdx(visibleIdx < apexIdx);
rightVisible = visibleIdx(visibleIdx > apexIdx);

% 左右末尾点：双曲线在图中最靠两边的可见点
leftTailIdx  = leftVisible(1);
rightTailIdx = rightVisible(end);

% 在左翼“末尾 -> 顶点”之间再分成 3 段，取 2 个中间点
leftOuterIdx = round(leftTailIdx + (apexIdx - leftTailIdx) / 3);
leftInnerIdx = round(leftTailIdx + 2 * (apexIdx - leftTailIdx) / 3);

% 在右翼“顶点 -> 末尾”之间再分成 3 段，取 2 个中间点
rightInnerIdx = round(apexIdx + (rightTailIdx - apexIdx) / 3);
rightOuterIdx = round(apexIdx + 2 * (rightTailIdx - apexIdx) / 3);

% 最终 7 个点
selectedIdx = [leftTailIdx, leftOuterIdx, leftInnerIdx, ...
               apexIdx, ...
               rightInnerIdx, rightOuterIdx, rightTailIdx];

selectedIdx = unique(selectedIdx, 'stable');   % 防止重复
selectedX = x(selectedIdx);
selectedT = ti_all(selectedIdx);

%% ===================== 5. 动画参数 =====================
ampThreshold = 0.08;    % 信号阈值
frontWidth = 2;         % 波前厚度（像素）
pauseTime = 0.10;       % 动画速度

fieldGamma = 0.45;      % 波前对比度增强
fieldGain = 1.8;        % 波前显示增益

saveGif = false;
gifName = 'reverse_time_circle_selected7.gif';

%% ===================== 6. 显示窗口 =====================
figure('Color', 'w', 'Position', [80, 100, 1320, 540]);

% ===================== 左图：原始双曲线 =====================
ax1 = subplot(1, 2, 1);

imagesc(x, t, dataShow);
set(ax1, 'YDir', 'reverse');
set(ax1, 'Color', 'w');
colormap(ax1, flipud(gray(256)));   % 0=白, 1=黑
caxis(ax1, [0 1]);

xlabel('x / m');
ylabel('t / ns');
title('原始理想双曲线图像（仅选7个代表性点）');
box on;
hold on;

% 动态横线：表示当前正在逆时传播的那一行
hLine = plot([x(1), x(end)], [t(end), t(end)], ...
    'r-', 'LineWidth', 1.5);

% 标出 7 个选中的代表性点
plot(selectedX, selectedT, 'ro', ...
    'MarkerSize', 7, ...
    'LineWidth', 1.4, ...
    'MarkerFaceColor', 'none');

% 顶点单独强调
plot(x(apexIdx), ti_all(apexIdx), 'rs', ...
    'MarkerSize', 8, ...
    'LineWidth', 1.4);

% 给 7 个点加标签
labels = {'左末尾', '左外中点', '左内中点', '顶点', ...
          '右内中点', '右外中点', '右末尾'};

for k = 1:numel(selectedIdx)
    text(selectedX(k)+0.012, selectedT(k)-0.22, labels{k}, ...
        'Color', 'r', 'FontSize', 10);
end

% ===================== 右图：逆时传播场 =====================
ax2 = subplot(1, 2, 2);

hImg = imagesc(x, z, zeros(Nt, Nx));
set(ax2, 'YDir', 'reverse');
set(ax2, 'Color', 'w');
colormap(ax2, flipud(gray(256)));   % 0=白, 1=黑
caxis(ax2, [0 1]);

xlabel('x / m');
ylabel('z / m');
title('7个代表性点的逆时半圆传播场');
box on;
hold on;

% 真实目标位置
plot(x0, z0, 'ro', ...
    'MarkerSize', 5, ...
    'LineWidth', 1.2, ...
    'MarkerFaceColor', 'none');

%% ===================== 7. 逆时半圆传播动画 =====================
for frame = 1:Nt
    
    field = zeros(Nt, Nx);
    
    % 左图动态横线：当前处理到的时间行
    currentRow = Nt - frame + 1;
    set(hLine, 'YData', [t(currentRow), t(currentRow)]);
    
    % 当前帧中，所有已经发射过的时间行继续传播
    for emitStep = 1:frame
        
        % 从原图最后一行开始逆时取出
        srcRow = Nt - emitStep + 1;
        
        % 已传播的时间步
        ageStep = frame - emitStep;
        
        % 当前传播半径
        radius = vm * ageStep * dt;
        
        % 当前这一行信号
        srcSignal = data(srcRow, :);
        
        % 只保留 7 个选中点，并要求该行幅值大于阈值
        srcCols = selectedIdx(srcSignal(selectedIdx) > ampThreshold);
        
        for k = 1:numel(srcCols)
            
            srcCol = srcCols(k);
            amp = srcSignal(srcCol);
            xSource = x(srcCol);
            
            % 半圆传播，只向下
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
            
            % 给波前一点厚度
            for w = -frontWidth:frontWidth
                rr = rowFront + w;
                valid2 = rr >= 1 & rr <= Nt;
                weight = exp(-(w^2) / 2);
                
                ind = sub2ind(size(field), rr(valid2), xCols(valid2));
                field(ind) = field(ind) + amp * weight;
            end
        end
    end
    
    %% ===================== 8. 显示增强 =====================
    if max(field(:)) > 0
        showField = field / max(field(:));
    else
        showField = field;
    end
    
    showField = showField .^ fieldGamma;
    showField = min(showField * fieldGain, 1);
    
    set(hImg, 'CData', showField);
    
    title(ax2, ['7个代表性点的逆时半圆传播场  |  当前逆时行数：', ...
        num2str(frame), '/', num2str(Nt)]);
    
    drawnow;
    pause(pauseTime);
    
    %% ===================== 9. 保存 GIF，可选 =====================
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

%% ===================== 10. 输出选中结果 =====================
disp('选中的7个代表性道索引为：');
disp(selectedIdx);

disp('对应的 x 位置 (m) 为：');
disp(selectedX);

disp('对应的双曲线峰值时刻 t (ns) 为：');
disp(selectedT);