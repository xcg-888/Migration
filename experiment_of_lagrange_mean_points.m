clear; close all; clc;

% 1. 基本参数设置
% 提取下半分支方程 (代码中 y 代表 t)
hyberbola_func = @(x) -2 - 2 * sqrt(36 + (x - 5).^2);

% 2. 生成无规律的随机端点
% 固定两端和顶点(x=5)，在中间随机插入 4 个点，并进行排序保证单调递增
% 左侧区间 [-10, 5] 随机划分
x_left_random = sort(-10 + 15 * rand(1, 4)); 
x_left = [-10, x_left_random, 5];    

% 右侧区间 [5, 20] 随机划分
x_right_random = sort(5 + 15 * rand(1, 4));
x_right = [5, x_right_random, 20];   

y_left = hyberbola_func(x_left);
y_right = hyberbola_func(x_right);

% 3. 计算无规律割线的斜率 f
slopes_left = diff(y_left) ./ diff(x_left);
slopes_right = diff(y_right) ./ diff(x_right);

% --- 画图准备 ---
x_full = linspace(-15, 25, 500); 
y_full = hyberbola_func(x_full);

figure('Color', 'w', 'Name', '无规律割线反演渐近线斜率');
plot(x_full, y_full, 'Color', [0.6 0.8 1], 'LineWidth', 3); hold on;
% 画出随机割线
plot(x_left, y_left, 'ro-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
plot(x_right, y_right, 'ro-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

% 标注顶点 (x=5, y=-14)
plot(5, -14, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
text(5, -13, '顶点 (5, -14)', 'HorizontalAlignment', 'center');

% 4 & 5. 打印数据，寻找中值点进行反演计算
fprintf('================= 左侧 5 条无规律割线反演计算 =================\n');
for i = 1:5
    x1 = x_left(i);   y1 = y_left(i);
    x2 = x_left(i+1);
    f = slopes_left(i);
    
    % 最大距离法求拉格朗日中值点
    dist_func = @(x) -(hyberbola_func(x) - (y1 + f * (x - x1)));
    x_r = fminbnd(dist_func, x1, x2);
    y_r = hyberbola_func(x_r); % y_r 就是此时的 t_r
    
    % 使用反演公式计算渐近线斜率 k: sqrt(abs( f * (t_r + 2) / (x_r - 5) ))
    k_inverse = sqrt(abs(f * (y_r + 2) / (x_r - 5)));
    
    % 打印区间的跨度(以显示其无规律性)和反演结果
    fprintf('区间跨度: %5.2f | f = %7.4f | 反演渐近线斜率 k = %7.4f\n', ...
        x2 - x1, f, k_inverse);
        
    plot(x_r, y_r, 'g*', 'MarkerSize', 8, 'LineWidth', 1.5);
end

fprintf('\n================= 右侧 5 条无规律割线反演计算 =================\n');
for i = 1:5
    x1 = x_right(i);   y1 = y_right(i);
    x2 = x_right(i+1);
    f = slopes_right(i);
    
    dist_func = @(x) -(hyberbola_func(x) - (y1 + f * (x - x1)));
    x_r = fminbnd(dist_func, x1, x2);
    y_r = hyberbola_func(x_r);
    
    k_inverse = sqrt(abs(f * (y_r + 2) / (x_r - 5)));
    
    fprintf('区间跨度: %5.2f | f = %7.4f | 反演渐近线斜率 k = %7.4f\n', ...
        x2 - x1, f, k_inverse);
        
    plot(x_r, y_r, 'g*', 'MarkerSize', 8, 'LineWidth', 1.5);
end

% 图形美化
title('使用随机无规律割线验证拉格朗日中值点并反演渐近线');
xlabel('X 轴');
ylabel('T 轴');
axis([-15 25 -45 0]);
grid on;