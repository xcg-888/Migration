function [Relative_permittivity1, Relative_permittivity2] = extract_lagrange_mean_points(data_out,TrackInterval,dt)
    % ==========================================
    % 【用户配置区】
    % 选择寻找拉格朗日中值点的方法：
    % 'distance' - 寻找矩形框内距离线段（割线）最远的点
    % 'midpoint' - 取线段中点 x，在该类上找对应的点（若存在多个取 y 的均值）
    % ==========================================
    method_choice = 'midpoint'; 
    
    % 物理常数
    c = 3e8; % 光速 (m/s)
    
    % 初始化输出，防止检测不到两条线时报错
    Relative_permittivity1 = NaN;
    Relative_permittivity2 = NaN;
    
    % 确保输入图像为逻辑二值化图像
    data_out = logical(data_out);
    
    % 对输入图像进行1个单位的腐蚀
    se_erode = strel('disk', 1); 
    data_out = imerode(data_out, se_erode);
    
    % ==========================================
    % 1. 自适应霍夫参数计算部分
    % ==========================================
    max_theta = 75; 
    
    stats = regionprops(data_out, 'BoundingBox');
    
    if isempty(stats)
        disp('未检测到任何有效连通域，程序终止。');
        return;
    end
    
    L_all = zeros(length(stats), 1);
    for k_idx = 1:length(stats)
        bb = stats(k_idx).BoundingBox; 
        L_all(k_idx) = sqrt(bb(3)^2 + bb(4)^2);
    end
    
    L_max = max(L_all);
    K = 2; 
    delta_theta_adaptive = (1/K) * rad2deg(atan(2/L_max)); 
    disp(['自适应角度步长 delta_theta_adaptive: ', num2str(delta_theta_adaptive)]);
    
    theta_range = -max_theta:delta_theta_adaptive:max_theta;
    theta_range = theta_range(theta_range >= -90 & theta_range < 90);
    
    % ==========================================
    % 2. 8联通分类与霍夫变换核心部分
    % ==========================================
    [L, n] = bwlabel(data_out, 8);
    fprintf('共检测到 %d 个8联通类。\n', n);
    
    all_peaks_record = []; 
    cell_theta = cell(n, 1);
    cell_rho = cell(n, 1);
    
    cols = ceil(sqrt(n));
    rows = ceil(n / cols);
    
    se_display = strel('disk', 4); 
    
    fig_classes = figure('Name', '每个类的独立图像 (Separated Classes)');
    fig_hough = figure('Name', '每个类的霍夫空间 (Hough Space)');
    
    for i = 1:n
        current_class_img = (L == i);
        
        figure(fig_classes);
        subplot(rows, cols, i);
        
        img_for_display = imdilate(current_class_img, se_display);
        
        % [修改点 1]：imshow 改为 imagesc
        imagesc(~img_for_display); 
        colormap(gca, gray); % 强制使用灰度图映射，确保黑白反转效果不变
        axis normal;          % 锁定真实的图像横纵比
        title(sprintf('Class %d', i));
        axis on;
        
        [H, theta, rho] = hough(current_class_img, 'Theta', theta_range);
        
        cell_theta{i} = theta;
        cell_rho{i} = rho;
        
        P = houghpeaks(H, 2); 
        
        for j = 1:size(P, 1)
            r_idx = P(j, 1); 
            c_idx = P(j, 2); 
            vote_value = H(r_idx, c_idx);
            all_peaks_record = [all_peaks_record; vote_value, r_idx, c_idx, i];
        end
        
        figure(fig_hough);
        subplot(rows, cols, i);
        
        % [修改点 2]：imshow 改为 imagesc，并传入 XData (theta) 和 YData (rho)
        imagesc(theta, rho, 1 - imadjust(rescale(H)));
        colormap(gca, gray); % 强制使用灰度图映射
        title(sprintf('Class %d', i));
        xlabel('\theta'), ylabel('\rho');
        axis on, axis normal, hold on;
        
        if ~isempty(P)
            plot(theta(P(:,2)), rho(P(:,1)), 's', 'color', 'red', 'LineWidth', 1.5);
        end
    end
    
    % ==========================================
    % 3. 全局极值筛选、绘图及拉格朗日中值点计算部分
    % ==========================================
    if isempty(all_peaks_record)
        disp('未检测到任何有效的直线峰值。');
        return;
    end
    
    [~, sort_idx] = sort(all_peaks_record(:, 1), 'descend');
    sorted_peaks = all_peaks_record(sort_idx, :);
    
    num_top = min(2, size(sorted_peaks, 1));
    top_global_peaks = sorted_peaks(1:num_top, :);
    
    figure('Name', sprintf('原图中的全局 Top 2 直线及中值点 (%s)', method_choice));
    
    % [修改点 3]：imshow 改为 imagesc
    imagesc(~data_out); 
    colormap(gca, gray); % 强制使用灰度图映射
    axis normal;          % 锁定原图横纵比
    hold on;
    
    calc_count = 1; 
    
    for k_idx = 1:num_top 
        vote_val = top_global_peaks(k_idx, 1);
        r_idx    = top_global_peaks(k_idx, 2);
        c_idx    = top_global_peaks(k_idx, 3);
        class_id = top_global_peaks(k_idx, 4);
        
        fprintf('\n=======================================================\n');
        fprintf('全局第 %d 大直线: 属于类 %d, 霍夫投票值为 %d\n', k_idx, class_id, vote_val);
        
        target_class_img = (L == class_id);
        theta = cell_theta{class_id};
        rho = cell_rho{class_id};
        
        % 获取当前类所有像素坐标
        [Y_class, X_class] = find(target_class_img);
        
        % 寻找该类的最高点（即图像中 Y 最小的点）求 x0
        min_Y = min(Y_class);
        idx_min_Y = (Y_class == min_Y);
        x0 = mean(X_class(idx_min_Y)); 
        
        lines = houghlines(target_class_img, theta, rho, [r_idx, c_idx], 'FillGap', 20, 'MinLength', 10);
        
        for line_idx = 1:length(lines)
            % 如果已经计算了两个介电常数，跳出循环
            if calc_count > 2
                break;
            end
            
            xy = [lines(line_idx).point1; lines(line_idx).point2];
            x1 = xy(1,1); y1 = xy(1,2);
            x2 = xy(2,1); y2 = xy(2,2);
            
            plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'blue');
            plot(x1, y1, 'x', 'LineWidth', 2, 'Color', 'red'); 
            plot(x2, y2, 'x', 'LineWidth', 2, 'Color', 'red');  
            
            % 提取线段斜率 |k|，也就是公式中的 f
            theta_line = lines(line_idx).theta; 
            k_slope = abs(cotd(theta_line)); 
            f = k_slope;
            
            lagrange_found = false; 
            
            % ========================================================
            % 根据用户选择的方法寻找中值点
            % ========================================================
            if strcmp(method_choice, 'distance')
                xmin = min(x1, x2); xmax = max(x1, x2);
                ymin = min(y1, y2); ymax = max(y1, y2);
                
                idx_in_box = (X_class >= xmin) & (X_class <= xmax) & ...
                             (Y_class >= ymin) & (Y_class <= ymax);
                X_box = X_class(idx_in_box);
                Y_box = Y_class(idx_in_box);
                
                if ~isempty(X_box)
                    A = y2 - y1;
                    B = x1 - x2;
                    C = x2*y1 - x1*y2;
                    denominator = sqrt(A^2 + B^2);
                    
                    if denominator > 0
                        distances = abs(A .* X_box + B .* Y_box + C) / denominator;
                        [~, max_idx] = max(distances);
                        lagrange_x = X_box(max_idx);
                        lagrange_y = Y_box(max_idx);
                        lagrange_found = true;
                        
                        plot(lagrange_x, lagrange_y, 'o', 'MarkerFaceColor', 'g', ...
                             'MarkerEdgeColor', 'k', 'MarkerSize', 6);
                        
                        fprintf('  -> [距离法] 线段 %d 端点: (%d,%d)-(%d,%d)\n', line_idx, x1, y1, x2, y2);
                    end
                end
                
            elseif strcmp(method_choice, 'midpoint')
                x_mid = round((x1 + x2) / 2);
                idx_mid = (X_class == x_mid);
                Y_at_xmid = Y_class(idx_mid);
                
                if ~isempty(Y_at_xmid)
                    lagrange_x = x_mid;
                    lagrange_y = round(mean(Y_at_xmid)); 
                    lagrange_found = true;
                    
                    plot(lagrange_x, lagrange_y, 'o', 'MarkerFaceColor', 'y', ...
                         'MarkerEdgeColor', 'k', 'MarkerSize', 6);
                    
                    fprintf('  -> [中点法] 线段 %d 端点: (%d,%d)-(%d,%d)\n', line_idx, x1, y1, x2, y2);
                else
                    fprintf('  -> [警告] 在中点 X=%d 处未找到该类的像素点。\n', x_mid);
                end
            end
            
            % ========================================================
            % 计算并分配相对介电常数 
            % ========================================================
            if lagrange_found
                x = lagrange_x;
                t = lagrange_y;
                
                % 【新增】：在图上标出该类的最高点 (顶点)，用品红色三角形标记
                plot(x0, min_Y, '^', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
                
                fprintf('     斜率 |f|: %.4f\n', f);
                fprintf('     类的最高点坐标 (x0, y0): (%.2f, %d)\n', x0, min_Y);
                fprintf('     拉格朗日点 (x, t): (%d, %d)\n', x, t);
                
                if abs(x - x0) > 0
                    eps_r = (t * dt * dt * c * c * f) / (4 * TrackInterval * TrackInterval * abs(x - x0));
                    
                    % 根据计数器分配给对应的输出变量
                    if calc_count == 1
                        Relative_permittivity1 = eps_r;
                        fprintf('     >>> 计算所得 Relative_permittivity1: %.4f <<<\n', eps_r);
                    elseif calc_count == 2
                        Relative_permittivity2 = eps_r;
                        fprintf('     >>> 计算所得 Relative_permittivity2: %.4f <<<\n', eps_r);
                    end
                    calc_count = calc_count + 1;
                else
                    fprintf('     [错误] |x - x0| 为 0，无法计算该线段的介电常数。\n');
                end
            end
            
        end
    end
    axis on;
    hold off;
    fprintf('\n=======================================================\n');
end