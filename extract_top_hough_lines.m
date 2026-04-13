function [Relative_permittivity,alpha] = extract_top_hough_lines(data_out,TrackInterval,dt)
    % 确保输入图像为逻辑二值化图像
    fprintf('\n=======================================================\n');
    data_out = logical(data_out);
    
    % ==========================================
    % 1. 自适应霍夫参数计算部分
    % ==========================================
    max_theta = 75; % 霍夫变换最大角度
    
    % 寻找连通域属性 (去除了 'Area' 的提取，因为不再需要面积过滤)
    stats = regionprops(data_out, 'BoundingBox');
    
    if isempty(stats)
        disp('未检测到任何有效连通域，程序终止。');
        return;
    end
    
    % 计算每个类外接矩形的对角线长度
    L_all = zeros(length(stats), 1);
    for k = 1:length(stats)
        bb = stats(k).BoundingBox; 
        L_all(k) = sqrt(bb(3)^2 + bb(4)^2);
    end
    
    % 提取最大对角线并计算自适应步长
    L_max = max(L_all);
    K = 2; % K倍过采样
    delta_theta_adaptive = (1/K) * rad2deg(atan(2/L_max)); 
    disp(['自适应角度步长 delta_theta_adaptive: ', num2str(delta_theta_adaptive)]);
    
    % 预先计算并限制全局 theta 范围，加入 [-90, 90) 的越界保护
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
    
    % 计算 subplot 的布局行数和列数
    cols = ceil(sqrt(n));
    rows = ceil(n / cols);
    
    % 【新增】：创建一个半径为2的圆盘结构元素，用于展示时的图像膨胀（使线宽变粗约2-4个像素）
    se_display = strel('disk', 4); 
    
    % 创建两个独立的图形窗口，并保存句柄
    fig_classes = figure('Name', '每个类的独立图像 (Separated Classes)');
    fig_hough = figure('Name', '每个类的霍夫空间 (Hough Space)');
    
    for i = 1:n
        current_class_img = (L == i);
        
        % 切换到“独立图像”窗口，绘制当前类
        figure(fig_classes);
        subplot(rows, cols, i);
        
        % 【修改】：为了可视化效果，专门生成膨胀后的图像用于展示
        img_for_display = imdilate(current_class_img, se_display);
        imshow(~img_for_display); % 反色展示膨胀后的图像，保持白底黑图风格
        
        title(sprintf('Class %d', i));
        axis on;
        
        % --- 计算霍夫变换（注意：这里依然使用的是未膨胀的原图 current_class_img，保证算法精度） ---
        [H, theta, rho] = hough(current_class_img, 'Theta', theta_range);
        [r , c] = size(H); 
        % H = [zeros(1,c);H(1:r-1,:)]+[H(2:r,:);zeros(1,c)]+[zeros(2,c);H(1:r-2,:)]+[H(3:r,:);zeros(2,c)]+H;
        H = [zeros(1,c);H(1:r-1,:)]+[H(2:r,:);zeros(1,c)]+H;
        cell_theta{i} = theta;
        cell_rho{i} = rho;
        
        P = houghpeaks(H, 2); 
        
        for j = 1:size(P, 1)
            r_idx = P(j, 1); % ρ
            c_idx = P(j, 2); % θ
            vote_value = H(r_idx, c_idx);
            all_peaks_record = [all_peaks_record; vote_value, r_idx, c_idx, i];
        end
        
        % 切换到“霍夫空间”窗口，绘制霍夫图
        figure(fig_hough);
        subplot(rows, cols, i);
        imshow(1 - imadjust(rescale(H)), 'XData', theta, 'YData', rho, ...
               'InitialMagnification', 'fit');
        title(sprintf('Class %d', i));
        xlabel('\theta'), ylabel('\rho');
        axis on, axis normal, hold on;
        
        if ~isempty(P)
            plot(theta(P(:,2)), rho(P(:,1)), 's', 'color', 'red', 'LineWidth', 1.5);
        end
    end
    
    % ==========================================
    % 3. 全局极值筛选与绘图部分
    % ==========================================
    if isempty(all_peaks_record)
        disp('未检测到任何有效的直线峰值。');
        return;
    end
    
    [~, sort_idx] = sort(all_peaks_record(:, 1), 'descend');
    sorted_peaks = all_peaks_record(sort_idx, :);
    
    num_top = min(2, size(sorted_peaks, 1));
    top_global_peaks = sorted_peaks(1:num_top, :);
    
    figure('Name', '原图中的全局 Top 2 直线');
    imagesc(~data_out); 
    colormap(gca, gray); % 强制使用灰度图映射
    axis normal;          % 锁定原图横纵比
    hold on;
    
    for k = 1:num_top
        vote_val = top_global_peaks(k, 1);
        r_idx    = top_global_peaks(k, 2);
        c_idx    = top_global_peaks(k, 3);
        class_id = top_global_peaks(k, 4);
        
        fprintf('全局第 %d 大直线: 属于类 %d, 霍夫投票值为 %d\n', k, class_id, vote_val);
        
        target_class_img = (L == class_id);
        theta = cell_theta{class_id};
        rho = cell_rho{class_id};
        
        lines = houghlines(target_class_img, theta, rho, [r_idx, c_idx], 'FillGap', 20, 'MinLength', 10);
        
        for line_idx = 1:length(lines)
            xy = [lines(line_idx).point1; lines(line_idx).point2];
            plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'blue');
            plot(xy(1,1), xy(1,2), 'x', 'LineWidth', 2, 'Color', 'red'); 
            plot(xy(2,1), xy(2,2), 'x', 'LineWidth', 2, 'Color', 'red');  
        end
    end
    axis on;
    hold off;
    

    beta = 1/0.7773; % 1/α系数
    % beta = 1/0.82;
    eclipse_line = 10.5; % 假想εr

    k_line = -cot(theta(top_global_peaks(:, 3))/180*pi);
    disp('data_out');
    % disp(['角度：' num2str(theta(top_global_peaks(:, 3))/180*pi)]);
    disp(['角度：' num2str(theta(top_global_peaks(:, 3)))]);
    disp(['斜率：' num2str(k_line)]);
    c = 3e8;
    % Relative_permittivity = mean((beta*abs(k_line1)*dt*Downsample_N*c/2/TrackInterval).^2);
    Relative_permittivity = mean(beta*abs(k_line)*dt*c/2/TrackInterval).^2;
    alpha = mean(abs(k_line))*dt*c/2/TrackInterval/sqrt(eclipse_line);

    fprintf('\n=======================================================\n');

end