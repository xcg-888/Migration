function metrics_table = compareMigrationMetrics(images_cell, names_cell)
    num_algs = length(images_cell);
    SCR_list = zeros(num_algs, 1);
    ISLR_list = zeros(num_algs, 1);
    PSLR_list = zeros(num_algs, 1);

    % 定义空间保护窗的半宽 (像素数)
    % 意思是：以最高点为中心，上下左右各延伸 15 个像素的方形区域都视为主瓣区
    % 注意：如果您的图像分辨率很高/低，可以适当调大/调小这个值 (例如 10~30)
    guard_win = 15; 

    for i = 1:num_algs
        img = images_cell{i};
        img(isnan(img)) = 0;
        
        if ~isreal(img)
            I = abs(img).^2;
        else
            I = abs(img).^2; 
        end
        
        P_tot = sum(I(:));
        if P_tot == 0
            SCR_list(i) = NaN; ISLR_list(i) = NaN; PSLR_list(i) = NaN; continue;
        end
        
        % 寻找峰值及其坐标
        [max_val, max_idx] = max(I(:));
        [r_max, c_max] = ind2sub(size(I), max_idx);
        
        % 划分 -3dB 区域计算 SCR 和 ISLR
        threshold_3dB = max_val * 0.5;
        main_lobe_mask = I >= threshold_3dB;
        P_tar = sum(I(main_lobe_mask));
        P_clutter = P_tot - P_tar;
        if P_clutter <= 0, P_clutter = eps; end
        
        SCR_list(i) = 10 * log10(P_tar / P_clutter);
        ISLR_list(i) = 10 * log10(P_tar / P_clutter); 
        
        % --- 修复后的 PSLR 计算 ---
        I_sidelobes = I;
        
        % 计算保护窗边界 (防止数组越界)
        r_start = max(1, r_max - guard_win);
        r_end   = min(size(I, 1), r_max + guard_win);
        c_start = max(1, c_max - guard_win);
        c_end   = min(size(I, 2), c_max + guard_win);
        
        % 将空间保护窗内的能量清零 (彻底挖去主瓣)
        I_sidelobes(r_start:r_end, c_start:c_end) = 0; 
        
        peak_sidelobe = max(I_sidelobes(:));
        
        if peak_sidelobe <= 0
            PSLR_list(i) = NaN;
        else
            PSLR_list(i) = 10 * log10(peak_sidelobe / max_val);
        end
    end

    metrics_table = table(SCR_list, ISLR_list, PSLR_list, ...
        'RowNames', names_cell, 'VariableNames', {'SCR_dB', 'ISLR_dB', 'PSLR_dB'});
        
    disp('========================================================');
    disp('                  偏移成像质量指标对比                  ');
    disp('========================================================');
    disp(metrics_table);
end