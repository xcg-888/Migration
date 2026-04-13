function [B_scan_final] = Direct_Coupled_Wave_Filtering(B_scan_image_edge,mode)
%% Direct_Coupled_Wave_Filtering - 探地雷达(GPR)二值图像形态学滤波
% 作者:许财国
% 功能描述:
%  “先串联、后并联”的混合滤波架构：
%      针对探地雷达 B-scan 边缘检测后的二值图像，利用连通域的几何特征剔除杂乱干扰波。
%   采用- 第一关 (串联)：全局面积过滤，预先滤除大量细碎噪点，降低后续计算量并防止基数污染。
%     - 第二关 (并联)：在剩余的较大目标中，分别计算“外接矩形斜率”与“基于顶点的左右对称度”。
%                    若目标在这两个指标中任意一个排名靠后（即满足任一剔除条件），则被滤除。
%   特别地，对称度计算采用了“顶点对齐截断”策略，并加入了顶点边缘防线，能有效识别并保留
%   双曲线（包括单侧长拖尾的双曲线），同时筛除斜向或单调的杂乱反射波。
%
% 输入参数:
%   B_scan_image_edge - (Matrix) 经过边缘检测（如Canny）及膨胀闭合处理后的逻辑二值图像。
%         mode        - (Matrix) 输入图像类型
% 输出参数:
%   B_scan_final      - (Matrix) 经过三重形态学特征筛选后，保留下来的有效目标二值图像。
% ==============================================================================

if strcmp(mode,'out') 
    % ==================== 可调参数设置（仿真数据） ====================
    area_drop_ratio = 0.00;     % 面积剔除比例：第一关，剔除全局连通域中面积最小的后 75%
    slope_drop_ratio = 0.00;    % 斜率剔除比例：第二关，剔除剩余目标中外接矩形高宽比（斜率）最小的后 50%
    symmetry_drop_ratio = 0.00; % 对称度剔除比例：第二关，剔除剩余目标中左右对称度最差的后 50%
    apex_edge_ratio = 0.00;     % 顶点边缘剔除阈值：若最高点（顶点）所处位置在当前目标宽度的前 20% 或后 20%，视为杂波
    % ================================================================
else
    % ==================== 可调参数设置（实测数据） ====================
    area_drop_ratio = 0.60;     % 面积剔除比例：第一关，剔除全局连通域中面积最小的后 75%
    slope_drop_ratio = 0.50;    % 斜率剔除比例：第二关，剔除剩余目标中外接矩形高宽比（斜率）最小的后 50%
    symmetry_drop_ratio = 0.50; % 对称度剔除比例：第二关，剔除剩余目标中左右对称度最差的后 50%
    apex_edge_ratio = 0.20;     % 顶点边缘剔除阈值：若最高点（顶点）所处位置在当前目标宽度的前 20% 或后 20%，视为杂波
    % ================================================================
end

% 1. 确保输入图像是逻辑类型 (0和1)，并获取图像尺寸用于后续坐标映射
B_scan_image_edge = logical(B_scan_image_edge);
[img_H, img_W] = size(B_scan_image_edge);

% 2. 对二值图像进行 8 连通聚类，将相邻的像素判定为一个目标
CC = bwconncomp(B_scan_image_edge, 8);
% 提取连通域属性：BoundingBox(外接矩形,用于算斜率), PixelIdxList(像素绝对索引), Area(总像素数)
stats = regionprops(CC, 'BoundingBox', 'PixelIdxList', 'Area');
num_clusters = length(stats);

% 初始化最终滤波图像（全黑背景）
B_scan_final = false(img_H, img_W);

if num_clusters == 0
    warning('图像中未找到任何连通区域。');
    return;
else
    % ==================== 第一关：串联面积滤波 ====================
    % 目的：快速去除零星噪点
    areas = [stats.Area];
    [~, sort_area_idx] = sort(areas, 'ascend'); % 按面积从小到大排序
    
    % 计算需要剔除的小面积目标数量
    drop_count_area = floor(num_clusters * area_drop_ratio);
    % 保存通过第一关的大面积目标索引
    keep_idx_area = sort_area_idx(drop_count_area + 1 : end); 
    num_stage1_kept = length(keep_idx_area);
    
    if num_stage1_kept == 0
        warning('所有目标均被面积滤波剔除。');
        return;
    end
    
    % ==================== 第二关：并联斜率与基于顶点的对称度滤波 ====================
    % 初始化第二关的特征数组（长度为通过第一关的目标数量）
    slopes_stage1 = zeros(1, num_stage1_kept);
    syms_stage1   = zeros(1, num_stage1_kept);
    
    for i = 1:num_stage1_kept
        global_idx = keep_idx_area(i); % 映射回全局索引
        pixel_idx = stats(global_idx).PixelIdxList;
        
        % --- 计算斜率特征 ---
        bbox = stats(global_idx).BoundingBox;
        width = bbox(3);
        height = bbox(4);
        if width > 0
            slopes_stage1(i) = height / width; % 用外接矩形的高宽比近似斜率
        else
            slopes_stage1(i) = inf;
        end
        
        % --- 计算对称度特征 (基于顶点的截断交并比算法) ---
        % 将一维线性索引还原为二维的行(r,快时间轴)和列(c,慢时间轴)坐标
        [r, c] = ind2sub([img_H, img_W], pixel_idx);
        
        % 1. 寻找双曲线顶点 (GPR图像中时间最小处，即行坐标r最小的像素)
        min_r = min(r);
        apex_c_list = c(r == min_r);
        c_axis = round(mean(apex_c_list)); % 如果有多个最高点，取其平均位置作为对称轴
        
        % 2. 计算当前目标的水平跨度及边界
        min_c = min(c);
        max_c = max(c);
        bbox_width = max_c - min_c + 1;
        
        % 3. 新增防线：判断最高点是否处于边缘危险区
        % 计算顶点在当前目标宽度内的相对位置比例
        relative_pos = (c_axis - min_c) / bbox_width;
        
        if relative_pos < apex_edge_ratio || relative_pos > (1 - apex_edge_ratio)
            % 顶点太靠两端 (属于单调斜线或被严重截断的残波)，直接赋予最低对称度 0.0
            syms_stage1(i) = 0.0;
        else
            % 4. 顶点位置合格，执行截断式对称度计算
            % 解决双曲线两翼不等长问题：切掉长尾巴，按较短一侧界定评估窗口
            dist_left = c_axis - min_c;
            dist_right = max_c - c_axis;
            R = min(dist_left, dist_right); % 有效对称评估半径
            
            if R <= 0
                % 半径极小的情况处理
                if bbox_width <= 2 
                    syms_stage1(i) = 1.0; % 极窄的垂直线，天然对称
                else
                    syms_stage1(i) = 0.0; % 纯水平/斜向线，不对称
                end
            else
                % 提取有效评估半径 [c_axis - R, c_axis + R] 内的像素
                valid_mask = (c >= c_axis - R) & (c <= c_axis + R);
                r_sym = r(valid_mask);
                c_sym = c(valid_mask);
                
                % 建立一个局部完美的居中画布，用于交并比(IoU)计算
                min_r_sym = min(r_sym);
                h_sym = max(r_sym) - min_r_sym + 1;
                w_sym = 2 * R + 1; % 画布宽度严格约束，使 c_axis 绝对居中
                
                local_sym_img = false(h_sym, w_sym);
                local_r = r_sym - min_r_sym + 1;
                local_c = c_sym - (c_axis - R) + 1; % 相对坐标平移映射
                local_idx = sub2ind([h_sym, w_sym], local_r, local_c);
                local_sym_img(local_idx) = true;
                
                % 计算图像的左右翻转交并比 (Intersection over Union)
                flipped_img = fliplr(local_sym_img); % 左右镜像
                intersection = sum(local_sym_img & flipped_img, 'all'); % 重叠像素数 (交集)
                union_pixels = sum(local_sym_img | flipped_img, 'all'); % 总占据像素数 (并集)
                
                if union_pixels > 0
                    syms_stage1(i) = intersection / union_pixels; % 对称度得分 [0, 1]
                else
                    syms_stage1(i) = 0;
                end
            end
        end
    end
    
    % ==================== 联合筛选：生成并联黑名单 ====================
    % 在第二关内部建立局部索引
    local_idx_all = 1:num_stage1_kept;
    
    % 分支 A：寻找斜率表现最差的批次
    [~, sort_slope_idx_local] = sort(slopes_stage1, 'ascend');
    drop_count_slope = floor(num_stage1_kept * slope_drop_ratio);
    drop_idx_slope_local = sort_slope_idx_local(1:drop_count_slope);
     
    % 分支 B：寻找对称度表现最差的批次
    [~, sort_sym_idx_local] = sort(syms_stage1, 'ascend');
    drop_count_sym = floor(num_stage1_kept * symmetry_drop_ratio);
    drop_idx_sym_local = sort_sym_idx_local(1:drop_count_sym);
      
    % 并联判定：取黑名单的并集（即满足任一剔除条件，就被无情淘汰）
    drop_idx_stage2_local = unique([drop_idx_slope_local, drop_idx_sym_local]);
    
    % 找出最终幸存的目标局部索引
    keep_idx_final_local = setdiff(local_idx_all, drop_idx_stage2_local);
    % 映射回原图的全局连通域索引
    final_keep_idx = keep_idx_area(keep_idx_final_local);
    
    % ==================== 生成最终图像 ====================
    for i = 1:length(final_keep_idx)
        B_scan_final(stats(final_keep_idx(i)).PixelIdxList) = true;
    end
    
    % ==================== 画图对比展示 ====================
    % 绘制 1x5 的过程展示图，直观对比各个过滤分支的效果
    figure('Name', '串并联混合滤波结果 (含顶点边缘剔除)', 'Position', [50, 100, 1800, 350]);
    
    % 1. 原图全貌
    subplot(1, 5, 1);
    imshow(~B_scan_image_edge);
    title(sprintf('1. 原始图像\n(共 %d 个)', num_clusters));
    hold on; for i=1:num_clusters, rectangle('Position', stats(i).BoundingBox, 'EdgeColor', '[0.8 0.8 0.8]', 'LineWidth', 0.5); end; hold off;
    
    % 2. 面积筛选结果
    subplot(1, 5, 2);
    img_area = false(img_H, img_W);
    for i=1:length(keep_idx_area), img_area(stats(keep_idx_area(i)).PixelIdxList) = true; end
    imshow(~img_area);
    title(sprintf('2. 面积达标\n(剩 %d 个)', num_stage1_kept));
    hold on; for i=1:length(keep_idx_area), rectangle('Position', stats(keep_idx_area(i)).BoundingBox, 'EdgeColor', 'g'); end; hold off;
    
    % 3. 展示分支A：面积合格 且 仅经过斜率筛选
    subplot(1, 5, 3);
    keep_branch_A_global = keep_idx_area(setdiff(local_idx_all, drop_idx_slope_local));
    img_branch_A = false(img_H, img_W);
    for i=1:length(keep_branch_A_global), img_branch_A(stats(keep_branch_A_global(i)).PixelIdxList) = true; end
    imshow(~img_branch_A);
    title(sprintf('3. 面积 + 斜率达标\n(剩 %d 个)', length(keep_branch_A_global)));
    hold on; for i=1:length(keep_branch_A_global), rectangle('Position', stats(keep_branch_A_global(i)).BoundingBox, 'EdgeColor', 'b'); end; hold off;
    
    % 4. 展示分支B：面积合格 且 仅经过对称度筛选
    subplot(1, 5, 4);
    keep_branch_B_global = keep_idx_area(setdiff(local_idx_all, drop_idx_sym_local));
    img_branch_B = false(img_H, img_W);
    for i=1:length(keep_branch_B_global), img_branch_B(stats(keep_branch_B_global(i)).PixelIdxList) = true; end
    imshow(~img_branch_B);
    title(sprintf('4. 面积 + 对称度达标\n(剩 %d 个)', length(keep_branch_B_global)));
    hold on; for i=1:length(keep_branch_B_global), rectangle('Position', stats(keep_branch_B_global(i)).BoundingBox, 'EdgeColor', 'm'); end; hold off;
    
    % 5. 最终结果：串联面积 + (斜率并联对称度)
    subplot(1, 5, 5);
    imshow(~B_scan_final);
    title(sprintf('5. 最终混合滤波结果\n(最终保留 %d 个)', length(final_keep_idx)));
    hold on; for i=1:length(final_keep_idx), rectangle('Position', stats(final_keep_idx(i)).BoundingBox, 'EdgeColor', 'r', 'LineWidth', 1.5); end; hold off;
end

end