function [mig_data_theta,mig_data] = ratate_migration(data, dt, dx, eps_r)
    % RATATE_MIGRATION   旋转偏移
    % 
    % 优化说明:
    % 1. 采用全孔径遍历，无空间截断，保证算法的完全普适性。
    % 2. 加入一维线性插值，消除阶梯状伪影。
    % 3. 保留克希霍夫积分理论中的“倾角因子”(Obliquity Factor)。
    %
    % 参数:
    %   data  : 二维雷达数据矩阵 (行: 时间 t, 列: 空间 x)
    %   dt    : 时间采样率 (单位必须是: a秒 s)
    %   dx    : 空间道间距 (单位: 米 ma)
    %   eps_r : 地下介质的相对介电常数 (无量纲)
    %
    % 返回:
    %   mig_data : 偏移后的空间域图像矩阵 (尺寸与原 data 保持一致)
    
    % --- 1. 物理参数与速度转换 ---
    c = 3e8;                        % 真空光速 (m/s)
    v_medium = c / sqrt(eps_r);     % 电磁波在介质中的真实传播速度 (m/s)
    v_mig = v_medium / 2;           % GPR双程走时等效偏移速度
    
    % 获取原始数据尺寸
    [Nt, Nx] = size(data);
  
    % 初始化偏移数据矩阵
    mig_data_theta = zeros(Nt, Nx);
    mig_data = zeros(Nt, Nx);
    
    % 预生成所有探地雷达的空间序列
    i_vec = 1:Nx; 
    
    % --- 2. 核心偏移累加过程 (全孔径) ---
    for xi = 1:Nx
        for ti = 1:Nt
            
            % 向量化计算当前像素点 (ti, xi) 到所有雷达天线位置的浮点时间索引
            t_float = sqrt( ((xi - i_vec) .* dx ./ (v_mig * dt)).^2 + ti^2 );
            
            % 全域积分
            valid_mask = (t_float >= 1) & (t_float <= Nt);
            
            % 提取有效数据索引
            valid_i = i_vec(valid_mask);
            valid_t = t_float(valid_mask);
            
            % 如果当前点周围没有有效数据，直接跳过
            if isempty(valid_i)
                continue;
            end
            
            % --- “对半分”处理：一维线性插值提取 ---
            t_floor = floor(valid_t);
            t_ceil  = min(Nt, t_floor + 1); 
            
            weight_ceil  = valid_t - t_floor;
            weight_floor = 1 - weight_ceil;
            
            % 二维坐标转一维线性索引加速提取 (idx = row + (col - 1)*total_rows)
            idx_floor = t_floor + (valid_i - 1) * Nt;
            idx_ceil  = t_ceil  + (valid_i - 1) * Nt;
            
            % 插值提取原始数据的振幅值
            extracted_vals = data(idx_floor) .* weight_floor + data(idx_ceil) .* weight_ceil;
            
            % 倾角因子 (Obliquity Factor): cos(theta) = ti / valid_t
            % 即使在全孔径下，这也是克希霍夫积分必需的物理权重
            obliquity_weight = ti ./ valid_t; 

            % 乘上倾角权重后，全孔径累加到当前像素点
            mig_data_theta(ti, xi) = sum(extracted_vals .* obliquity_weight);
            mig_data(ti, xi) = sum(extracted_vals);
        end
    end
end