function image_out = Kirchhoff_Migration(data_in, TrackInterval, dt, epsilon_r, Radius)
% 柯希霍夫偏移 (Kirchhoff Migration) MATLAB 实现 - 【含目标半径修正版】
% 作用：将具有物理半径的管状目标的回波能量，精确聚焦到其几何中心。
%
% 输入:
%   data_in       : 原始 B-scan 雷达数据矩阵 (Size: Nt x Nx)
%   TrackInterval : 走航道间距 dx (单位: 米)
%   dt            : 时间轴采样间隔 (单位: 纳秒 ns)  <--- 强烈建议用 ns
%   epsilon_r     : 相对介电常数
%   Radius        : 目标物理半径 (单位: 米)
%
% 输出:
%   image_out     : 偏移成像后的矩阵 (Size: Nt x Nx)

    % --- 1. 设定物理参数 ---
    c = 3e8; % 光速 (单位: 米/纳秒 m/ns，防止极小浮点数截断误差)
    v = c / sqrt(epsilon_r); % 介质中的电磁波速
    v_m = v / 2; % 偏移速度 (单程波速)
    
    [Nt, Nx] = size(data_in);
    x_vec = (0 : Nx-1) * TrackInterval; 
    image_out = zeros(Nt, Nx);
    
    % --- 2. 预处理：沿时间轴做一阶差分 ---
    D_prime = diff(data_in, 1, 1); 
    D_prime = [D_prime; zeros(1, Nx)]; 

    % --- 3. 主循环遍历 ---
    for ix = 1:Nx
        x0 = x_vec(ix);
        for iz = 1:Nt
            % 将行索引 iz 转换为目标【几何中心】的物理深度 z0 (米)
            z0 = (iz - 1) * dt * v_m; 
            
            % 如果假设的中心深度还不如半径大，说明管子露在地面上了，没有物理意义
            if z0 <= Radius
                continue; 
            end
            
            % --- 4. 向量化的孔径循环 ---
            dx_vec = x_vec - x0; 
            
            % 【修正处 1：几何关系】
            % a. 先计算天线到管线几何中心的绝对距离
            R_center = sqrt(dx_vec.^2 + z0^2);
            % b. 减去管线半径，得到实际发生反射的最短路径距离
            R_surface = R_center - Radius;
            
            % 安全保护：如果计算出的表面距离小于0（天线插进管子里了），则剔除
            R_surface(R_surface < 0) = NaN; 
            
            % 计算单向走时 t 数组 (单位: ns)
            t_vec = R_surface / v_m;
            
            % 走时转索引 (带小数)
            j_exact = (t_vec / dt) + 1; 
            
            % 剔除超出数据时间窗口或含有 NaN 的无效索引
            valid_mask = (j_exact >= 1) & (j_exact <= Nt) & ~isnan(j_exact);
            j_valid = j_exact(valid_mask);
            i_valid = find(valid_mask); 

            % 以上步骤实现了双曲线轨迹寻址
            
            if isempty(i_valid)
                continue;
            end
            
            % --- 5. 提取数据与一维线性插值 ---
            j_floor = floor(j_valid);           % 向下取整
            j_ceil = min(j_floor + 1, Nt);      % 向上取整
            weight_ceil = j_valid - j_floor;    % 插值比例（上）
            weight_floor = 1 - weight_ceil;     % 插值比例（下）
            
            idx_floor = sub2ind([Nt, Nx], j_floor, i_valid);
            idx_ceil  = sub2ind([Nt, Nx], j_ceil, i_valid);
            
            A = D_prime(idx_floor) .* weight_floor + D_prime(idx_ceil) .* weight_ceil;
            
            % --- 6. 计算权重并加权累加 ---
            % 权重分母应当依然使用到中心的距离 R_center 来计算球面发散，
            % 因为波束的角度扩散是以中心为参考系的。
            W = z0 ./ (R_center(valid_mask).^2); 
            
            sum_energy = sum(W .* A .* TrackInterval); 
            image_out(iz, ix) = sum_energy;
        end
    end
end