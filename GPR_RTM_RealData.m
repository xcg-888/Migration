function [Image_Final, Image_Raw] = GPR_RTM_RealData(Bscan, dt, dx, dz, v_model, f0, offset)
% GPR_RTM_REALDATA 针对真实探地雷达 B-scan 数据的二维逆时偏移成像
%
% 输入参数:
%   Bscan   : 预处理后的雷达数据矩阵 (大小: nt 时间采样点 x nx 空间道数)
%             注意：必须提前去直达波、零点校正、并应用适度的增益(AGC)。
%   dt      : 雷达数据的时间采样间隔 (单位: 秒, 例如 0.1e-9)
%   dx      : 空间道间距 (移动步长) (单位: 米)
%   dz      : 深度网格间距 (单位: 米)
%   v_model : 速度模型矩阵 (大小: nz x nx)。(单位: m/s, 例如真空 3e8)
%             注意：v_model 的列数必须与 Bscan 的道数(nx)一致。
%   f0      : 雷达天线中心频率 (单位: Hz, 例如 400e6)，用于重构震源波场。
%   offset  : 收发天线偏移距 (单位: 米)。若为单天线零偏移距，设为 0。
%
% 输出参数:
%   Image_Final : 经过照明补偿的最终偏移剖面
%   Image_Raw   : 未补偿的原始叠加剖面

    [nt, nx] = size(Bscan);
    [nz, nx_model] = size(v_model);
    
    if nx ~= nx_model
        error('错误: v_model 的列数必须与 Bscan 的道数一致。');
    end

    % 1. CFL 稳定性条件检查
    v_max = max(v_model(:));
    dt_max = 1 / (v_max * sqrt(1/dx^2 + 1/dz^2));
    if dt > dt_max * 0.99
        warning('输入的时间步长 dt (%.3e) 接近或超过 CFL 稳定极限 (%.3e)。如果发生数值频散或发散，请对 Bscan 进行升采样或减小速度模型上限。', dt, dt_max);
    end

    % 2. 吸收边界设置 (海绵层 Sponge Boundary)
    npad = 20; % 吸收层厚度(网格数)
    damp = ones(nz, nx);
    for i = 1:npad
        decay = exp(-(0.015 * (npad - i + 1))^2);
        damp(1:i, :) = damp(1:i, :) .* decay;         % 顶
        damp(end-i+1:end, :) = damp(end-i+1:end, :) .* decay; % 底
        damp(:, 1:i) = damp(:, 1:i) .* decay;         % 左
        damp(:, end-i+1:end) = damp(:, end-i+1:end) .* decay; % 右
    end

    % 3. 重构震源子波 (Ricker)
    t_vec = (0:nt-1) * dt;
    t0 = 1.2 / f0; 
    src = (1 - 2 * (pi * f0 * (t_vec - t0)).^2) .* exp(-(pi * f0 * (t_vec - t0)).^2);

    % 初始化全局成像矩阵
    Image_Raw = zeros(nz, nx);
    Illum_Global = zeros(nz, nx);
    
    disp('开始逐道(Shot-Profile)逆时偏移计算...');
    h = waitbar(0, 'RTM Processing...');
    
    % 4. 逐道循环计算 (实际工程中，可将此 for 改为 parfor 以开启并行池加速)
    for trace_idx = 1:nx
        
        % 确定当前道的 Tx 和 Rx 在网格中的 X 索引
        % 假设当前道中心位置为 trace_idx
        offset_grids = round(offset / dx);
        idx_Tx = max(1, min(nx, trace_idx - round(offset_grids/2)));
        idx_Rx = max(1, min(nx, trace_idx + round(offset_grids/2)));
        
        % 假设天线贴地，深度索引为 2 (留出顶层边界)
        depth_idx = 2; 
        
        % --- A. 震源波场正演 ---
        S_field = zeros(nz, nx, nt); % 为当前道分配正演波场内存
        P_past = zeros(nz, nx); P_now = zeros(nz, nx);
        
        for n = 1:nt
            d2P_dx2 = (circshift(P_now, [0, -1]) - 2*P_now + circshift(P_now, [0, 1])) / dx^2;
            d2P_dz2 = (circshift(P_now, [-1, 0]) - 2*P_now + circshift(P_now, [1, 0])) / dz^2;
            
            P_next = 2*P_now - P_past + (v_model.^2 * dt^2) .* (d2P_dx2 + d2P_dz2);
            P_next(depth_idx, idx_Tx) = P_next(depth_idx, idx_Tx) + src(n) * (v_model(depth_idx,idx_Tx)^2 * dt^2);
            P_next = P_next .* damp;
            
            S_field(:, :, n) = P_now;
            P_past = P_now; P_now  = P_next;
        end
        
        % --- B. 检波器波场逆时反传与互相关成像 ---
        R_past = zeros(nz, nx); R_now  = zeros(nz, nx);
        Image_Local = zeros(nz, nx);
        Illum_Local = zeros(nz, nx);
        
        trace_data = Bscan(:, trace_idx); % 取出当前 A-scan
        
        for n = nt:-1:1
            d2R_dx2 = (circshift(R_now, [0, -1]) - 2*R_now + circshift(R_now, [0, 1])) / dx^2;
            d2R_dz2 = (circshift(R_now, [-1, 0]) - 2*R_now + circshift(R_now, [1, 0])) / dz^2;
            
            R_next = 2*R_now - R_past + (v_model.^2 * dt^2) .* (d2R_dx2 + d2R_dz2);
            
            % 将真实雷达数据作为二次震源逆序注入 Rx 位置
            R_next(depth_idx, idx_Rx) = R_next(depth_idx, idx_Rx) + trace_data(n) * (v_model(depth_idx,idx_Rx)^2 * dt^2);
            R_next = R_next .* damp;
            
            % 提取当前时刻正演波场
            S_now = S_field(:, :, n);
            
            % 零延迟互相关成像与照明度累加
            Image_Local = Image_Local + S_now .* R_now;
            Illum_Local = Illum_Local + S_now.^2;
            
            R_past = R_now; R_now  = R_next;
        end
        
        % --- C. 叠加到全局剖面 ---
        Image_Raw = Image_Raw + Image_Local;
        Illum_Global = Illum_Global + Illum_Local;
        
        % 更新进度条
        if mod(trace_idx, 10) == 0
            waitbar(trace_idx/nx, h, sprintf('RTM Processing... %d / %d', trace_idx, nx));
        end
    end
    close(h);
    
    % 5. 照明度补偿 (去除强浅层能量掩盖)
    epsilon = max(Illum_Global(:)) * 1e-6; % 稳定项
    Image_Final = Image_Raw ./ (Illum_Global + epsilon);
    
    % 绘制结果
    figure('Name', 'RTM Result', 'Position', [100, 100, 800, 600]);
    subplot(2,1,1);
    imagesc(Bscan); colormap('gray'); title('Input B-scan Data');
    xlabel('Trace Number'); ylabel('Time Sample');
    
    subplot(2,1,2);
    imagesc(Image_Final); colormap('gray'); title('RTM Compensated Image');
    xlabel('Trace Number (X)'); ylabel('Depth Grid (Z)');
    % 自动截断异常极值以增强对比度
    caxis([-std(Image_Final(:))*3, std(Image_Final(:))*3]); 
    
    disp('RTM 成像完成。');
end