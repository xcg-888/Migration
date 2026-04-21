function [migrated_img] = phase_shift_migration(data_in, dx, dt, epsilon_r)
    % PHASE_SHIFT_MIGRATION 相移偏移算法 (Gazdag)
    % 输入:
    %   data_in   - 二维雷达原始数据 (行: 时间 t, 列: 空间位置 x)
    %   dx        - 空间道间距 (单位: m)
    %   dt        - 时间采样间隔 (单位: 秒 s)  <-- 强烈建议统一步伐至国际标准单位
    %   epsilon_r - 相对介电常数
    
    c = 3e8; % 光速 (单位: 米/秒 m/s)
    v = c / sqrt(epsilon_r); % 介质中的电磁波速 (m/s)
    [Nt, Nx] = size(data_in);
    
    % --- 1. 构建频率轴 w 和水平波数轴 kx ---
    df = 1 / (Nt * dt);
    f = (0:Nt-1)' * df;
    f(f >= 1/(2*dt)) = f(f >= 1/(2*dt)) - 1/dt; % 修正奈奎斯特边界的取等号更严谨
    w = 2 * pi * f; 
    
    dk = 1 / (Nx * dx);
    k_x = (0:Nx-1) * dk;
    k_x(k_x >= 1/(2*dx)) = k_x(k_x >= 1/(2*dx)) - 1/dx;
    k_x = 2 * pi * k_x; 
    
    % --- 2. 转换到 f-k 域 ---
    P_wk = fft2(data_in); 
    
    % --- 3. 准备向下延拓的参数 ---
    v_esm = v / 2; % 爆炸反射面等效速度
    dz = v_esm * dt; 
    Nz = Nt;         
    
    [Kx_grid, W_grid] = meshgrid(k_x, w);
    
    % --- 4. 倏逝波(Evanescent Waves)处理 ---
    inner_term = (W_grid ./ v_esm).^2 - Kx_grid.^2;
    evanescent_mask = (inner_term < 0);
    
    kz = zeros(size(inner_term));
    kz(~evanescent_mask) = sqrt(inner_term(~evanescent_mask));
    
    % --- 关键修复：加入 sign(W_grid) 保证共轭对称性 ---
    K = exp(1j * sign(W_grid) .* kz * dz);
    K(evanescent_mask) = 0; 
    
    % --- 5. 逐层向下延拓与成像 ---
    migrated_img = zeros(Nz, Nx);
    P_current = P_wk; 
    
    for iz = 1:Nz
        P_current = P_current .* K;
        P_kx_t0 = sum(P_current, 1); 
        migrated_img(iz, :) = real(ifft(P_kx_t0)); 
    end
    
    migrated_img = migrated_img / max(abs(migrated_img(:)));
end