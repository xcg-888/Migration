function [mig_data, z_axis] = stolt_migration(data, dt, dx, eps_r)
    % STOLT_MIGRATION   优化版 Stolt 偏移 (含边界补零与高阶插值)
    % 
    % 参数:
    %   data  : 二维雷达数据矩阵 (行: 时间 t, 列: 空间 x)
    %   dt    : 时间采样率 (单位必须是: 秒 s)
    %   dx    : 空间道间距 (单位: 米 m)
    %   eps_r : 地下介质的相对介电常数 (无量纲)
    %
    % 返回:
    %   mig_data : 偏移后的空间域图像矩阵 (尺寸与原 data 保持一致)
    %   z_axis   : 对应的深度物理坐标轴 (单位: 米 m)

    % 1. 物理参数与速度转换
    c = 3e8;                        % 真空光速 (m/s)
    v_medium = c / sqrt(eps_r);     % 电磁波在介质中的真实传播速度 (m/s)
    
    % GPR 通常记录双程走时，因此偏移用的等效速度是真实波速的一半
    v_mig = v_medium / 2;           
    
    % 获取原始数据尺寸
    [Nt_orig, Nx_orig] = size(data);
    
    % =========================================================================
    % 补丁 A: 边界补零 (Zero-padding) 以消除 2D FFT 的空间与时间周期性环绕假象
    % =========================================================================
    pad_t = Nt_orig;               % 时间维度向下补 1 倍的零 (模拟深层无信号)
    pad_x = round(Nx_orig / 2);    % 空间维度左右各补 0.5 倍的零 (扩展空间孔径)
    
    Nt = Nt_orig + pad_t;          % 补零后的新时间维度
    Nx = Nx_orig + 2 * pad_x;      % 补零后的新空间维度
    
    % 创建扩充矩阵并将原始数据放在正上方中央
    data_padded = zeros(Nt, Nx);
    data_padded(1:Nt_orig, pad_x+1 : pad_x+Nx_orig) = data;
    
    % 2. 定义频率 (omega) 和水平波数 (kx) 轴 (注意：基于补零后的 Nt, Nx)
    df = 1 / (Nt * dt);
    f = (-(Nt/2):(Nt/2-1)) * df; 
    omega = 2 * pi * f(:);          % 列向量 (rad/s)
    
    dkx = 1 / (Nx * dx);
    kx = (-(Nx/2):(Nx/2-1)) * dkx;
    kx = 2 * pi * kx;               % 行向量 (rad/m)
    
    % 3. 2D FFT：将扩充后的数据转换到频率-波数域
    E_wkx = fftshift(fft2(data_padded));
    
    % 4. 定义垂直波数 (kz) 轴
    kz = omega / v_mig; 
    
    % 计算物理深度轴 (仅计算原始矩阵大小所对应的深度)
    dz = v_mig * dt; 
    z_axis = (0:Nt_orig-1) * dz;
    
    % 5. 预分配映射后的矩阵
    E_kzkx = zeros(Nt, Nx);
    
    % 6. 核心步骤：逐列进行 Stolt 映射和插值
    for i = 1:Nx
        kx_val = kx(i);
        
        % 频散关系：计算目标 kz 网格在原始 omega 上的对应点
        w_mapped = v_mig * sign(kz) .* sqrt(kx_val^2 + kz.^2);
        
        % 雅可比变换缩放因子
        jacobian_scale = (v_mig^2 * kz) ./ (w_mapped + eps);
        
        % =====================================================================
        % 补丁 B: 高阶插值 (Spline) 替代线性插值 (Linear)，抑制插值带来的数值振荡
        % =====================================================================
        mapped_column = interp1(omega, E_wkx(:, i), w_mapped, 'spline', 0);
        
        % 应用雅可比能量缩放
        E_kzkx(:, i) = mapped_column .* jacobian_scale;
    end
    
    % 7. 2D 逆 FFT：返回空间深度域
    mig_data_full = real(ifft2(ifftshift(E_kzkx)));
    
    % 8. 裁剪回原始尺寸
    % 砍掉之前我们在两边和下方补充的空缺数据，只保留原始雷达扫描的核心区域
    mig_data = mig_data_full(1:Nt_orig, pad_x+1 : pad_x+Nx_orig);
    
end