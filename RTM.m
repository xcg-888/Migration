function [R_image, x, z] = RTM(data, t, ep_model, x_max, z_max, dx, dz, dx_trace, n_c, npml)
% GPR_RTM - 封装好的探地雷达逆时偏移(RTM)主函数
%
% 输入参数:
%   data     - GPR原始数据矩阵 (B-scan)
%   t        - 时间向量 (单位: 纳秒 ns)
%   ep_model - 真实相对介电常数模型矩阵 (大小必须为 [nz, nx])
%   x_max    - 模型最大水平距离 (m)
%   z_max    - 模型最大深度 (m)
%   dx       - 水平网格步长 (m)
%   dz       - 垂直网格步长 (m)
%   dx_trace - 雷达道间距 (m)
%   n_c      - 需要切除的直达波/延迟时间采样点数
%   npml     - [可选] PML吸收边界层数 (默认: 40)

    % 检查可选参数 npml
    if nargin < 10
        npml = 40; 
    end
    
    % --- 1. 数据预处理 ---
    data = BGRm(data);      % 背景去除
    
    % 将纳秒转换为秒
    t = t .* 1e-9;          
    
    % 安全地切除直达波延迟 (兼容 t 是行向量或列向量的情况)
    if n_c > 0 && n_c < size(data, 1)
        data(1:n_c, :) = 0; % 一定要保证这里是置零而不是删掉整行！
    end

    % --- 2. 建立空间网格 ---
    x = 0:dx:x_max;
    z = 0:dz:z_max;
    nx = size(x, 2);
    nz = size(z, 2);
    
    % 校验输入的介电常数矩阵大小是否匹配
    if any(size(ep_model) ~= [nz, nx])
        error('输入的 ep_model 尺寸 [%d, %d] 与网格划分尺寸 [%d, %d] 不匹配！', ...
            size(ep_model, 1), size(ep_model, 2), nz, nx);
    end
    
    % --- 3. 设置物理参数矩阵 ---
    % 【致命错误修复】：恢复半速模型 (Exploding Reflector Model)
    % 零偏移距雷达记录的是双程走时。为了用标准波方程做偏移，必须将全域波速减半。
    % 因为 v = c / sqrt(ep)，速度减半等效于将相对介电常数 ep 乘 4！
    % 这同时也大大放宽了 FDTD 的 CFL 稳定性限制。
    ep = (ep_model .* 4)';
    
    mu = ones(nz, nx)';
    sig = 0.001 * ones(nz, nx)';
    
    % 插值到更细的 FDTD 网格
    x2 = min(x):dx/2:max(x);
    z2 = min(z):dx/2:max(z);
    ep2 = gridinterp(ep, x, z, x2, z2, 'nearest');
    mu2 = gridinterp(mu, x, z, x2, z2, 'nearest');
    sig2 = gridinterp(sig, x, z, x2, z2, 'nearest');
    
    % --- 4. 扩展 PML 边界 ---
    [ep3, x3, z3] = padgrid(ep2, x2, z2, 2*npml);
    [mu3, x3, z3] = padgrid(mu2, x2, z2, 2*npml);
    [sig3, x3, z3] = padgrid(sig2, x2, z2, 2*npml);
    
    % --- 5. 设置收发天线位置矩阵 ---
    srcx = (0:dx_trace:x_max)';
    srcz = zeros(size(srcx));
    recx = srcx; % 零偏移距配置
    recz = srcz;
    srcloc = [srcx, srcz];
    recloc = [recx, recz];
    
    % --- 6. 运行核心 RTM 算法 ---
    disp('开始运行有限差分波场逆推 (FDTD RTM)...');
    tic;
    % 调用你的核心底层函数 GPR_RTM
    [R_pro, ~, ~, ~, ~] = GPR_RTM(ep3, mu3, sig3, x3, z3, srcloc, recloc, t, npml, data);
    disp(['FDTD 计算完毕，总运行时间 = ', num2str(toc/3600), ' 小时']);
    
    % --- 7. 提取结果 ---
    [nz_pro, nx_pro] = size(R_pro);
    % 截取掉外围的 PML 边界，只返回核心成像区域，转置回常规 [深度 x 距离] 矩阵
    R_image = R_pro(npml:nz_pro-npml, npml:nx_pro-npml);
end