function migration_space = Hyperbolic_Diffraction_Summation_total(Ydata,dt,dx,q,comp,Downsample_N,R)
%% 函数核心功能：双曲线绕射叠加（全图网格遍历版）
%% 输入
% Ydata：B-scan预处理后的数据
% dt：快时间采样时间间隔，单位ns
% dx：道间距，单位m
% q：速度参数的扫描范围，当介电常数已知时，可以直接设置q = epsilon*(1 / 0.3)^2;
% comp：时间补偿因子，单位ns, 如过不做补偿可设置为0
% Downsample_N：下采样阶数
% R：目标半径
%% 输出
% migration_space：参数空间，其中的极大值点反映了双曲线对应的[时间，雷达道]
%% 读取参数
    s = Ydata;
    m = size(s,1);
    n = size(s,2);
    migration_space = zeros(m,n); % 参数统计后的矩阵
    
%% 参数空间投票（全图遍历）
    % 逻辑：遍历图像上的每一个网格，假设它都是潜在的顶点，然后按理论双曲线轨迹拾取原图能量
    for i = 1:m      % 遍历所有行（快时间轴/深度）
        for j = 1:n  % 遍历所有列（雷达道/空间位置）
            
            % 设定当前假设顶点(i,j)的偏移孔径范围，防止溢出
            low = max(j-75, 1);    
            high = min(j+75, n);   
            
            % 使用局部变量累加当前顶点的能量，这比直接在矩阵上累加(migration_space(i,j)=...)速度快很多
            temp_energy_sum = 0; 
            
            for k = low:high  % 遍历孔径内的每一道
                
                % 第二代公式（考虑管线半径，不考虑发射源与接受点的间距）：
                % i: 当前假设顶点的行坐标
                % j: 当前假设顶点的列坐标
                % k: 当前正在积分（拾取能量）的列坐标
                row = round(  (sqrt(((i*dt*Downsample_N-comp)/2/sqrt(q)+R)^2+(k*dx-j*dx)^2)-R) *2*sqrt(q)/dt/Downsample_N+comp/dt);
                
                % 边界保护：如果计算出的双曲线轨迹点在图像范围内
                if (row > 0) && (row <= m) 
                    % 拾取原图中该轨迹点上的绝对能量并累加
                    temp_energy_sum = temp_energy_sum + abs(s(row,k));
                end
            end
            
            % 将累加好的能量赋给当前假设的顶点
            migration_space(i,j) = temp_energy_sum;
            
        end
    end
end