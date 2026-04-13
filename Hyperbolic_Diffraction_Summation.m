function migration_space = Hyperbolic_Diffraction_Summation(Ydata,dt,dx,q,comp,Downsample_N,R)
%% 函数核心功能：使用Hough变换检测双曲线，也可以理解为模板匹配
%% 输入
% Ydata：B-scan预处理后的二值化数据
% dt：快时间采样时间间隔，单位ns
% dx：道间距，单位m
% q：速度参数的扫描范围，当介电常数已知时，可以直接设置q = epsilon*(1 / 0.3)^2;
% comp：时间补偿因子，单位ns, 如过不做补偿可设置为0
% Downsample_N：下采样阶数
% R：目标半径
%% 输出
% hough_space：参数空间，其中的极大值点反映了双曲线对应的[时间，雷达道，介电常数]

%% 读取参数
    s = Ydata;
    m = size(s,1);
    n = size(s,2);
    migration_space=zeros(m,n); % 参数统计后的矩阵
    [rows,cols] = find(s); % B-scan每个点的行和列索引
    ecount = length(rows);

%% 参数空间投票
    % 逻辑：先假设扫描到的每一个点都是顶点，然后按照双曲线的形状累加双曲线的值，如果该点真的是顶点，那么会出现一个极值
    for i=1:ecount  % 扫描B-scan中的每个点
        low = max(cols(i)-75,1);    % 防溢出处理
        high = min(cols(i)+75,n);   % 防溢出处理
        for j = low:high
            % 第一代公式（不考虑管线半径，，不考虑发射源与接受点的间距）：
            % row = round(  sqrt((rows(i)*dt*Downsample_N-comp)^2/4+q(k)*(j*dx-cols(i)*dx)^2)  *2/dt/Downsample_N+comp/dt); % 假设该点为双曲线顶点，计算每个A-scan对应的采样点
            % row = round( sqrt((rows(i)*dt-comp)^2 + 4*q(k)*(j*dx-cols(i)*dx)^2)/dt+comp/dt ); % 假设该点为双曲线顶点，计算每个A-scan对应的采样点

            % 第二代公式（考虑管线半径，不考虑发射源与接受点的间距）：
            row = round(  (sqrt(((rows(i)*dt*Downsample_N-comp)/2/sqrt(q)+R)^2+(j*dx-cols(i)*dx)^2)-R) *2*sqrt(q)/dt/Downsample_N+comp/dt);
            
            % 第三代公式（考虑管线半径，考虑发射源与接受点的间距）
            % distance_io = 0.04 ; % 发射源与接收源的间距
            % v = 1/sqrt(q(k));
            % 
            % t1 = 1/v*(  sqrt( ((rows(i)*dt*Downsample_N-comp)/2*v+R)^2 + (j*dx-cols(i)*dx)^2 )  -R  );   % 前半程时间
            % 
            % sin_a = (rows(i)*dt*Downsample_N*v/2+R)/(t1*v+R);
            % cos_a = (j*dx-cols(i)*dx)/(t1*v+R);
            % 
            % t2 = 1/v*sqrt(  (t1*v*sin_a)^2  +  (j*dx-cols(i)*dx+distance_io-R*cos_a)^2  );  %后半程时间
            % row = round(   (t1+t2)/dt/Downsample_N+comp/dt   );

            % 模板匹配
            if (row > 0) && (row <= m) 
                migration_space(rows(i),cols(i))=migration_space(rows(i),cols(i)) + abs(s(row,j));

            % % 模板范围匹配，即双曲线变粗
            % range_N = 2;
            % for range_down = range_N:-1:1
            %     if row - range_down > 0
            %         for range_down_n = 1:range_down
            %             hough_space(rows(i),cols(i),k) = hough_space(rows(i),cols(i),k) + s(row-range_down_n,j);
            %         end
            %         break;
            %     end
            % end
            % 
            % for range_up = range_N:-1:1
            %     if row + range_up < m
            %         for range_up_n = 1:range_up
            %             hough_space(rows(i),cols(i),k) = hough_space(rows(i),cols(i),k) + s(row+range_up_n,j);
            %         end
            %         break;
            %     end
            % end

            end
        end
    end

%% 双曲线检验——以(150,200)为顶点画双曲线
    % for k=1:length(q) % 扫描介质参数空间
    %     for j = 1:n
    %         % row = round(  sqrt((5*dt*Downsample_N-comp)^2/4+q(k)*(j*dx-150*dx)^2)  *2/dt/Downsample_N+comp/dt); % 假设该点为双曲线顶点，计算每个A-scan对应的采样点
    %         % row = round( sqrt((rows(i)*dt-comp)^2 + 4*q(k)*(j*dx-cols(i)*dx)^2)/dt+comp/dt ); % 假设该点为双曲线顶点，计算每个A-scan对应的采样点        
    %         row = round(  (sqrt(((440*dt*Downsample_N-comp)/2/sqrt(q(k))+R)^2+(j*dx-125*dx)^2)-R) *2*sqrt(q(k))/dt/Downsample_N+comp/dt);
    %         if (row > 0) && (row <= m) 
    %             hough_space(row,j,k)=hough_space(row,j,k) + 100;
    %         end
    %     end
    % end
   
end

