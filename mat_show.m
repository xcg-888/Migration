% load('UnetIn.mat');  
load('Concrete_Mid_Water_Rebar_0_z_Ez_gaploss33.mat');  

% 2. 图像化显示（核心代码，适配GPR数据特点）
figure('Color','white');  % 创建白色背景图窗（避免默认灰色背景）
imagesc(Ez, [-1, 1]);     % 显示二维数据，强制显示范围为[-1,1]（对应Python归一化）
colormap('gray');            % 用jet伪彩色图（和Python的cmap='jet'一致，清晰区分数值差异）
colorbar;                 % 添加色条（显示颜色对应的数据值，关键！）
title('U-Net Input (GPR Ez Data)');  % 图像标题
xlabel('Trace Index (道数)');         % X轴标签（GPR探测道数）
ylabel('Time Index (时间/深度)');     % Y轴标签（探测时间，对应地下深度）


load('UnetIn.mat');  

% 2. 图像化显示（核心代码，适配GPR数据特点）
figure('Color','white');  % 创建白色背景图窗（避免默认灰色背景）
imagesc(Ez, [-1, 1]);     % 显示二维数据，强制显示范围为[-1,1]（对应Python归一化）
colormap('gray');            % 用jet伪彩色图（和Python的cmap='jet'一致，清晰区分数值差异）
colorbar;                 % 添加色条（显示颜色对应的数据值，关键！）
title('U-Net Output (GPR Ez Data)');  % 图像标题
xlabel('Trace Index (道数)');         % X轴标签（GPR探测道数）
ylabel('Time Index (时间/深度)');     % Y轴标签（探测时间，对应地下深度）