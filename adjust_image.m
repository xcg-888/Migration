function img_final = adjust_image(gpr_data,min,max)
% min = 0.02;
% max = 0.7;
% 1. 取绝对值（只看能量大小，忽略正负）
img_amp = abs(gpr_data);

% 2. 归一化到 0-1 之间
img_norm = mat2gray(img_amp);

% 3. 【核心步骤】非线性增强与去噪
% imadjust(图片, [低阈值 高阈值], [输出范围])
% 下面代码的意思是：把 10% 以下的弱信号直接变成 0（去灰雾），把 50% 以上的信号直接变成 1（增强目标）
img_adjust = imadjust(img_norm, [min max], []); 

% 4. 颜色反转（黑底变白底，白目标变黑目标）
img_final = imcomplement(img_adjust); 

end