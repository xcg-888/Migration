function [skeleton, iter_num] = zhang_suen_thinning(bw_img)
    % ZHANG_SUEN_THINNING 实现 Zhang-Suen 并行骨架细化算法
    % 输入：bw_img - 二值图像（逻辑矩阵或数值矩阵，前景=1，背景=0）
    % 输出：skeleton - 细化后的骨架图像（逻辑矩阵）
    %       iter_num - 迭代次数
    
    % 步骤1：输入预处理（确保是二值逻辑矩阵，前景为1，背景为0）
    if ~islogical(bw_img)
        bw_img = logical(bw_img);  % 转为逻辑矩阵
    end
    skeleton = bw_img;  % 初始化骨架为原始图像
    iter_num = 0;       % 迭代次数计数
    change = true;      % 标记是否有像素被删除（迭代终止条件）
    
    % 定义8-邻域索引（以中心像素p1为(0,0)，周围8个像素的相对坐标）
    % p1(中心) p2(上) p3(右上) p4(右) p5(右下) p6(下) p7(左下) p8(左) p9(左上)
    neighbors = [-1,0; -1,1; 0,1; 1,1; 1,0; 1,-1; 0,-1; -1,-1];  % 8个邻域的(y,x)偏移
    
    % 步骤2：迭代细化（直到无像素可删除）
    while change
        change = false;  % 重置“是否变化”标志
        
        % -------------------------- 第1轮迭代：标记要删除的像素 --------------------------
        delete_mask1 = false(size(skeleton));  % 第1轮删除掩码（标记要删除的像素）
        % 遍历图像内部像素（避免边界溢出，不处理边缘1行1列）
        for i = 2:size(skeleton,1)-1
            for j = 2:size(skeleton,2)-1
                if skeleton(i,j)  % 仅处理前景像素
                    % 步骤a：计算8-邻域像素值（p2~p9）
                    p = zeros(1,8);  % 存储p2~p9的像素值（1=前景，0=背景）
                    for k = 1:8
                        y = i + neighbors(k,1);
                        x = j + neighbors(k,2);
                        p(k) = skeleton(y,x);
                    end
                    
                    % 步骤b：计算条件1-4
                    N = sum(p);  % 条件2：8-邻域前景像素数（2≤N≤6）
                    S = 0;       % 条件3：连通数S（0→1的跳变数）
                    for k = 1:7
                        S = S + (p(k) == 0 && p(k+1) == 1);
                    end
                    S = S + (p(8) == 0 && p(1) == 1);  % 最后一个跳变（p8→p2）
                    
                    % 条件4（第1轮）：p2∨p4∨p6∨p8 = 0（偶数位置：p2(1)、p4(3)、p6(5)、p8(7)）
                    cond4 = ~(p(1) || p(3) || p(5) || p(7));
                    
                    % 满足所有条件则标记删除
                    if (N >= 2 && N <= 6) && (S == 1) && cond4
                        delete_mask1(i,j) = true;
                        change = true;
                    end
                end
            end
        end
        skeleton(delete_mask1) = false;  % 第1轮：批量删除标记的像素
        
        % -------------------------- 第2轮迭代：标记要删除的像素 --------------------------
        delete_mask2 = false(size(skeleton));  % 第2轮删除掩码
        for i = 2:size(skeleton,1)-1
            for j = 2:size(skeleton,2)-1
                if skeleton(i,j)  % 仅处理前景像素
                    % 计算8-邻域像素值（p2~p9）
                    p = zeros(1,8);
                    for k = 1:8
                        y = i + neighbors(k,1);
                        x = j + neighbors(k,2);
                        p(k) = skeleton(y,x);
                    end
                    
                    % 计算条件1-4
                    N = sum(p);
                    S = 0;
                    for k = 1:7
                        S = S + (p(k) == 0 && p(k+1) == 1);
                    end
                    S = S + (p(8) == 0 && p(1) == 1);
                    
                    % 条件4（第2轮）：p2∨p3∨p6∨p7 = 0（p2(1)、p3(2)、p6(5)、p7(6)）
                    cond4 = ~(p(1) || p(2) || p(5) || p(6));
                    
                    % 满足所有条件则标记删除
                    if (N >= 2 && N <= 6) && (S == 1) && cond4
                        delete_mask2(i,j) = true;
                        change = true;
                    end
                end
            end
        end
        skeleton(delete_mask2) = false;  % 第2轮：批量删除标记的像素
        
        iter_num = iter_num + 1;  % 迭代次数+1
    end
    
    % 步骤3：输出结果（确保骨架为单像素宽度）
    skeleton = logical(skeleton);
end

% =========================================================================
% 测试代码：运行 Zhang-Suen 算法并显示结果
% =========================================================================
clc; clear; close all;

% 1. 读取并预处理图像（替换为你的图像路径）
img_path = 'test_shape.png';  % 测试图像（建议用二值图，如手写数字、几何形状）
img = imread(img_path);
if size(img,3) == 3  % 彩色图像转灰度图
    img_gray = rgb2gray(img);
else
    img_gray = img;
end
bw_img = im2bw(img_gray, 0.5);  % 二值化（阈值0.5可调整，前景=1，背景=0）
bw_img = imcomplement(bw_img);  % 若目标是黑色，反转图像（确保目标为1）

% 2. 运行 Zhang-Suen 细化算法
[skeleton, iter_num] = zhang_suen_thinning(bw_img);

% 3. 显示结果
figure('Color','white','Position',[100,100,800,300]);
subplot(1,3,1); imshow(img_gray); title('原始灰度图');
subplot(1,3,2); imshow(bw_img); title('二值图像');
subplot(1,3,3); imshow(skeleton); title(['Zhang-Suen 骨架（迭代',num2str(iter_num),'次）']);

% 4. 可选：保存骨架结果
% imwrite(uint8(skeleton)*255, 'skeleton_result.png');