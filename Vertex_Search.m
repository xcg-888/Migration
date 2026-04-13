%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    顶点搜寻函数                                                                             
% target_area:在hough_space同位置，保留大于阈值的元素，小于阈值的元素清零                   
% target_area_binary：在hough_space同位置，大于阈值的元素置一，小于阈值的元素置零           
% T1_threshold：霍夫空间显示阈值                                                                     
% T2_threshold：顶点搜寻阈值                                                                      
% row：顶点的快时间轴坐标                                                                
% col：顶点的慢时间轴坐标
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 函数
function [target_area,row,col] = Vertex_Search(temp)
    T1_threshold = 0.05; 
    T2_threshold = 0.70;
    target_area = zeros(size(temp)); 
    target_area_binary = zeros(size(temp));
    T1 = max(max(temp))*T1_threshold;     % 阈值等于全局最大值的T1_threshold倍
    T2 = max(max(temp))*T2_threshold;     % 阈值等于全局最大值的T2_threshold倍
    [rows,cols] = find(temp>T1); % 大于阈值为目标，并记录坐标
    len = size(rows,1);         % 大于阈值的点的个数
    for j = 1:len
        target_area(rows(j),cols(j)) = temp(rows(j),cols(j));   
        target_area_binary(rows(j),cols(j)) = 1;
    end
    
    % [row,col]代表了每个聚类极大值点的坐标
    [L,V] = bwlabel(target_area_binary);   % 8连通矩阵检索，聚类 L：聚类矩阵  V：类的数目
    row = zeros(1,V);
    col = zeros(1,V);
    for v = 1:V
        target_area_temp = target_area.*(L==v);
        [tr,tc] = find(target_area_temp == max(max(target_area_temp)));
        if temp(tr(1,1),tc(1,1)) > T2
            row(v) = tr(1,1);
            col(v) = tc(1,1);
        end
    end
    % 找出每个类里面，值最大的点
 
    row = flip(row);
    col = flip(col);
end
