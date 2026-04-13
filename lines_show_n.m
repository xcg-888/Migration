function lines_show_n(B_scan_image_Mean_cancel,B_scan_image_edge,lines,max_N)
% 非膨胀
%%%%%%%%%%%%%%%%%%%%%%%%%% 在均值对消后的图像上显示结果 %%%%%%%%%%%%%%%%%%%%%%%%%%
figure; imagesc(B_scan_image_Mean_cancel);colormap('gray'); hold on;        % 显示原图并保持
max_len = norm(lines(1).point1 - lines(1).point2);
xy_long = [lines(1).point1; lines(1).point2];
disp('非膨胀');
% for k = 1:length(lines)
for k = 1:max_N
    % 绘制每条直线段（红色，线宽2）
    % plot([lines(k).point1(1), lines(k).point2(1)], [lines(k).point1(2), lines(k).point2(2)], 'LineWidth', 2, 'Color', 'r');
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

    % Plot beginnings and ends of lines
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    disp(['线段' num2str(k) '-x1:' num2str(xy(1,1))]);
    disp(['线段' num2str(k) '-x2:' num2str(xy(2,1))]);
    % Determine the endpoints of the longest line segment
    len = norm(lines(k).point1 - lines(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','cyan');
title('霍夫变换检测到的直线（非膨胀）');

%%%%%%%%%%%%%%%%%%%%%%%%%% 在边缘检测后的图像上显示结果 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; imagesc(B_scan_image_edge);colormap('gray'); hold on;        % 显示原图并保持
max_len = norm(lines(1).point1 - lines(1).point2);
xy_long = [lines(1).point1; lines(1).point2];
for k = 1:max_N
    % 绘制每条直线段（红色，线宽2）
    % plot([lines(k).point1(1), lines(k).point2(1)], [lines(k).point1(2), lines(k).point2(2)], 'LineWidth', 2, 'Color', 'r');
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

    % Plot beginnings and ends of lines
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

    % Determine the endpoints of the longest line segment
    len = norm(lines(k).point1 - lines(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','cyan');
title('边缘检测-霍夫变换检测到的直线（非膨胀）');

end