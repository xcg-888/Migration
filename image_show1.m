function image_show1(col,row,B_scan_image_down,B_scan_image_down_Mean_cancel,target_area,q_n,dt,Downsample_N,R,dx,q)
%% 第二代双曲线
    comp = 0;
    L = 150;% 双曲线展示窗口大小的一半
    epsilon = 2.5:0.1:14;

    figure;
    imagesc(B_scan_image_down);colormap('gray');hold on;   
    scatter(col,row,20,'r','filled');   
    title(['第二代-下采样-figure:',num2str(q_n),'  epsilon:',num2str(epsilon(q_n))]);
    % 绘制双曲线
    for i = 1:length(col)
        min_x = max(col(i)-L,1);
        max_x = min(col(i)+L,size(B_scan_image_down,2));
        for j = min_x:max_x
            y = round(  (sqrt(((row(i)*dt*Downsample_N-comp)/2/sqrt(q(q_n))+R)^2+(j*dx-col(i)*dx)^2)-R) *2*sqrt(q(q_n))/dt/Downsample_N+comp/dt);
            if (y > 0) && (y <= size(B_scan_image_down,1)) 
                scatter(j,y,2,'y','filled');
            end
        end
    end


    figure;
    imagesc(B_scan_image_down_Mean_cancel);colormap('gray');hold on;   
    scatter(col,row,20,'r','filled');   
    title(['第二代-均值对消-figure:' num2str(q_n),'  epsilon:',num2str(epsilon(q_n))]);
    % 绘制双曲线
    for i = 1:length(col)
        min_x = max(col(i)-L,1);
        max_x = min(col(i)+L,size(B_scan_image_down_Mean_cancel,2));
        for j = min_x:max_x
            y = round(  (sqrt(((row(i)*dt*Downsample_N-comp)/2/sqrt(q(q_n))+R)^2+(j*dx-col(i)*dx)^2)-R) *2*sqrt(q(q_n))/dt/Downsample_N+comp/dt);
            if (y > 0) && (y <= size(B_scan_image_down_Mean_cancel,1)) 
                scatter(j,y,2,'y','filled');
            end
        end
    end


    figure;
    imagesc(target_area);colormap('gray');hold on;   
    scatter(col,row,20,'r','filled');   
    title(['第二代-霍夫空间-figure:' num2str(q_n),'  epsilon:',num2str(epsilon(q_n))]);
    % 绘制双曲线
    for i = 1:length(col)
        min_x = max(col(i)-L,1);
        max_x = min(col(i)+L,size(target_area,2));
        for j = min_x:max_x
            y = round(  (sqrt(((row(i)*dt*Downsample_N-comp)/2/sqrt(q(q_n))+R)^2+(j*dx-col(i)*dx)^2)-R) *2*sqrt(q(q_n))/dt/Downsample_N+comp/dt);
            if (y > 0) && (y <= size(target_area,1)) 
                scatter(j,y,2,'y','filled');
            end
        end
    end

end

