function image_show2(col,row,B_scan_image_down,B_scan_image_down_Mean_cancel,target_area,q_n,dt,Downsample_N,R,dx,q)
%% 第三代双曲线
    comp = 0;
    L = 150;% 双曲线展示窗口大小的一半
    epsilon = 2.5:0.1:14;
    
    figure;
    imagesc(B_scan_image_down);colormap('gray');hold on;   
    scatter(col,row,20,'r','filled');   
    title(['第三代-下采样-figure:',num2str(q_n),'  epsilon:',num2str(epsilon(q_n))]);
    % 绘制双曲线
    for i = 1:length(col)
        min_x = max(col(i)-L,1);
        max_x = min(col(i)+L,size(B_scan_image_down,2));
        for j = min_x:max_x  
            distance_io = 0.04 ; % 发射源与接收源的间距
            v = 1/sqrt(q(q_n));

            t1 = 1/v*(  sqrt( ((row(i)*dt*Downsample_N-comp)/2*v+R)^2 + (j*dx-col(i)*dx)^2 )  -R  );   % 前半程时间

            sin_a = (row(i)*dt*Downsample_N*v/2+R)/(t1*v+R);
            cos_a = (j*dx-col(i)*dx)/(t1*v+R);

            t2 = 1/v*sqrt(  (t1*v*sin_a)^2  +  (j*dx-col(i)*dx+distance_io-R*cos_a)^2  );  %后半程时间
            y = round(   (t1+t2)/dt/Downsample_N+comp/dt   );
            
            if (y > 0) && (y <= size(B_scan_image_down,1)) 
                scatter(j,y,2,'y','filled');
            end
        end
    end


    figure;
    imagesc(B_scan_image_down_Mean_cancel);colormap('gray');hold on;   
    scatter(col,row,20,'r','filled');   
    title(['第三代-均值对消-figure:' num2str(q_n),'  epsilon:',num2str(epsilon(q_n))]);
    % 绘制双曲线
    for i = 1:length(col)
        min_x = max(col(i)-L,1);
        max_x = min(col(i)+L,size(B_scan_image_down_Mean_cancel,2));
        for j = min_x:max_x
            distance_io = 0.04 ; % 发射源与接收源的间距
            v = 1/sqrt(q(q_n));

            t1 = 1/v*(  sqrt( ((row(i)*dt*Downsample_N-comp)/2*v+R)^2 + (j*dx-col(i)*dx)^2 )  -R  );   % 前半程时间

            sin_a = (row(i)*dt*Downsample_N*v/2+R)/(t1*v+R);
            cos_a = (j*dx-col(i)*dx)/(t1*v+R);

            t2 = 1/v*sqrt(  (t1*v*sin_a)^2  +  (j*dx-col(i)*dx+distance_io-R*cos_a)^2  );  %后半程时间
            y = round(   (t1+t2)/dt/Downsample_N+comp/dt   );
            
            if (y > 0) && (y <= size(B_scan_image_down_Mean_cancel,1)) 
                scatter(j,y,2,'y','filled');
            end
        end
    end


    figure;
    imagesc(target_area);colormap('gray');hold on;   
    scatter(col,row,20,'r','filled');   
    title(['第三代-霍夫空间-figure:' num2str(q_n),'  epsilon:',num2str(epsilon(q_n))]);
    % 绘制双曲线
    for i = 1:length(col)
        min_x = max(col(i)-L,1);
        max_x = min(col(i)+L,size(target_area,2));
        for j = min_x:max_x
            distance_io = 0.04 ; % 发射源与接收源的间距
            v = 1/sqrt(q(q_n));

            t1 = 1/v*(  sqrt( ((row(i)*dt*Downsample_N-comp)/2*v+R)^2 + (j*dx-col(i)*dx)^2 )  -R  );   % 前半程时间

            sin_a = (row(i)*dt*Downsample_N*v/2+R)/(t1*v+R);
            cos_a = (j*dx-col(i)*dx)/(t1*v+R);

            t2 = 1/v*sqrt(  (t1*v*sin_a)^2  +  (j*dx-col(i)*dx+distance_io-R*cos_a)^2  );  %后半程时间
            y = round(   (t1+t2)/dt/Downsample_N+comp/dt   );
                       
            if (y > 0) && (y <= size(target_area,1)) 
                scatter(j,y,2,'y','filled');
            end
        end
    end

end