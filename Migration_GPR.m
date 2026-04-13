clear;
close all;
clc;

% python -m gprMax /WorkSpaceP1/XuCaiGuo/updata/model.in -n 336 -gpu
% python -m tools.outputfiles_merge /WorkSpaceP1/XuCaiGuo/updata/model

% python -m gprMax /workspace/xcg/updata/model.in -n 336 -gpu
% python -m tools.outputfiles_merge /workspace/xcg/updata/model

%% 读取文件
DZT_path = 'D:\gprMax\数据库\2024_11_1\ZSX__019.DZT';
lte_path = 'D:\硕士文件\gprMAX相关\实测数据\lte良好数据整合\ltefile596.lte';
mat_path = 'C:\Users\86187\Desktop\single_crack_1.6GHz_slice1.mat';
out_path = 'D:\gprMax\model1\a\model1a2DBscan8\model1a_merged.out';
png_path = 'D:\gprMax\model1\a\model1a2DBscan8\model1a_merged.out'; 
image_path = 'C:\My_workfile\毕业论文\第一个创新点钢筋网仿真in文件\新结果图\实测实验\霍夫变换\';

%% 参数设置
image_save = "off";  % 是否保存图片
mode = 'out';        % 选择输入什么后缀的文件

%% 霍夫变换（内含图像预处理）
if strcmp(mode,'lte')
    [hough_space1,hough_space2,B_scan_image_down_Mean_cancel,B_scan_image_down,Downsample_N,Relative_permittivity1,Relative_permittivity2,R,dx,dt,alpha1,alpha2,n,epsilon_r] = generate_hough_space(lte_path,mode,image_path,image_save);
elseif strcmp(mode,'DZT')
    [hough_space1,hough_space2,B_scan_image_down_Mean_cancel,B_scan_image_down,Downsample_N,Relative_permittivity1,Relative_permittivity2,R,dx,dt,alpha1,alpha2,n,epsilon_r] = generate_hough_space(DZT_path,mode,image_path,image_save);
elseif strcmp(mode,'mat')
    [hough_space1,hough_space2,B_scan_image_down_Mean_cancel,B_scan_image_down,Downsample_N,Relative_permittivity1,Relative_permittivity2,R,dx,dt,alpha1,alpha2,n,epsilon_r] = generate_hough_space(mat_path,mode,image_path,image_save);
elseif strcmp(mode,'out')
    [hough_space1,hough_space2,B_scan_image_down_Mean_cancel,B_scan_image_down,Downsample_N,Relative_permittivity1,Relative_permittivity2,R,dx,dt,alpha1,alpha2,n,epsilon_r] = generate_hough_space(out_path,mode,image_path,image_save);
elseif strcmp(mode,'png')
    [hough_space1,hough_space2,B_scan_image_down_Mean_cancel,B_scan_image_down,Downsample_N,Relative_permittivity1,Relative_permittivity2,R,dx,dt,alpha1,alpha2,n,epsilon_r] = generate_hough_space(png_path,mode,image_path,image_save);
end

%% 读取霍夫空间参数
k = size(hough_space1,3);  % 霍夫空间第三维度的尺寸
epsilon = 2.5:0.1:14;
q = epsilon/(3e8)^2;

%% 顶点搜寻，双曲线绘制
% for i = 1:k
%     [target_area,row,col] = Vertex_Search(hough_space1(:,:,i));
%     col = col(col ~= 0);
%     row = row(row ~= 0);
%     figure;
%     imagesc(target_area);colormap('gray');hold on;   
%     scatter(col,row,20,'r','filled');   
%     title(['figure:',num2str(i),'（epsilon:',num2str(epsilon(i)),'）']);
%     L = 150;
%     comp = 0;
%     for n = 1:length(col)
%         min_x = max(col(n)-L,1);
%         max_x = min(col(n)+L,size(target_area,2));
%         for j = min_x:max_x
%             y = round(  (sqrt(((row(n)*dt*Downsample_N-comp)/2/sqrt(q(k))+R)^2+(j*dx-col(n)*dx)^2)-R) *2*sqrt(q(k))/dt/Downsample_N+comp/dt);
%             if (y > 0) && (y <= size(target_area,1)) 
%                 scatter(j,y,2,'y','filled');
%             end
%         end
%     end
% end

%% 匹配度最高的模板
    % sum_hough_space = sum(hough_space , [1,2]);
    % max_match_q =  find(sum_hough_space == max(sum_hough_space));

    max_para1=max(hough_space1,[],1);
    max_para1=max(max_para1,[],2);
    [~,max_match_q] = max(squeeze(max_para1));
    [target_area,row1,col1] = Vertex_Search(hough_space1(:,:,max_match_q));
    col1 = col1(col1 ~= 0);
    row1 = row1(row1 ~= 0);
    image_show1(col1,row1,B_scan_image_down,B_scan_image_down_Mean_cancel,target_area,max_match_q,dt,Downsample_N,R,dx,q);

    max_para2=max(hough_space2,[],1);
    max_para2=max(max_para2,[],2);
    [~,max_match_q2] = max(squeeze(max_para2));
    [target_area2,row2,col2] = Vertex_Search(hough_space2(:,:,max_match_q2));
    col2 = col2(col2 ~= 0);
    row2 = row2(row2 ~= 0);
    image_show2(col2,row2,B_scan_image_down,B_scan_image_down_Mean_cancel,target_area2,max_match_q2,dt,Downsample_N,R,dx,q);

%% 数据结果展示
    % disp(['花费时间是：',num2str(toc),'s']);
    % disp(['请查阅figure',num2str(max_match_q),'左右的图片']);
    % disp(['下采样阶数：',num2str(Downsample_N)]);
    disp(['第二代双曲线反演的介电常数：',num2str(epsilon(max_match_q))]);
    disp(['第三代双曲线反演的介电常数：',num2str(epsilon(max_match_q2))]);
    % disp(['渐近线线反演的介电常数(膨胀)：',num2str(Relative_permittivity2)]);
    % disp(['渐近线线反演的介电常数(均值对消)：',num2str(Relative_permittivity)]);
    % disp(['alpha1:',num2str(alpha1)]);
    % disp(['alpha2:',num2str(alpha2)]);
    % disp(['厚度反演出的相对介电常数:',num2str(epsilon_r)]);
    % disp(['渐近线线反演的介电常数(data_out)：',num2str(Relative_permittivity1)]);
    % disp(['渐近线线反演的介电常数(非膨胀)：',num2str(Relative_permittivity2)]);
    % disp(['alpha1:',num2str(alpha1)]);
    % disp(['alpha2:',num2str(alpha2)]);

%% 将最佳模板附近的双曲线在霍夫空间中展示出来,进行二次筛选
% % 第二代公式霍夫变换
% nums = 5;
% min_k1 = max(max_match_q-nums,1);
% max_k1 = min(max_match_q+nums,k);
% for i = min_k1:max_k1
%     [target_area1n,row1n,col1n] = Vertex_Search(hough_space1(:,:,i));
%     col1n = col1n(col1n ~= 0);
%     row1n = row1n(row1n ~= 0);
%     figure;
%     imagesc(target_area1n);colormap('gray');hold on;   
%     scatter(col1n,row1n,20,'r','filled');   
%     title(['第二代-figure:',num2str(i),'（epsilon:',num2str(epsilon(i)),'）']);
%     L = 75;
%     comp = 0;
%     for n = 1:length(col1n)
%         min_x = max(col1n(n)-L,1);
%         max_x = min(col1n(n)+L,size(target_area1n,2));
%         for j = min_x:max_x
%             y = round(  (sqrt(((row1n(n)*dt*Downsample_N-comp)/2/sqrt(q(i))+R)^2+(j*dx-col1n(n)*dx)^2)-R) *2*sqrt(q(i))/dt/Downsample_N+comp/dt);
%             if (y > 0) && (y <= size(target_area1n,1)) 
%                 scatter(j,y,2,'y','filled');
%             end
%         end
%     end
% end
% 
% % 第三代公式霍夫变换
% min_k2 = max(max_match_q2-5,1);
% max_k2 = min(max_match_q2+5,k);
% for i = min_k2:max_k2
%     [target_area2n,row2n,col2n] = Vertex_Search(hough_space2(:,:,i));
%     col2n = col2n(col2n ~= 0);
%     row2n = row2n(row2n ~= 0);
%     figure;
%     imagesc(target_area2n);colormap('gray');hold on;   
%     scatter(col2n,row2n,20,'r','filled');   
%     title(['第三代-figure:',num2str(i),'（epsilon:',num2str(epsilon(i)),'）']);
%     L = 75;
%     comp = 0;
%     for n = 1:length(col2n)
%         min_x = max(col2n(n)-L,1);
%         max_x = min(col2n(n)+L,size(target_area2n,2));
%         for j = min_x:max_x
%             distance_io = 0.04 ; % 发射源与接收源的间距
%             v = 1/sqrt(q(i));
% 
%             t1 = 1/v*(  sqrt( ((row2n(n)*dt*Downsample_N-comp)/2*v+R)^2 + (j*dx-col2n(n)*dx)^2 )  -R  );   % 前半程时间
% 
%             sin_a = (row2n(n)*dt*Downsample_N*v/2+R)/(t1*v+R);
%             cos_a = (j*dx-col2n(n)*dx)/(t1*v+R);
% 
%             t2 = 1/v*sqrt(  (t1*v*sin_a)^2  +  (j*dx-col2n(n)*dx+distance_io-R*cos_a)^2  );  %后半程时间
%             y = round(   (t1+t2)/dt/Downsample_N+comp/dt   );
%             if (y > 0) && (y <= size(target_area2n,1)) 
%                 scatter(j,y,2,'y','filled');
%             end
%         end
%     end
% end

%% 将最佳模板附近的双曲线在均值对消图域展示出来,进行二次筛选
% % 第二代公式霍夫变换
% nums = 5;
% min_k1 = max(max_match_q-nums,1);
% max_k1 = min(max_match_q+nums,k);
% for i = min_k1:max_k1
%     figure;
%     imagesc(B_scan_image_down_Mean_cancel);colormap('gray');hold on;   
%     scatter(col1,row1,20,'r','filled');   
%     title(['第二代-figure:',num2str(i),'（epsilon:',num2str(epsilon(i)),'）']);
%     L = 150;
%     comp = 0;
%     for n = 1:length(col1)
%         min_x = max(col1(n)-L,1);
%         max_x = min(col1(n)+L,size(B_scan_image_down_Mean_cancel,2));
%         for j = min_x:max_x
%             y = round(  (sqrt(((row1(n)*dt*Downsample_N-comp)/2/sqrt(q(i))+R)^2+(j*dx-col1(n)*dx)^2)-R) *2*sqrt(q(i))/dt/Downsample_N+comp/dt);
%             if (y > 0) && (y <= size(B_scan_image_down_Mean_cancel,1)) 
%                 scatter(j,y,2,'y','filled');
%             end
%         end
%     end
% end
% 
% % 第三代公式霍夫变换
% min_k2 = max(max_match_q2-5,1);
% max_k2 = min(max_match_q2+5,k);
% for i = min_k2:max_k2
%     figure;
%     imagesc(B_scan_image_down_Mean_cancel);colormap('gray');hold on;   
%     scatter(col2,row2,20,'r','filled');   
%     title(['第三代-figure:',num2str(i),'（epsilon:',num2str(epsilon(i)),'）']);
%     L = 150;
%     comp = 0;
%     for n = 1:length(col2)
%         min_x = max(col2(n)-L,1);
%         max_x = min(col2(n)+L,size(B_scan_image_down_Mean_cancel,2));
%         for j = min_x:max_x
%             distance_io = 0.04 ; % 发射源与接收源的间距
%             v = 1/sqrt(q(i));
% 
%             t1 = 1/v*(  sqrt( ((row2(n)*dt*Downsample_N-comp)/2*v+R)^2 + (j*dx-col2(n)*dx)^2 )  -R  );   % 前半程时间
% 
%             sin_a = (row2(n)*dt*Downsample_N*v/2+R)/(t1*v+R);
%             cos_a = (j*dx-col2(n)*dx)/(t1*v+R);
% 
%             t2 = 1/v*sqrt(  (t1*v*sin_a)^2  +  (j*dx-col2(n)*dx+distance_io-R*cos_a)^2  );  %后半程时间
%             y = round(   (t1+t2)/dt/Downsample_N+comp/dt   );
%             if (y > 0) && (y <= size(B_scan_image_down_Mean_cancel,1)) 
%                 scatter(j,y,2,'y','filled');
%             end
%         end
%     end
% end
% 
% 
