%% 图像中直线检测
function [I_out,abs_K,Relative_permittivity]=HoughLineDetect(Img,N,thresholdValue)
%input: 
%       img:输入图像；
%       thresholdValue：hough阈值；
%output: 
%       I_out:检测直线结果，二值图像；
    if  ~exist( 'thresholdValue', 'var' )
        thresholdValue = 20;
    end
    
    if length(size(Img))>2
        I_gray=rgb2gray(Img); %如果输入图像是彩色图，需要转灰度图
    else
        I_gray=Img;
    end
    [rows,cols]=size(I_gray);   %图像大小
    rho_max=floor(sqrt(rows^2+cols^2))+1; %由图像坐标算出ρ最大值，结果取整并加1，作为极坐标系最大值
    AccArray=zeros(rho_max,180);       %初始化极坐标系的数组
    Theta=0:pi/180:pi;                           %定义θ数组，范围从0-180度


    for n=1:rows
        for m=1:cols
            if I_gray(n,m)==1
                for k=10:170
                    %hough变换方程求ρ值
                    rho=(n*cos(Theta(k)))+(m*sin(Theta(k)));
                    %为了防止ρ值出现负数，将ρ值与ρ最大值的和的一半作为ρ的坐标值
                    rho_int=round(rho/2+rho_max/2);
                    %在极坐标中标识点，相同点累加
                    AccArray(rho_int,k)=AccArray(rho_int,k)+1;
                end
            end
        end
    end

    %利用hough变换提取直线
    K=1;                             %存储数组计数器
    % 使用find获取满足条件的rho（行索引）和theta（列索引）
    [case_accarray_rho, case_accarray_theta] = find(AccArray >= thresholdValue);
    % 此时K的值等于找到的元素个数，即 length(case_accarray_rho)
    K = length(case_accarray_rho);
   
    %将直线提取出来,输出图像数组I_out
    I_out=zeros(rows,cols);
    for n=1:rows
        for m=1:cols
             if I_gray(n,m)==1
                 for k=1:180
                    rho=(n*cos(Theta(k)))+(m*sin(Theta(k)));
                    rho_int=round(rho/2+rho_max/2);
                    for a=1:K
                        if rho_int==case_accarray_rho(a)&&k==case_accarray_theta(a)
                            I_out(n,m)=1; 
                        end
                    end
                 end
             end
        end
    end

    [~, max_theta] = find(AccArray == max(max(AccArray)));

    dt = 4.7173e-12;
   
    c = 3e8;
    dx = 0.01;
    abs_K = abs(-cot(max_theta));
    Relative_permittivity = (abs_K*dt*N*c/2/dx).^2;

    % figure,imshow(Img);title('输入图像');
    % figure,imshow(BW);title('edge处理后的边界图');
    % figure,imshow(I_out);title('Hough变换检测出的直线');

    figure,imagesc(I_out);  
    colormap(gray);  
    title('Hough变换检测出的直线');
end