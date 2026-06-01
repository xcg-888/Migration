%%PS3D
% clc;
% clear all;
%% 参数设置                                                              
flag='Cylinder';                                                      
M=89; %%M is the number of B-Scan
N=89; %%N is the trace number of a B-Scan

%% 天线接收BScan数据读取
timeFilePath ='D:\Files\MATLAB\GPR\Migration\3D\time_89.mat'; 
load(timeFilePath);
%% 读取数据
% timeFilePath ='E:\project\Cylinder_0.06_0.02\21_89\time_89.mat'; 
% load(timeFilePath);
imag = zeros(2598,89,89);
for m=1:89
    Data.samp = load(['D:\Files\MATLAB\GPR\Migration\3D\gen_out89\L_Two_',num2str(m),'.mat']);
%     Data.samp.data(:,1:end)=Data.samp.data(:,1:end)-mean(Data.samp.data);  % 均值滤波
%     imag(:,1:end,m)=imag(:,1:end,m)-mean(imag(:,:,m));% 均值滤波
    imag(:,1:size(Data.samp.data,2),m) = Data.samp.data;
end

void_imag = zeros(2598,89,89);
for m=1:89
    Data.samp = load(['D:\Files\MATLAB\GPR\Migration\3D\void_out1\Cylinder_void_',num2str(m),'.mat']);
%     Data.samp.data(:,1:end)=Data.samp.data(:,1:end)-mean(Data.samp.data);  % 均值滤波
%     void_imag(:,1:end,m)=void_imag(:,1:end,m)-mean(void_imag(:,:,m));% 均值滤波
    void_imag(:,1:size(Data.samp.data,2),m) = Data.samp.data;
end

imag = imag - void_imag;
% imag = [imag;zeros(10001-size(imag,1),89,89)];%%减小频谱栅栏效应
t=time;
tic
%% PS算法处理
C = 3e8/sqrt(6);
f0 = 5e9;
Lamda = C/f0;

DX=0.005;%%道间距
DY=0.005;%%道间距

SampSize=[0.5,0.001, 0.5,0.001, 0.1,0.001];                               %天线扫面区域大小
Xrange = SampSize(1);                                                     % X方向成像场景大小
Yrange = SampSize(3);                                                     % Y方向成像场景大小
Zrange = SampSize(5);                                                     % Z方向成像场景大小
dx = SampSize(2);
dy = SampSize(4);
dz = SampSize(6);

fft_imag = fftn(imag);%三维傅里叶变换   x,y,t傅里叶变换
omega_tem = 2*pi/(size(imag,1)*1.9258e-12)*([1:floor(size(imag,1)/2),(size(imag,1)-(ceil(size(imag,1)/2)+1:size(imag,1)))]).'*ones(1,N);
omega=repmat(omega_tem,1,1,M);

Kx_tem=2*pi/(N*DX)*([1:floor(N/2),(N-(floor(N/2)+1:N))]);
Kx=repmat(repmat(Kx_tem,size(imag,1),1),1,1,M);
Ky_tem=reshape(2*pi/(N*DY)*([1:floor(M/2),(M-(floor(M/2)+1:M))]),[1,1,M]);
Ky=repmat(Ky_tem,size(imag,1),N,1);

Kz = sqrt(4*omega.^2/C^2-Kx.^2-Ky.^2);
temp = 4*omega.^2/C^2-Kx.^2-Ky.^2;

Kz(floor(size(Kz,1)/2):end,:,:) = -Kz(floor(size(Kz,1)/2):end,:,:);
Kz(1,:)=0;
% Kz=reshape(Kz,[89,89,2598]);

fft_imag(temp<0) = 0;
fft_imag(4*omega.^2/C^2>((pi/dx)^2+(pi/dy)^2+(pi/dz)^2)) = 0;

fSceneZ_array=dz:dz:0.1;                                                
nSceneZ=length(fSceneZ_array);                                          
PS_Image=zeros(500,500,100);

[x,y]=meshgrid(0.03:0.005:0.47);
[xs,ys] = meshgrid(dx:dx:0.5);

for deep=1:nSceneZ   %遍历每个深度z
    Res=fft_imag.*exp(1i*Kz*fSceneZ_array(deep));%根据当前深度fSceneZ_array(deep)计算相位因子exp(1iKzfSceneZ_array(deep))并将其应用于imag
%     Res=reshape(Res,[2598,89,89]);
    a=sum(Res);
    temp=reshape(a,N,M).';
%     temp=reshape(sum(Res),N,M).';
    temp_ifft=ifft2(temp);
    PS_Image(:,:,deep) = interp2(x,y,abs(temp_ifft),xs,ys,'cubic');%沿着kx、ky作2D_FFT   使用插值方法将temp插值到指定的网格点上
end
toc
PS_Image=permute(PS_Image,[2,1,3]);
% save(strcat(strcat('3D',flag),'ScenePS_89_cube.mat'),'PS_Image');

%  PS_Image(isnan(PS_Image)) = 0;

%%
% load E:\project\Cylinder_0.06_0.02\21_89\3DCylinderScenePS_89_cube.mat;
% PS_Image(:,:,1:20)=NaN;
PS_Image=PS_Image/max(max(max(PS_Image)));
PS_Image(PS_Image<0.4)=NaN;

[x,y,z] = meshgrid(0.001:0.001:0.5,0.001:0.001:0.5,0.001:0.001:0.1);
xs = 0.001:0.001:0.5;
ys = 0.001:0.001:0.5;
zs = 0.001:0.001:0.1;
h = slice(x,y,z,PS_Image(:,:,1:1:end),xs,ys,zs);
xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)');
axis([0 ,0.5,0,0.5,0,0.1]);
% set(h,'FaceColor','interp','EdgeColor','none');
% set(gcf,'unit','centimeters','position',[2,2,12,6]);
set(gca,'linewidth',1,'fontsize',7.5,'fontname','Times New Roman');
set(gca,'ZDir','reverse');
view(45,32);
shading interp
camproj perspective;box on;
colormap ('jet');colorbar();
axis equal;grid on;
% saveas(gcf,'plot_PS3D_cylinder_2cm-89-89-equal.png');%%
%%圆柱体21-89结束
% toc

% % %%L形开始_PS
% clc;
% clear all;
% %% 参数设置                                                              
% flag='L_';                                                      
% M=89; %%M is the number of B-Scan
% N=89; %%N is the trace number of a B-Scan
% 
% %% 天线接收BScan数据读取
% timeFilePath = 'E:\project\L型0.5_0.5_0.15\89_89\gen_out89\time_L_Two.mat'; 
% load(timeFilePath);
% imag = zeros(2598,89,89);
% for m=1:89
%     Data.samp = load(['E:\project\L型0.5_0.5_0.15\89_89\gen_out89\L_Two_',num2str(m),'.mat']);
%     Data.samp.data(:,1:end)=Data.samp.data(:,1:end)-mean(Data.samp.data);  % 均值滤波
%     imag(:,m,1:end)=imag(:,m,1:end)-mean(imag(:,m,:));% 均值滤波
%     imag(:,m,1:size(Data.samp.data,2)) = Data.samp.data;
% end
% t=time;
% tic
% %% PS算法处理
% C = 3e8/sqrt(6);
% f0 = 5e9;
% Lamda = C/f0;
% 
% DX=0.005;%%道间距
% DY=0.005;%%道间距
% 
% SampSize=[0.5,0.001, 0.5,0.001, 0.1,0.001];                               %天线扫面区域大小
% Xrange = SampSize(1);                                                     % X方向成像场景大小
% Yrange = SampSize(3);                                                     % Y方向成像场景大小
% Zrange = SampSize(5);                                                     % Z方向成像场景大小
% dx = SampSize(2);
% dy = SampSize(4);
% dz = SampSize(6);
% 
% fft_imag = fftn(imag);%三维傅里叶变换   x,y,t傅里叶变换
% % d=ceil(size(imag,3)/2);
% % b=ceil(size(imag,3)/2)+1;
% % c=(size(imag,3)-(b:size(imag,3)));
% % omega_tem = 2*pi/(size(imag,3)*1.9258e-12)*([1:d,c]).'*ones(1,N);
% omega_tem = 2*pi/(size(imag,1)*1.9258e-12)*([1:floor(size(imag,1)/2),(size(imag,1)-(ceil(size(imag,1)/2)+1:size(imag,1)))]).'*ones(1,N);
% omega=repmat(omega_tem,1,1,M);
% 
% Kx_tem=2*pi/(N*DX)*([1:floor(N/2),(N-(floor(N/2)+1:N))]);%1*89
% Kx=repmat(repmat(Kx_tem,size(imag,1),1),1,1,M);%2598*89*89
% % e=(M-(floor(M/2)+1:M));
% % f=[1:floor(M/2),(M-(floor(M/2)+1:M))];
% % g=2*pi/(N*DY)*([1:floor(M/2),(M-(floor(M/2)+1:M))])';
% Ky_tem=reshape(2*pi/(N*DY)*([1:floor(M/2),(M-(floor(M/2)+1:M))]),[1,1,M]);%1*1*89
% Ky=repmat(Ky_tem,size(imag,1),N,1);%2598*89*89
% 
% Kz = sqrt(4*omega.^2/C^2-Kx.^2-Ky.^2);
% temp = 4*omega.^2/C^2-Kx.^2-Ky.^2;
% 
% Kz(floor(size(Kz,1)/2):end,:,:) = -Kz(floor(size(Kz,1)/2):end,:,:);
% Kz(1,:)=0;
% % Kz=reshape(Kz,[89,89,2598]);
% 
% fft_imag(temp<0) = 0;
% fft_imag(4*omega.^2/C^2>((pi/dx)^2+(pi/dy)^2+(pi/dz)^2)) = 0;
% 
% fSceneZ_array=dz:dz:0.1;                                                
% nSceneZ=length(fSceneZ_array);                                          
% PS_Image=zeros(500,500,100);
% 
% [x,y]=meshgrid(0.03:0.005:0.47);
% [xs,ys] = meshgrid(dx:dx:0.5);
% 
% for deep=1:nSceneZ   %遍历每个深度z
%     Res=fft_imag.*exp(1i*Kz*fSceneZ_array(deep));%根据当前深度fSceneZ_array(deep)计算相位因子exp(1iKzfSceneZ_array(deep))并将其应用于imag
% %     Res=reshape(Res,[2598,89,89]);
%     a=sum(Res);
%     temp=reshape(a,N,M).';
% %     temp=reshape(sum(Res),N,M).';
%     temp_ifft=ifft2(temp);
%     PS_Image(:,:,deep) = interp2(x,y,abs(temp_ifft),xs,ys,'cubic');%沿着kx、ky作2D_FFT   使用插值方法将temp插值到指定的网格点上
% end
% toc
% % save(strcat(strcat('3D',flag),'PS_89_cube.mat'),'PS_Image');
% 
% % %%
% % clear
% % load E:\project\L型0.5_0.5_0.15\89_89\3DL_PS_89_cube.mat;
% PS_Image(:,:,1:30)=NaN;
% PS_Image=PS_Image/max(max(max(PS_Image)));
% PS_Image(PS_Image<0.4)=NaN;
% 
% [x,y,z] = meshgrid(0.001:0.001:0.5,0.001:0.001:0.5,0.001:0.001:0.1);
% xs = 0.001:0.001:0.5;
% ys = 0.001:0.001:0.5;
% zs = 0.001:0.001:0.1;
% h = slice(x,y,z,PS_Image(:,:,end:-1:1),xs,ys,zs);
% xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)');
% axis([0 ,0.5,0,0.5,0,0.1]);
% set(h,'FaceColor','interp','EdgeColor','none');
% % set(gcf,'unit','centimeters','position',[2,2,12,6]);
% set(gca,'linewidth',1,'fontsize',12,'fontname','Times New Roman');
% 
% % camproj perspective;box on;
% view(12,20);
% colormap ('jet');colorbar();
% axis equal;grid on;
% % saveas(gcf,'plot_PS3D_L_89-89-equal.png');
% %%L形89-89结束