%%%FK3D
clc
clear
%% 鍙傛暟璁剧疆                                                                                                                  
M=89;    %%M is the number of B-Scan
N=89; %%N is the trace number of a B-Scan

%% 读取数据
timeFilePath ='D:\Files\MATLAB\GPR\Migration\3D\time_89.mat'; 
load(timeFilePath);
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
% imag = permute(imag,[1,3,2]);
% imag = [imag;zeros(7000-size(imag,1),89,89)];%%减小频谱栅栏效应
% for kk=1:89
%     imag(:,:,kk) = imag(:,:,kk) - mean(imag(:,:,kk));
% end
% save('F:\project_50\L型0.5_0.5_0.15\RTM_codes\imag_clean','imag');
tic
%% FK成像处理
C = 3e8/sqrt(6);
f0 = 5e9;
Lamda = C/f0;
t_buchang = 1/f0; 
DX=0.005;%%道间距
DY=0.005;%%道间距

SampSize=[0.5,0.001, 0.5,0.001, 0.1,0.001];                              
Xrange = SampSize(1);                                                     % X方向成像场景大小
Yrange = SampSize(3);                                                     % Y方向成像场景大小
Zrange = SampSize(5);                                                     % Z方向成像场景大小
dx = SampSize(2);
dy = SampSize(4);
dz = SampSize(6);

zAxis=dz:dz:0.1;
fSceneZ_array=dz:dz:0.1;                                                  %%成像深度
nSceneZ=length(fSceneZ_array);                                            %%成像场景深度个数

fft_imag = fftn(imag);
omega_tem = 2*pi/(size(imag,1)*1.9258e-12)*([1:floor(size(imag,1)/2),(size(imag,1)-(ceil(size(imag,1)/2)+1:size(imag,1)))]).'*ones(1,N);%1*2598
omega=repmat(omega_tem,1,1,M);

Kx_tem=2*pi/(N*DX)*([1:floor(N/2),(N-(floor(N/2)+1:N))]);%1*89
Kx=repmat(repmat(Kx_tem,size(imag,1),1),1,1,M);

Ky_tem=reshape(2*pi/(N*DY)*([1:floor(M/2),(M-(floor(M/2)+1:M))]),[1,1,M]);%1*89
Ky=repmat(Ky_tem,size(imag,1),N,1);

Kz_tem=2*pi/Zrange*(1:nSceneZ);%1*100
Kz=sqrt(4*omega.^2/C^2- Kx.^2-Ky.^2);
temp = 4*omega.^2/C^2-Kx.^2-Ky.^2;
Kz(floor(size(Kz,1)/2):end,:,:) = -Kz(floor(size(Kz,1)/2):end,:,:);
Kz(1,:)=0;

fft_imag(temp<0) = 0;
fft_imag(4*omega.^2/C^2>((pi/dx)^2+(pi/dy)^2+(pi/dz)^2)) = 0;

conv_fft_imag=zeros(N,M,nSceneZ);%89*89*100
for kx=1:N  %插值运算
    for ky=1:M
        index=floor(size(imag,1)*1.9258e-12*C/(4*pi)*sqrt((Kx_tem(1,kx))^2+(Ky_tem(1,ky))^2+Kz_tem.^2)+1);%索引在后续中访问特定元素
        index(51:end)=-index(50:-1:1)+2444;%调整数据结果，后半部分依赖于前半部分
%         Kz = Kz(:,1:89);index = index(:,1:89);
        ecn=(Kx_tem(1,kx))^2+(Ky_tem(1,ky))^2+Kz_tem.^2;
        mval=C/2*Kz_tem./sqrt((Kx_tem(1,kx))^2+(Ky_tem(1,ky))^2+Kz_tem.^2+eps).*(fft_imag(index,kx,ky)).';%公式推导部分
        conv_fft_imag(kx,ky,:)=reshape(mval,[1,1,nSceneZ]);
    end
end
res=abs(ifftn(conv_fft_imag));
FK_Image=zeros(500,500,100);
[x,y] = meshgrid(0.03:0.005:0.47);
[xs,ys] = meshgrid(dx:dx:0.5);
for deep=1:nSceneZ
    FK_Image(:,:,deep) = interp2(x,y,res(:,:,deep),xs,ys,'cubic');
end
toc
%  FK_Image(isnan(FK_Image)) = 0;


% load E:\project\代码\3DCylinderSceneFK_89try2.mat;
% FK_Image(:,:,1:25)=NaN;
FK_Image=FK_Image/max(max(max(FK_Image)));
FK_Image(FK_Image<0.4)=NaN;
[x,y,z] = meshgrid(0.001:0.001:0.5,0.001:0.001:0.5,0.001:0.001:0.1);
xs = 0.001:0.001:0.5;
ys = 0.001:0.001:0.5;
zs = 0.001:0.001:0.1;
h = slice(x,y,z,FK_Image(:,:,1:1:end),xs,ys,zs);
xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)');
axis([0 ,0.5,0,0.5,0,0.1]);
xlim([0,0.5])
ylim([0,0.5])
zlim([0,0.1])
% set(h,'FaceColor','interp','EdgeColor','none');
set(gca,'ZDir','reverse');
set(gca,'linewidth',1,'fontsize',7.5,'fontname','Times New Roman');
shading interp
camproj perspective;box on;
% set(gcf,'LineWidth',2);
view(45,32);
colormap ('jet');colorbar();
axis equal;grid on;
% saveas(gcf,'E:\project_50\Cylinder_0.06_0.02\21_89\plot_FK3D_cylinder.png');%%
%%FK3D_圆柱体
% toc


% %%%FK3D_L_89-89
% clc
% clear
% %% 鍙傛暟璁剧疆                                                                                                               
% M=89;    %%M is the number of B-Scan
% N=89; %%N is the trace number of a B-Scan
% 
% %% 读取数据
% timeFilePath ='E:\project_50\L型0.5_0.5_0.15\89_89\gen_out89\time_L_Two.mat';
% load(timeFilePath);
% imag = zeros(2598,89,89);
% for m=1:89
%     Data.samp = load(['E:\project_50\L型0.5_0.5_0.15\89_89\gen_out89\L_Two_',num2str(m),'.mat']);
%     Data.samp.data(:,1:end)=Data.samp.data(:,1:end)-mean(Data.samp.data);  % 均值滤波
%     imag(:,1:end,m)=imag(:,1:end,m)-mean(imag(:,:,m));% 均值滤波
%     imag(:,1:size(Data.samp.data,2),m) = Data.samp.data;
% end
% void_imag = zeros(2598,89,89);
% for m=1:89
%     Data.samp = load(['E:\project_50\Cylinder_0.06_0.02\void\void_out1\Cylinder_void_',num2str(m),'.mat']);
%     Data.samp.data(:,1:end)=Data.samp.data(:,1:end)-mean(Data.samp.data);  % 均值滤波
%     void_imag(:,1:end,m)=void_imag(:,1:end,m)-mean(void_imag(:,:,m));% 均值滤波
%     void_imag(:,1:size(Data.samp.data,2),m) = Data.samp.data;
% end
% 
% imag = imag - void_imag;
% 
% t=time;
% tic
% %% FK成像处理
% C = 3e8/sqrt(6);
% f0 = 5e9;
% Lamda = C/f0;
% 
% DX=0.005;%%道间距
% DY=0.005;%%道间距
% 
% SampSize=[0.5,0.001, 0.5,0.001, 0.1,0.001];                              
% Xrange = SampSize(1);                                                     % X方向成像场景大小
% Yrange = SampSize(3);                                                     % Y方向成像场景大小
% Zrange = SampSize(5);                                                     % Z方向成像场景大小
% dx = SampSize(2);
% dy = SampSize(4);
% dz = SampSize(6);
% 
% zAxis=dz:dz:0.1;
% fSceneZ_array=dz:dz:0.1;                                                  %%成像深度
% nSceneZ=length(fSceneZ_array);                                            %%成像场景深度个数
% 
% fft_imag = fftn(imag);
% omega_tem = 2*pi/(size(imag,1)*1.9258e-12)*([1:floor(size(imag,1)/2),(size(imag,1)-(ceil(size(imag,1)/2)+1:size(imag,1)))]).'*ones(1,N);%1*2598
% omega=repmat(omega_tem,1,1,M);
% 
% Kx_tem=2*pi/(N*DX)*([1:floor(N/2),(N-(floor(N/2)+1:N))]);%1*89
% Kx=repmat(repmat(Kx_tem,size(imag,1),1),1,1,M);
% 
% Ky_tem=reshape(2*pi/(N*DY)*([1:floor(M/2),(M-(floor(M/2)+1:M))]),[1,1,M]);%1*89
% Ky=repmat(Ky_tem,size(imag,1),N,1);
% 
% Kz_tem=2*pi/Zrange*(1:nSceneZ);%1*100
% Kz=sqrt(4*omega.^2/C^2- Kx.^2-Ky.^2);
% temp = 4*omega.^2/C^2-Kx.^2-Ky.^2;
% Kz(floor(size(Kz,1)/2):end,:,:) = -Kz(floor(size(Kz,1)/2):end,:,:);
% Kz(1,:)=0;
% 
% fft_imag(temp<0) = 0;
% fft_imag(4*omega.^2/C^2>((pi/dx)^2+(pi/dy)^2+(pi/dz)^2)) = 0;
% fft_imag=fft_imag.*exp(1i*Kz*0.1);
% conv_fft_imag=zeros(N,M,nSceneZ);%89*89*100
% for kx=1:N  %插值运算
%     for ky=1:M 
%         index=floor(size(imag,1)*1.9258e-12*C/(4*pi)*sqrt((Kx_tem(1,kx))^2+(Ky_tem(1,ky))^2+Kz_tem.^2)+1);
%         index(51:end)=-index(50:-1:1)+2444;
% %         Kz = Kz(:,1:89);index = index(:,1:89);
%         mval=C/2*Kz_tem./sqrt((Kx_tem(1,kx))^2+(Ky_tem(1,ky))^2+Kz_tem.^2+eps).*(fft_imag(index,kx,ky)).';
%         conv_fft_imag(kx,ky,:)=reshape(mval,[1,1,nSceneZ]);
%     end
% end
% res=abs(ifftn(conv_fft_imag));
% FK_Image=zeros(500,500,100);
% [x,y] = meshgrid(0.03:0.005:0.47);
% [xs,ys] = meshgrid(dx:dx:0.5);
% for deep=1:nSceneZ
%     FK_Image(:,:,deep) = interp2(x,y,res(:,:,deep),xs,ys,'cubic');
% %     imagesc(FK_Image(:,:,deep));
% %     colorbar
% end
% toc
% % save(strcat(strcat('3D',flag),'FK_89_cube_try1.mat'),'FK_Image');
% 
% 
% 
% % load E:\project\代码\3DL_FK_89_cube_try1.mat;
% % FK_Image(:,:,1:40)=NaN;
% FK_Image=FK_Image/max(max(max(FK_Image)));
% FK_Image(FK_Image<0.3)=NaN;
% 
% [x,y,z] = meshgrid(0.001:0.001:0.5,0.001:0.001:0.5,0.001:0.001:0.1);
% xs = 0.001:0.001:0.5;
% ys = 0.001:0.001:0.5;
% zs = 0.001:0.001:0.1;
% h = slice(x,y,z,FK_Image(:,:,end:-1:1),xs,ys,zs);
% xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)');
% axis([0 ,0.5,0,0.5,0,0.1]);
% set(h,'FaceColor','interp','EdgeColor','none');
% set(gca,'linewidth',1,'fontsize',12,'fontname','Times New Roman');
% set(gca,'ZDir','reverse')
% camproj perspective;box on;
% view(12,20);
% colormap ('jet');colorbar();
% axis equal;grid on;
% saveas(gcf,'E:\project_50\L型0.5_0.5_0.15\89_89\plot_FK3D_L.png');%%