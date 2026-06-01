clc
clear
%%%L型-21-89
%% 参数设置                                                                %%可以为Ex，Ey，Ez
flag='L_Two';                                                          %%读取三维场景数据，可以为SingleBar、SteelMesh
M=21; %%M is the number of B-Scan
N=89; %%N is the trace number of a B-Scan

%% 天线接收BScan数据读取
% timeFilePath = 'E:\project\L型0.5_0.5_0.15\form_out\time_L_Two.mat';  
% load(timeFilePath);
imag = zeros(21,89,2598);
for m=35:55
    Data.samp = load(['F:\Layered_RM\order\single_layer\3D\form_out\L_Two_',num2str(m),'.mat']);
    Data.samp.data(1:end,:)=Data.samp.data(1:end,:)-mean(Data.samp.data);  % 均值滤波
    imag(m,1:size(Data.samp.data, 2), :) = Data.samp.data';
end
% t=time;

tic
%% 不同深度成像处理    
Epsilon=6;
C=3e8/sqrt(Epsilon);
%     BP成像处理
%     天线网格位置
x_line=(0.2:0.005:0.3).'*ones(1,N);
y_line_s=ones(M,1)*(0.03:0.005:0.47);
y_line_r=ones(M,1)*(0.035:0.005:0.475);
z_line=0.1*ones(M,N);

dx = 0.001;
dy = 0.001;
xAxis = 0.001:dx:0.5;
yAxis = 0.001:dy:0.5;%%成像区间
xMatrix = ones(length(yAxis),1)*xAxis;
yMatrix = yAxis'*ones(1,length(xAxis));

fSceneZ_array=0.001:0.001:0.1;                                               %%成像深度序列
nSceneZ=length(fSceneZ_array);                                             %%成像场景深度个数
bpImage = zeros(length(yAxis),length(xAxis),nSceneZ);
for deep=1:nSceneZ
    zMatrix = fSceneZ_array(deep)*ones(length(yAxis),length(xAxis));
    for i1=1:M    %每一个通道(B_scan)的N个点
        for i2 = 1:N
           signalTemp =  imag(i1,i2,:);%获取当前天线的位置处的所有迭代次数
           rangeMatrix_s = sqrt((xMatrix-x_line(i1,i2)).^2+(yMatrix-y_line_s(i1,i2)).^2+(zMatrix-z_line(i1,i2)).^2);%天线到每一点的距离
           rangeMatrix_r = sqrt((xMatrix-x_line(i1,i2)).^2+(yMatrix-y_line_r(i1,i2)).^2+(zMatrix-z_line(i1,i2)).^2);
           nRangeCell = round((rangeMatrix_s+rangeMatrix_r)/(C*1.9258e-12));%距离转化为时间点的索引
           nRangeCell(nRangeCell<1) = 1;
           nRangeCell(nRangeCell>size(imag, 3)) = size(imag, 3);
           bpImage(:,:,deep) = bpImage(:,:,deep)+signalTemp(nRangeCell);%.*exp(1j*4*pi*rangeMatrix/Lamda)
        end
    end
end
% bpImage=abs(bpImage)/max(max(max(abs(bpImage))));
% toc
% save(strcat(strcat('3D',flag),'_SceneBp_L_21-89.mat'),'bpImage');

% %%
% load E:\project_50\L型0.5_0.5_0.15\21_89\3DL_Two_Bp_L_Two.mat;
% bpImage(:,:,90:end)=NaN;
bpImage=bpImage/max(max(max(bpImage)));
bpImage(bpImage<0.88)=NaN;

[x,y,z] = meshgrid(0.001:0.001:0.5,0.001:0.001:0.5,0.001:0.001:0.1);
xs = 0.001:0.001:0.5;
ys = 0.001:0.001:0.5;
zs = 0.001:0.001:0.1;
h = slice(x,y,z,bpImage(:,:,end:-1:1),xs,ys,zs);
xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)');
axis([0,0.5,0,0.5,0,0.1]);
set(h,'FaceColor','interp','EdgeColor','none');
set(gca,'linewidth',1,'fontsize',12,'fontname','Times New Roman');
set(gca,'ZDir','reverse')
camproj perspective;box on;
view(12,20);
colormap ('jet');colorbar();
axis equal;grid on;
% saveas(gcf,'plot_BP3D_cylinder_2cm-89-89.png');%%