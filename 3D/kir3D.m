%% 参数设置                                                                %%可以为Ex，Ey，Ez
flag='Cylinder';                                                          %%读取三维场景数据，可以为SingleBar、SteelMesh
M=89;    %%M is the number of B-Scan
N=89; %%N is the trace number of a B-Scan

%% 天线接收BScan数据读取"/WorkSpaceP2/SXT/SinglePointScene0_/web_change/web_small/time.mat"
timeFilePath = 'D:\Files\MATLAB\GPR\Migration\3D\time_89.mat';  
load(timeFilePath);
imag = zeros(89,89,2598);
for m=1:89
    Data.samp = load(['D:\Files\MATLAB\GPR\Migration\3D\gen_out89\L_Two_',num2str(m),'.mat']);
    Data.samp.data(1:end,:)=Data.samp.data(1:end,:)-mean(Data.samp.data);  % 均值滤波
    imag(m,1:end,:)=imag(m,1:end,:)-mean(imag(m,:,:));% 均值滤波
    imag(m,1:size(Data.samp.data, 2), :) = Data.samp.data';
end
t=time;
%% 不同深度成像处理    
Epsilon=6;
C=3e8/sqrt(Epsilon);
%     BP成像处理
%     天线网格位置
x_line=(0.03:0.005:0.47).'*ones(1,N);
y_line_s=ones(M,1)*(0.03:0.005:0.47);
y_line_r=ones(M,1)*(0.035:0.005:0.475);
z_line=0.1*ones(M,N);

dx = 0.001;
dy = 0.001;
xAxis = 0.001:dx:0.5;%%
yAxis = 0.001:dy:0.5;%%成像区间
xMatrix = ones(length(yAxis),1)*xAxis;
yMatrix = yAxis'*ones(1,length(xAxis));

fSceneZ_array=0.001:0.001:0.1;                                               %%成像深度序列
nSceneZ=length(fSceneZ_array);                                             %%成像场景深度个数
KircImage = zeros(length(yAxis),length(xAxis),nSceneZ);
for deep=1:nSceneZ
    zMatrix = fSceneZ_array(deep)*ones(length(yAxis),length(xAxis));
    for i1=1:M
        for i2 = 1:N
            vector_X = xMatrix-x_line(i1,i2);
            vector_Y_s = yMatrix-y_line_s(i1,i2);
            vector_Y_r = yMatrix-y_line_r(i1,i2);
            vector_Z = zMatrix-z_line(i1,i2);
            %%向量的模值
            rangeMatrix_s = sqrt(vector_X.^2+vector_Y_s.^2+vector_Z.^2);
            rangeMatrix_r = sqrt(vector_X.^2+vector_Y_r.^2+vector_Z.^2);
            [row,column]=find(rangeMatrix_s==0);
            for j=1:length(row)
            rangeMatrix_s(row(j),column(j))=eps;
            end
            %%数据平面法线方向
            vector_S = [0,0,-1];
            %%两个向量之间的夹角余弦值
            cos_theta = (vector_X*vector_S(1)+vector_Y_s*vector_S(2)+vector_Z*vector_S(3))./((rangeMatrix_s+rangeMatrix_r)/2)/(sqrt(vector_S*vector_S.')); 
           signalTemp =  imag(i1,i2,:);
           signalTemp1 =  diff(imag(i1,i2,:))/1.9258e-12;
           nRangeCell = round((rangeMatrix_r+rangeMatrix_s)/(C*1.9258e-12));
           nRangeCell(nRangeCell<1) = 1;
           nRangeCell(nRangeCell>size(imag , 3)-1) = size(imag , 3)-1;
           KircImage(:,:,deep) = KircImage(:,:,deep)+(cos_theta./rangeMatrix_s/C).*signalTemp1(nRangeCell)+...
               (cos_theta./(rangeMatrix_r.^2)).*signalTemp(nRangeCell); %KircImage(:,:,deep):第一维与第二维的全部数据，第三维的deep数据
        end
    end
end

% save(strcat(strcat('3D',flag),'_SceneKir_89_cube.mat'),'KircImage');

% load E:\project_50\L型0.5_0.5_0.15\21_89\3DCylinder_SceneKir_L_Two.mat;
% imageField=readNPY('D:\postgraduation\偏移成像\成像_副本\python_code\imageField.npy');
% imageField=imageField/max(max(max(imageField)));
% imageField(imageField<0.2)=NaN;
KircImage=KircImage/max(max(max(KircImage)));
KircImage(KircImage<0.8)=NaN;
[x,y,z] = meshgrid(0.001:0.001:0.5,0.001:0.001:0.5,0.001:0.001:0.1);
xs = 0.001:0.001:0.5;
ys = 0.001:0.001:0.5;
zs = 0.001:0.001:0.1;
h = slice(x,y,z,KircImage(:,:,end:-1:1),xs,ys,zs);
xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)');
axis([0,0.5,0,0.5,0,0.1]);
set(h,'FaceColor','interp','EdgeColor','none');

set(gca,'linewidth',1,'fontsize',12,'fontname','Times New Roman');

camproj perspective;box on;
view(12,20);
colormap ('jet');colorbar();
axis equal;grid on;