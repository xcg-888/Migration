clear all;
clc;
close all;
%---------------------

%% 输入参数
filename = 'ltefile327_p02.lte';
[TrackInterval,dt,Ydata] = readB_scan(filename);
Ydata(:,1:20) = [];
% 提取数据及文头
t = (1-1:1024-1).*dt.*1e9;
Er = 8; %介质介电常数   
D = t.*(0.3/sqrt(Er))./2; %深度
Ydata = Ydata./max(max(abs(Ydata))); %归一化
[ny, nx] = size(Ydata);
x = (0:nx-1);  dx = 0.5; % 步进，单位[cm]
X = x.*dx; % 距离，单位[cm]
t_delay= 1.157;  %零时校准[ns]
t_max=10; %设置显示的时窗[ns]
D_max=t_max.*(0.3/sqrt(Er))./2;

colorbar_amp=0.1; %显示颜色图尺度
max_gain=1; %显示增益系数，1表示没有显示增益，越大杂波越强

%% 去背景
Ydata = Ydata-repmat(sum(Ydata,2)/nx,1,nx); % 均值对消
%% 零时校正
index_cut_start=find(t>=t_delay,1,'first'); %起始截取下标
Ydata(1:index_cut_start,:)=0;
Ydata=circshift(Ydata,[-1*index_cut_start,0]);

%%
figure(1)
imagesc(X,t,Ydata,[-colorbar_amp colorbar_amp])
% imagesc(x,tsamp,Ydata,[-1000000 1000000])
colormap gray
% colorbar
xlabel('Distance/cm')
% xlabel('Trace')
ylabel('Time/ns')
xlim([0,X(end)])
ylim([0,t_max])
% ylim([3,7.5])
% ylim([0 t_max-t_delay]);
set(gca,'Fontsize',18)
% set(gca,'YDir','reverse');
% set(gca,'xtick',[],'xticklabel',[])
title(num2str(filename))
% title(' ')
set(gcf,'unit','centimeters','position',[2 2 10 10.5])



