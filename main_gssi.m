function [dx, dt,Ydata] = main_gssi(file)

clc;
%-------------
%% 输入参数
filename = file;
Data = readgssi(filename);
% 提取数据及文头
Ydata = Data.samp;  

Ydata = Ydata./max(max(abs(Ydata))); %归一化
Tmax =  Data.head.range*1e-9; %最大时间，Data.head.range：在P_1中是整个时间剖面的大小
% Tmax =  Data.head.range*1e-9+Data.head.position*1e-9;
[ny, nx] = size(Ydata);
dt = Tmax/ny; %
t = (1-1:ny-1).*dt.*1e9; %时间
Er = 1; %介质介电常数   
D = t.*(0.3/sqrt(Er))./2; %深度
x = (0:nx-1);  dx = 0.5; % 步进，单位[cm]
X = x.*dx; % 距离，单位[cm]

gain_p = 0;  % 是否增益，1：增益，0：不增益
gain_t=5.*1e-9; % 指数增益时窗[ns]
t_delay= 0.2504;  %零时校准[ns]
t_max=3; %设置显示的时窗[ns]
D_max=t_max.*(0.3/sqrt(Er))./2;

colorbar_amp=0.5; %显示颜色图尺度
max_gain=1; %显示增益系数，1表示没有显示增益，越大杂波越强

ntrace = 16; %查看某道波形
%% 去背景
% Ydata = Ydata-repmat(sum(Ydata,2)/nx,1,nx); % 均值对消
% nn = 50; % 滑动窗口
% [Ydata, nn] = rmswbackgr(Ydata ,nn); % 滑动去背景
%% 零时校正
index_cut_start=find(t>=t_delay,1,'first'); %起始截取下标
Ydata(1:index_cut_start,:)=0;
Ydata=circshift(Ydata,[-1*index_cut_start,0]);

%%
% figure(1)
% imagesc(X,t,Ydata,[-colorbar_amp colorbar_amp])
% colormap gray
% xlabel('Distance/cm')
% ylabel('Time/ns')
% ylim([0,t_max])
% set(gca,'Fontsize',18)
% title(num2str(filename(1:end-8)))
% set(gcf,'unit','centimeters','position',[2 2 10 10.5])
% 
% figure(2)
% imagesc(X,D, Ydata,[-colorbar_amp colorbar_amp])
% colormap gray
% xlabel('Distance/cm')
% ylabel('Depth/m')
% ylim([0,D_max])
% set(gca,'Fontsize',18)
% title(num2str(filename(1:end-8)))
% set(gcf,'unit','centimeters','position',[2 2 10 10.5])


% figure(3)
% plot(t,Ydata(:,ntrace))
% xlabel('Time/ns')
% ylabel('Amplitude')
% xlim([0,t_max])
% 
% set(gca,'Fontsize',18)

end


% % % function [dx, dt,Ydata] = main_gssi(file)
% % % % clear all;
% % % clc;
% % % close all;
% % % %-------------
% % % %% 输入参数
% % % filename = file;
% % % Data = readgssi(filename);
% % % % 提取数据及文头
% % % Ydata = Data.samp;  
% % % % Ydata(1034:end,:) = [];
% % % Ydata = Ydata./max(max(abs(Ydata))); %归一化
% % % Tmax =  Data.head.range*1e-9; %最大时间，Data.head.range：在P_1中是整个时间剖面的大小
% % % % Tmax =  Data.head.range*1e-9+Data.head.position*1e-9;
% % % [ny, nx] = size(Ydata);
% % % dt = Tmax/ny; %
% % % t = (1-1:ny-1).*dt.*1e9; %时间
% % % Er = 1; %介质介电常数   
% % % D = t.*(0.3/sqrt(Er))./2; %深度
% % % x = (0:nx-1);  dx = 0.5; % 步进，单位[cm]
% % % X = x.*dx; % 距离，单位[cm]
% % % 
% % % gain_p = 0;  % 是否增益，1：增益，0：不增益
% % % gain_t=5.*1e-9; % 指数增益时窗[ns]
% % % t_delay= 0.2504;  %零时校准[ns]
% % % t_max=3; %设置显示的时窗[ns]
% % % D_max=t_max.*(0.3/sqrt(Er))./2;
% % % 
% % % colorbar_amp=0.5; %显示颜色图尺度
% % % max_gain=1; %显示增益系数，1表示没有显示增益，越大杂波越强
% % % 
% % % ntrace = 16; %查看某道波形
% % % %% 去背景
% % % % Ydata = Ydata-repmat(sum(Ydata,2)/nx,1,nx); % 均值对消
% % % % nn = 50; % 滑动窗口
% % % % [Ydata, nn] = rmswbackgr(Ydata ,nn); % 滑动去背景
% % % %% 零时校正
% % % index_cut_start=find(t>=t_delay,1,'first'); %起始截取下标
% % % Ydata(1:index_cut_start,:)=0;
% % % Ydata=circshift(Ydata,[-1*index_cut_start,0]);
% % % %% Gain
% % % % g=(1:ny).';
% % % % g(round(gain_t/dt):end)=g(round(gain_t/dt));
% % % % % Gend=round(gain_t/dt);
% % % % % g(1:Gend)=exp(t(1:Gend)*0.008).*(t(1:Gend)+dt);
% % % % % g(Gend+1:end)=g(Gend+1:end)*g(Gend);
% % % %   if gain_p == 0
% % % %      g(:) = 1;
% % % %   end
% % % % Ydata=Ydata.*repmat(g,1,nx);
% % % % Ydata = Ydata.*g; % 
% % % % figure 
% % % % tr = 10;
% % % % plot(Ydata1(:,tr),'r')
% % % % hold on
% % % % plot(Ydata(:,tr),'--')
% % % % figure 
% % % % plot(g(:,tr),'--')
% % % %%
% % % figure(1)
% % % imagesc(X,t,Ydata,[-colorbar_amp colorbar_amp])
% % % % imagesc(x,tsamp,Ydata,[-1000000 1000000])
% % % colormap gray
% % % % colorbar
% % % xlabel('Distance/cm')
% % % % xlabel('Trace')
% % % ylabel('Time/ns')
% % % % xlim([0,X(end)])
% % % ylim([0,t_max])
% % % % ylim([3,7.5])
% % % % ylim([0 t_max-t_delay]);
% % % set(gca,'Fontsize',18)
% % % % set(gca,'YDir','reverse');
% % % % set(gca,'xtick',[],'xticklabel',[])
% % % title(num2str(filename(1:end-8)))
% % % % title(' ')
% % % set(gcf,'unit','centimeters','position',[2 2 10 10.5])
% % % figure(2)
% % % imagesc(X,D, Ydata,[-colorbar_amp colorbar_amp])
% % % colormap gray
% % % % colorbar
% % % xlabel('Distance/cm')
% % % % xlabel('Trace')
% % % ylabel('Depth/m')
% % % % xlim([0,X(end)])
% % % ylim([0,D_max])
% % % % ylim([3,7.5])
% % % % ylim([0 t_max-t_delay]);
% % % set(gca,'Fontsize',18)
% % % % set(gca,'YDir','reverse');
% % % % set(gca,'xtick',[],'xticklabel',[])
% % % title(num2str(filename(1:end-8)))
% % % % title(' ')
% % % set(gcf,'unit','centimeters','position',[2 2 10 10.5])
% % % 
% % % 
% % % figure(3)
% % % 
% % % plot(t,Ydata(:,ntrace))
% % % xlabel('Time/ns')
% % % ylabel('Amplitude')
% % % xlim([0,t_max])
% % % % ylim([-4e6,6e6])
% % % set(gca,'Fontsize',18)
% % % % view(90,90)
% % % % title('Trace25')
% % % % title(['trace:',num2str(ntrace)])
% % % % set(gcf,'unit','centimeters','position',[2 2 12 10])
% % % % legend('正向','反向')% 
% % % end
