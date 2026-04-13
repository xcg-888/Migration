% 本程序用于处理一维(单条测线)数据。
% 使用时先修改'FolderName'变量的值为存放原始数据的文件夹名，
% 再修改基本参数（“参数输入”小节部分），然后运行即可。
% 本程序对原始数据进行了去背景、带通滤波、零时校正处理，最终显示单条测线的剖面图。
% clear all,clc,close all
function [data_t, dx, dt] = readVectorNet_data()
%% 参数输入
[name, path]=uigetfile;
FolderName = path;
fmin=0.1;  %数据采集起始频率[GHz]
fmax=2.9;     %数据采集终止频率[GHz]
samples=501;%频率采样点数

dx=0.025;    %道间距[m]
traces=63; %数据道数
x=dx*(1:traces);

colorbar_amp=0.02; %显示颜色图尺度，建议取值0~0.1之间
max_gain=1; %显示增益系数，1表示没有显示增益，越大杂波越强

f_low_stop=fmin;  % 带通滤波器的第1个频率点[GHz]
f_low_pass=fmin+0.8;    % 带通滤波器的第2个频率点[GHz]
f_high_pass=fmax-0.4; % 带通滤波器的第3个频率点[GHz]
f_high_stop=fmax;   % 带通滤波器的第4个频率点[GHz]

fftcoef=16;      % ifft点数因子,,值越大则时域点数越多

t_delay=0;  %零时校准[ns]
t_max=30; %设置显示的时窗[ns]
%% 数据读取
data_f=zeros(samples,traces);
for k=1:traces
    str= strcat (FolderName,'/',int2str(k) , '.txt') ;
    datatemp=load(str);
    for m=2:samples+1
        data_f(m-1,k)=datatemp(m,1)+1j*datatemp(m,2);
    end
end
%% 带通滤波 Band pass filter
f=linspace(fmin,fmax,samples);
df=f(2)-f(1);%频域采样间隔
[win]=mm_hamming(f_low_stop,f_low_pass,f_high_pass,f_high_stop,fmin,df,samples);
win=win';
win=repmat(win,[1,traces]);
data_f=data_f.*win;
%% 逆傅里叶变换 ifft
nop=2^(nextpow2(samples))*fftcoef; %ifft点数
data_t=real(ifft(data_f,nop));
dt=1/(nop*df);         %时间间隔
t=(0:(nop-1))*dt;  %时域时间轴
%% Gain
g=(1:nop).';
g(round(max_gain/dt):end)=g(round(max_gain/dt));
data_t=data_t.*repmat(g,1,traces);
%% 零时校正
index_cut_start=find(t>=t_delay,1,'first'); %起始截取下标
data_t(1:index_cut_start,:)=0;
data_t=circshift(data_t,[-1*index_cut_start,0]);
data_t = data_t(1:round(t_max./dt),:);
end
%% 显示剖面图
% x=dx*(1:traces);
% figure,
% subplot(1,2,1)
% set(gca,'fontsize',16); 
% imagesc(x,t,data_t,[-colorbar_amp colorbar_amp]);set(gca,'YDir','reverse');
% colormap(gray);colorbar;colorbar('FontSize',13);
% ylim([0 t_max-t_delay]);
% xlabel('Position [m]','fontsize',18);set(gca,'fontsize',18);
% ylabel('Time [ns]','fontsize',18);set(gca,'fontsize',18);
% 
% % figure,
% % plot(t(1:1000),data_t(1:1000,413))
% subplot(1,2,2)
% % figure,
% set(gca,'fontsize',16); 
% imagesc(x,t,data_t-repmat(mean(data_t,2),1,traces),[-colorbar_amp colorbar_amp]);set(gca,'YDir','reverse');
% colormap(gray);colorbar;colorbar('FontSize',13);
% ylim([0 t_max-t_delay]);
% xlabel('Position [m]','fontsize',18);set(gca,'fontsize',18);
% ylabel('Time [ns]','fontsize',18);set(gca,'fontsize',18);
% set(gcf,'unit','centimeters','position',[2 2 30 10])