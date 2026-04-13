%读取LTD雷达实测原始B扫图像函数
%输入：格式为LTE的文件路径
%输出：B-san 图像原始数据 道间距 采样间隔
%% 
function [TrackInterval,dt,Ydata] = readB_scan(fileName)
% clear;
% clc ;
% close all;
FileHead=1024;
% [name, path]=uigetfile('D:/study/GPR_yi/sim/B_scan_imaging_interpreting/*.lte','打开雷达探测数据');
% fileName=strcat(path,name);
% if exist('fileName')==0
%     disp('File Not Exist');
% end  

fid =fopen(fileName);
% fseek(fid,4,'bof');
% fseek(fid,5,'bof');
fseek(fid,0,'bof');
tag=fread(fid,1,'int16');%读出标志位,2600为256,扫速、时窗等用整形，其余设备用float型
fseek(fid,4,'bof');
samPointNum=fread(fid,1,'int16');%读出采样点
fseek(fid,0,'bof');
p1=ftell(fid);
fseek(fid,0,'eof');
p2=ftell(fid);
lenperL=(p2-p1-FileHead)/(samPointNum*2);%从头文件中读出总的道数

fseek(fid,18,'bof');
PulseInterval = fread(fid,1,'float32');%读出脉冲间隔（第一个间隔）

fseek(fid,109,'bof');
Interval2 = fread(fid,1,'short');%读出第二个间隔

TrackInterval = PulseInterval*Interval2*1e-2;%道间距
if tag==256
    fseek(fid,26,'bof');
    timeWindow=fread(fid,1,'int32');%从头文件中读出时窗
    fseek(fid,10,'bof');
    scanspeed=fread(fid,1,'int32');%从头文件中读出扫速
else
    fseek(fid,26,'bof');
    timeWindow=fread(fid,1,'float32');%从头文件中读出时窗
    fseek(fid,10,'bof');
    scanspeed=fread(fid,1,'float32');%从头文件中读出扫速
end
dt = timeWindow/samPointNum*1e-9;         %采样间隔
fseek(fid,30,'bof');
freq=fread(fid,1,'int16');%从头文件中读出频率
fseek(fid,FileHead,'bof');
Ydata=fread(fid,[samPointNum,lenperL],'int16');%读出原始数据
fclose(fid);
end