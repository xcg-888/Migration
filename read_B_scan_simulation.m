%% 读取仿真数据
function [field, time, traces,dt] = read_B_scan_simulation()
[filename, pathname] = uigetfile('*.out', 'Select gprMax output file to plot B-scan', 'MultiSelect', 'on');
fullfilename = fullfile(pathname, filename);% 弹出文件选择框并进行选择

if fullfilename ~= 0
    % 读取HDF5文件属性
    header.iterations = double(h5readatt(fullfilename, '/', 'Iterations')); % 仿真运行的总时间步数
    header.dt = h5readatt(fullfilename, '/', 'dt'); % 仿真的时间步长
    header.srcsteps = 0.01;
    % header.srcsteps = h5readatt(fullfilename, '/rxs', 'step')

    % 读取电场数据（Ez分量）
    fieldpath = strcat('/rxs/rx1/', 'Ez'); % 连接字符串，数据集路径
    field = h5read(fullfilename, fieldpath)'; % 转置为时间步×扫描道格式

    % 生成时间轴和扫描道索引
    time = linspace(0, (header.iterations - 1) * header.dt, header.iterations)'; % linspace(start,end,n)
    traces = 1:size(field, 2);
    dt =  header.dt;
%     X = traces*header.srcsteps;
end
end