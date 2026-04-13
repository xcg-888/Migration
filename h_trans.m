function epsilon_r = h_trans(data_out, dt)
%% --- 基于已知厚度反演相对介电常数 (双模式合并版) ---
% 输入:
%   data_out: 数据矩阵
%   dt: 采样间隔 (程序会自动识别秒或纳秒)

% ================= 用户参数设置 (请在此处修改) =================
known_depth_m = 0.255
% --- 选择校正模式 ---
% 1 = 【自动模式】: 程序自动找零点，你只需点 1 次目标 (快捷)
% 2 = 【手动模式】: 你手动点 2 次 (第一次点表面，第二次点目标) -> 最稳妥，推荐！
correction_mode = 2;       
% ==========================================================

% --- 1. 单位自动检测与转换 (防止计算结果为0) ---
if dt < 1e-4
    dt = dt * 1e9; % 如果单位是秒，转为 ns
    fprintf('提示: 已自动将 dt 单位转换为纳秒 (%.4f ns)\n', dt);
end

% 构建基础时间轴
[num_samples, num_traces] = size(data_out);
time_window_ns = dt * num_samples;
base_time_axis = linspace(0, time_window_ns, num_samples);

% 变量初始化
final_time_axis = base_time_axis;
t_zero_val = 0;

% --- 2. 根据模式处理时间轴 ---
if correction_mode == 1
    % === 自动模式 ===
    disp('【当前模式：自动零时矫正】');
    % 计算平均道寻找直达波
    mean_trace = mean(abs(data_out), 2);
    search_limit = floor(num_samples * 0.2); 
    if search_limit < 5, search_limit = num_samples; end
    [~, zero_idx] = max(mean_trace(1:search_limit));
    
    t_zero_val = base_time_axis(zero_idx);
    
    % 平移时间轴 (使得 0ns 对应地表)
    final_time_axis = base_time_axis - t_zero_val;
    fprintf('程序判定地表直达波位置: %.4f ns\n', t_zero_val);
    
else
    % === 手动模式 ===
    disp('【当前模式：手动双点拾取】');
    disp('完全保留原始时间轴，不做自动平移。');
    final_time_axis = base_time_axis;
end

% --- 3. 绘制图像 ---
figure;
imagesc(1:num_traces, final_time_axis, data_out);
colormap('gray'); 

% 根据模式设置标题和参考线
if correction_mode == 1
    ylabel('校正后走时 (ns)');
    title(['自动模式: 请点击 1 次【目标反射】 (地表已对齐0ns)']);
    hold on; yline(0, 'r-', 'LineWidth', 1.5, 'Label', 'Surface (Auto)'); hold off;
else
    ylabel('原始走时 (ns)');
    title(['厚度反演']);
end

% --- 4. 交互操作 ---
disp('-------------------------------------------------------');
if correction_mode == 1
    disp('请在图中点击：');
    disp(' >> 底层反射波)');
    [~, y_clicks] = ginput(1);
    
    t_surface_picked = 0;          % 自动模式下，表面就是0
    t_target_picked = y_clicks(1); % 读取点击值
    t_travel = t_target_picked;    % 走时即为点击值
    
else
    disp('请按顺序在图中点击：');
    disp(' 1. 介质表面 / 空气直达波 (通常是最上方的强黑/白线)');
    disp(' 2. 底层反射波');
    [~, y_clicks] = ginput(2);
    
    t_surface_picked = y_clicks(1);
    t_target_picked = y_clicks(2);
    
    % 计算差值
    t_travel = t_target_picked - t_surface_picked;
end
disp('-------------------------------------------------------');

% --- 5. 安全检查 ---
if t_travel <= 0
    error('错误: 计算出的走时 <= 0。请检查点击顺序 (手动模式需先点上方的地表，再点下方的目标)。');
end

% --- 6. 计算 ---
% 公式: er = ((c * t) / (2 * h))^2
c_vacuum = 0.3; 
epsilon_r = ((c_vacuum * t_travel) / (2 * known_depth_m))^2;
v_medium = c_vacuum / sqrt(epsilon_r);

% --- 7. 结果输出 ---
fprintf('\n=== 计算结果 ===\n');
fprintf('已知深度 (h): %.4f m\n', known_depth_m);
if correction_mode == 2
    fprintf('手动拾取地表时刻: %.4f ns\n', t_surface_picked);
    fprintf('手动拾取目标时刻: %.4f ns\n', t_target_picked);
end
fprintf('有效双程走时 (t): %.4f ns\n', t_travel);
fprintf('--------------------------\n');
fprintf('计算波速 (v): %.4f m/ns\n', v_medium);
fprintf('>> 相对介电常数 (εr): %.4f <<\n', epsilon_r);
fprintf('========================\n');

% --- 8. 标记结果 ---
hold on;
if correction_mode == 1
    yline(t_target_picked, 'g--', 'LineWidth', 1.5, 'Label', 'Target');
else
    yline(t_surface_picked, 'r--', 'LineWidth', 1.5, 'Label', 'Surface');
    yline(t_target_picked, 'g--', 'LineWidth', 1.5, 'Label', 'Target');
end
plot(xlim, [t_target_picked, t_target_picked], 'g*');
hold off;

end