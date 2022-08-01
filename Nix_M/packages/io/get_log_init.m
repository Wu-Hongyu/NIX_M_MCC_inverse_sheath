function [ log_name ] = get_log_init( path, prefix_str )
% 生成初始日志，返回日志文件名
% 在path中生成log文件，以前缀和时间自动命名
now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
log_name=[prefix_str datestr(now,'yyyy_mm_dd_HH_MM') '.log'];
% 1min内连续初始化日志，会导致连续写入同一log文件
log_name=[path '\' log_name];

% 重定向方式log到文件
diary(log_name)
disp(now_str)
diary off
end

