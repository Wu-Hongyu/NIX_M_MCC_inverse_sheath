function [ log_name ] = get_log_init( in_var, prefix_str )
% 生成初始日志，返回日志文件名
% 在不同当前工作目录下，可以在./other中生成log文件，以前缀和时间自动命名
now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
log_name=[prefix_str datestr(now,'yyyy_mm_dd_HH_MM') '.log'];
% 判断当前文件夹，尽可能在./other中生成log文件
current_folder = pwd; 
if strcmp(current_folder(end-4:end),'Nix_M')
    log_name=['./others/' log_name];
elseif strcmp(current_folder(end-7:end),'packages')
    log_name=['../others/' log_name]; % ..为上级目录（即Nix_M）路径
end

% 重定向方式
diary(log_name)
disp(now_str)
disp(in_var)
diary off

% % 写文件方式
% log_str_cell{1}=[now_str '\r\n\r\n'];
% log_str_cell{2}=['dt=' num2str(simulation.dt) ';  时间步长\r\n'];
% log_str_cell{3}=['EndTime=' num2str(simulation.end_time) ';  仿真时间总长\r\n'];
% log_str_cell{4}=['AllTimeSteps=' num2str(simulation.all_timesteps) ';  总循环次数\r\n'];
% log_str_cell{5}=['Lx=' num2str(simulation.Lx) ';  仿真区域长度\r\n'];
% log_str_cell{6}=['NG=' num2str(simulation.num_grid_point) ';  网格格点数目\r\n'];
% log_str_cell{7}=['dx=' num2str(simulation.dx) ';  空间步长\r\n'];
% log_str_cell{8}=['n0=' num2str(simulation.n0) ';  等离子体密度\r\n'];
% log_str_cell{9}=['PeN=' num2str(simulation.num0_macro_e) ';  初始电子数目 Particle e Num\r\n'];
% log_str_cell{10}=['PHN=' num2str(simulation.num0_macro_H) ';  初始离子数目Particle H Num\r\n'];
% log_str_cell{11}=['weight=' num2str(simulation.weight) ';  每个宏粒子代表实际粒子的权重\r\n'];
% log_str_cell{12}=['Te=' num2str(simulation.Te) ';  1eV的温度，单位K\r\n'];
% log_str_cell{13}=['TH=' num2str(simulation.TH) ';  1eV的温度，单位K\r\n'];
% fp=fopen(log_name,'a');
% for logi=1:13
    % fprintf(log_str_cell{logi})
    % fprintf(fp,log_str_cell{logi});
% end
% fclose(fp);
end

