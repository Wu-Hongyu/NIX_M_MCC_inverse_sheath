function [ log_name ] = get_log_init( in_var, prefix_str )
% ���ɳ�ʼ��־��������־�ļ���
% �ڲ�ͬ��ǰ����Ŀ¼�£�������./other������log�ļ�����ǰ׺��ʱ���Զ�����
now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
log_name=[prefix_str datestr(now,'yyyy_mm_dd_HH_MM') '.log'];
% �жϵ�ǰ�ļ��У���������./other������log�ļ�
current_folder = pwd; 
if strcmp(current_folder(end-4:end),'Nix_M')
    log_name=['./others/' log_name];
elseif strcmp(current_folder(end-7:end),'packages')
    log_name=['../others/' log_name]; % ..Ϊ�ϼ�Ŀ¼����Nix_M��·��
end

% �ض���ʽ
diary(log_name)
disp(now_str)
disp(in_var)
diary off

% % д�ļ���ʽ
% log_str_cell{1}=[now_str '\r\n\r\n'];
% log_str_cell{2}=['dt=' num2str(simulation.dt) ';  ʱ�䲽��\r\n'];
% log_str_cell{3}=['EndTime=' num2str(simulation.end_time) ';  ����ʱ���ܳ�\r\n'];
% log_str_cell{4}=['AllTimeSteps=' num2str(simulation.all_timesteps) ';  ��ѭ������\r\n'];
% log_str_cell{5}=['Lx=' num2str(simulation.Lx) ';  �������򳤶�\r\n'];
% log_str_cell{6}=['NG=' num2str(simulation.num_grid_point) ';  ��������Ŀ\r\n'];
% log_str_cell{7}=['dx=' num2str(simulation.dx) ';  �ռ䲽��\r\n'];
% log_str_cell{8}=['n0=' num2str(simulation.n0) ';  ���������ܶ�\r\n'];
% log_str_cell{9}=['PeN=' num2str(simulation.num0_macro_e) ';  ��ʼ������Ŀ Particle e Num\r\n'];
% log_str_cell{10}=['PHN=' num2str(simulation.num0_macro_H) ';  ��ʼ������ĿParticle H Num\r\n'];
% log_str_cell{11}=['weight=' num2str(simulation.weight) ';  ÿ�������Ӵ���ʵ�����ӵ�Ȩ��\r\n'];
% log_str_cell{12}=['Te=' num2str(simulation.Te) ';  1eV���¶ȣ���λK\r\n'];
% log_str_cell{13}=['TH=' num2str(simulation.TH) ';  1eV���¶ȣ���λK\r\n'];
% fp=fopen(log_name,'a');
% for logi=1:13
    % fprintf(log_str_cell{logi})
    % fprintf(fp,log_str_cell{logi});
% end
% fclose(fp);
end

