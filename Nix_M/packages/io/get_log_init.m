function [ log_name ] = get_log_init( path, prefix_str )
% ���ɳ�ʼ��־��������־�ļ���
% ��path������log�ļ�����ǰ׺��ʱ���Զ�����
now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
log_name=[prefix_str datestr(now,'yyyy_mm_dd_HH_MM') '.log'];
% 1min��������ʼ����־���ᵼ������д��ͬһlog�ļ�
log_name=[path '\' log_name];

% �ض���ʽlog���ļ�
diary(log_name)
disp(now_str)
diary off
end

