function [ path_list ] = get_path_init(project_name)
% ��ʼ��·�������س���·��
% ����Nix_M�ľ���·��
current_folder = pwd;
% ��ִ�б�������˵���Ѵ���Nix_M�������ļ�����
if strcmp(current_folder(end-5:end),'\Nix_M')
    path_list.root=current_folder;
else
    k = strfind(current_folder,'\Nix_M\');
    if 1==length(k)
        path_list.root=[current_folder(1:k) 'Nix_M'];
    else
        error('�޷�ȷ��Nix_M/·��')
    end
end
% ����ģ�������ļ���
path_list.packages=[path_list.root '\packages'];
addpath(genpath(path_list.packages)) % ���./packages�������ļ��е�·��
assert(7==exist(path_list.packages,'dir')) % ����ļ���
path_list.cross_sections=[path_list.packages '\cross_section'];
% % ʹ�ð��������ļ���
% path_list.examples=[path_list.root '\examples'];
addpath(genpath([path_list.root '\examples'])) % ���./packages�������ļ��е�·��
% assert(7==exist(path_list.packages,'dir')) % ����ļ���
% �����ļ������ļ���
path_list.others=[path_list.root '\others'];
make_dir(path_list.others)
path_list.temp=[path_list.others '\temp'];
make_dir(path_list.temp)
switch project_name
    case 'now' % ʹ��ʱ��Ϊ����ļ�����
        path_list.project_name=['project' datestr(now,'yyyy_mm_dd')];
%         path_list.project_name=datestr(now,30)
    otherwise % ʹ�������ַ���Ϊ����ļ�����
        path_list.project_name=project_name;
end
path_list.output=[path_list.others '\' path_list.project_name];
make_dir(path_list.output)
end

function make_dir(path)
% ����ļ��У����������½�
if 0==exist(path,'dir') % ����ļ������ļ���
    [status, msg]=mkdir(path);
    if ~status
        error(['[ERROR] mkdir ' path ' failed: ' msg])
    end
end
end