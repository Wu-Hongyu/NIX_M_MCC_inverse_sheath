function [ path_list ] = get_path_init(project_name)
% 初始化路径，返回常用路径
% 返回Nix_M的绝对路径
current_folder = pwd;
% 能执行本函数，说明已处于Nix_M或其子文件夹内
if strcmp(current_folder(end-5:end),'\Nix_M')
    path_list.root=current_folder;
else
    k = strfind(current_folder,'\Nix_M\');
    if 1==length(k)
        path_list.root=[current_folder(1:k) 'Nix_M'];
    else
        error('无法确定Nix_M/路径')
    end
end
% 功能模块所在文件夹
path_list.packages=[path_list.root '\packages'];
addpath(genpath(path_list.packages)) % 添加./packages及其子文件夹到路径
assert(7==exist(path_list.packages,'dir')) % 检查文件夹
path_list.cross_sections=[path_list.packages '\cross_section'];
% % 使用案例所在文件夹
% path_list.examples=[path_list.root '\examples'];
addpath(genpath([path_list.root '\examples'])) % 添加./packages及其子文件夹到路径
% assert(7==exist(path_list.packages,'dir')) % 检查文件夹
% 其他文件所在文件夹
path_list.others=[path_list.root '\others'];
make_dir(path_list.others)
path_list.temp=[path_list.others '\temp'];
make_dir(path_list.temp)
switch project_name
    case 'now' % 使用时间为输出文件夹名
        path_list.project_name=['project' datestr(now,'yyyy_mm_dd')];
%         path_list.project_name=datestr(now,30)
    otherwise % 使用输入字符串为输出文件夹名
        path_list.project_name=project_name;
end
path_list.output=[path_list.others '\' path_list.project_name];
make_dir(path_list.output)
end

function make_dir(path)
% 检查文件夹，不存在则新建
if 0==exist(path,'dir') % 输出文件所在文件夹
    [status, msg]=mkdir(path);
    if ~status
        error(['[ERROR] mkdir ' path ' failed: ' msg])
    end
end
end