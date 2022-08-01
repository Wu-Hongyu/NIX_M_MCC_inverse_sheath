%% Main function to generate tests
function tests = test_get_path_init
% test get_path_init
tests = functiontests(localfunctions);
end

%% Test Functions
function test_different_pwd(testCase)
% test different pwd
% �ȶ�ʵ�ʴ洢·����Ԥ�ڴ洢·��
project_name='test000';
path_nix_m = pwd;
excepted_path_list.root=path_nix_m;
excepted_path_list.packages=[path_nix_m '\packages'];
excepted_path_list.cross_sections=[path_nix_m '\packages\cross_section'];
excepted_path_list.others=[path_nix_m '\others'];
excepted_path_list.temp=[path_nix_m '\others\temp'];
excepted_path_list.project_name=project_name;
excepted_path_list.output=[path_nix_m '\others\' project_name];
% pwd=Nix_M
actual_path_list = get_path_init(project_name);
verifyEqual(testCase, actual_path_list, excepted_path_list)
% pwd=Nix_M/packages
cd('./packages')
actual_path_list = get_path_init(project_name);
verifyEqual(testCase, actual_path_list, excepted_path_list)
% pwd=Nix_M�ϼ�Ŀ¼
cd('../../')
e=get_exception('empty','');
try
get_path_init(project_name);
catch e
end
verifyEqual(testCase,e.message,'�޷�ȷ��Nix_M/·��')
cd(path_nix_m)
rmdir([path_nix_m '\others\' project_name])
end

function test_make_dir(testCase)
% test ����ָ���ļ���
% �ȶ�ʵ���ļ�������Ԥ���ļ�����
project_name='to_be_deleted';
path_nix_m = pwd;
get_path_init(project_name);
% �г�Nix_M���ļ�������
listing = dir(path_nix_m);
flag = [listing(:).isdir];
listing = {listing(flag).name}';
% listing(ismember(listing,{'.','..'})) = []; % Ӧ��Linux����.��..Ŀ¼
verifyEqual(testCase,sum(strcmp( listing , 'packages' )),1)
verifyEqual(testCase,sum(strcmp( listing , 'others' )),1)
% verifyEqual(testCase,sum(strcmp( listing , 'examples' )),1)
% �г�./other���ļ�������
listing = dir([path_nix_m '\others']);
flag = [listing(:).isdir];
listing = {listing(flag).name}';
verifyEqual(testCase,sum(strcmp( listing , 'temp' )),1)
verifyEqual(testCase,sum(strcmp( listing , project_name )),1)
rmdir([path_nix_m '\others\' project_name])
end

function test_make_dir2(testCase)
% test ʹ��ʱ��Ϊ����ļ�����
path_nix_m = pwd;
project_name=['project' datestr(now,'yyyy_mm_dd')];
flag_dir_already_exist=false;
if 7==exist([path_nix_m '\others\' project_name],'dir')
    flag_dir_already_exist=true;
end
get_path_init('now');
% �г�Nix_M���ļ�������
listing = dir([path_nix_m '\others']);
flag = [listing(:).isdir];
listing = {listing(flag).name}';
verifyEqual(testCase,sum(strcmp( listing , project_name )),1)
if ~flag_dir_already_exist
    warning('������0�����Ҳ��Կ��ܵ��²���ʧ�ܣ�������')
    rmdir([path_nix_m '\others\' project_name])
end
end

%% Optional file fixtures  
function setupOnce(testCase)  % do not change function name
% ����ǰ����Nix_MΪ��ǰ�ļ���
current_folder = pwd;
testCase.TestData.origPath =current_folder;
if strcmp(current_folder(end-5:end),'\Nix_M')
    path_root=current_folder;
else
    k = strfind(current_folder,'\Nix_M\');
    if 1==length(k)
        path_root=[current_folder(1:k) '\Nix_M'];
    else
        error('�޷�ȷ��Nix_M/·��')
    end
end
cd(path_root);
% ����ǰ��ɾ��test000��һ��������ʱ�ļ���
path_dir=[path_root '\others\test000'];
if 7==exist(path_dir,'dir')
    rmdir(path_dir,'s') %ɾ���ǿ��ļ���
end
end

function teardownOnce(testCase)  % do not change function name
% change back to original path, for example
cd(testCase.TestData.origPath)
end

%% Optional fresh fixtures  
function setup(testCase)  % do not change function name
% open a figure, for example
end

function teardown(testCase)  % do not change function name
% close figure, for example
end
