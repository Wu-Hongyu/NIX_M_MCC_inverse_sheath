%% Main function to generate tests
function tests = test_get_path_init
% test get_path_init
tests = functiontests(localfunctions);
end

%% Test Functions
function test_different_pwd(testCase)
% test different pwd
% 比对实际存储路径与预期存储路径
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
% pwd=Nix_M上级目录
cd('../../')
e=get_exception('empty','');
try
get_path_init(project_name);
catch e
end
verifyEqual(testCase,e.message,'无法确定Nix_M/路径')
cd(path_nix_m)
rmdir([path_nix_m '\others\' project_name])
end

function test_make_dir(testCase)
% test 生成指定文件夹
% 比对实际文件夹名与预期文件夹名
project_name='to_be_deleted';
path_nix_m = pwd;
get_path_init(project_name);
% 列出Nix_M子文件夹名称
listing = dir(path_nix_m);
flag = [listing(:).isdir];
listing = {listing(flag).name}';
% listing(ismember(listing,{'.','..'})) = []; % 应对Linux风格的.与..目录
verifyEqual(testCase,sum(strcmp( listing , 'packages' )),1)
verifyEqual(testCase,sum(strcmp( listing , 'others' )),1)
% verifyEqual(testCase,sum(strcmp( listing , 'examples' )),1)
% 列出./other子文件夹名称
listing = dir([path_nix_m '\others']);
flag = [listing(:).isdir];
listing = {listing(flag).name}';
verifyEqual(testCase,sum(strcmp( listing , 'temp' )),1)
verifyEqual(testCase,sum(strcmp( listing , project_name )),1)
rmdir([path_nix_m '\others\' project_name])
end

function test_make_dir2(testCase)
% test 使用时间为输出文件夹名
path_nix_m = pwd;
project_name=['project' datestr(now,'yyyy_mm_dd')];
flag_dir_already_exist=false;
if 7==exist([path_nix_m '\others\' project_name],'dir')
    flag_dir_already_exist=true;
end
get_path_init('now');
% 列出Nix_M子文件夹名称
listing = dir([path_nix_m '\others']);
flag = [listing(:).isdir];
listing = {listing(flag).name}';
verifyEqual(testCase,sum(strcmp( listing , project_name )),1)
if ~flag_dir_already_exist
    warning('在晚上0点左右测试可能导致测试失败，请无视')
    rmdir([path_nix_m '\others\' project_name])
end
end

%% Optional file fixtures  
function setupOnce(testCase)  % do not change function name
% 测试前先以Nix_M为当前文件夹
current_folder = pwd;
testCase.TestData.origPath =current_folder;
if strcmp(current_folder(end-5:end),'\Nix_M')
    path_root=current_folder;
else
    k = strfind(current_folder,'\Nix_M\');
    if 1==length(k)
        path_root=[current_folder(1:k) '\Nix_M'];
    else
        error('无法确定Nix_M/路径')
    end
end
cd(path_root);
% 测试前先删除test000这一测试用临时文件夹
path_dir=[path_root '\others\test000'];
if 7==exist(path_dir,'dir')
    rmdir(path_dir,'s') %删除非空文件夹
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
