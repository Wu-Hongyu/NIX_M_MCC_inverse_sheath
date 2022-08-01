%% Main function to generate tests
function tests = test_log
% test get_log_init
tests = functiontests(localfunctions);
end

%% Test Functions
function test_get_log_init(testCase)
% test get_log_init
project_name='test000';
prefix_str='test';
path_list=get_path_init(project_name);
path_nix_m=path_list.root;
log_name=get_log_init( path_list.output, prefix_str );
now_str=datestr(now,'yyyy_mm_dd_HH');
warning('在跨小时测试可能导致测试失败，请无视')
A = dir([path_nix_m '\others\' project_name '\' prefix_str now_str '*.log']);
verifyFalse(testCase, isempty(A)) %A非空，存在该文件
verifyEqual(testCase, A.bytes, 21) %检查写入内容
verifyEqual(testCase, [A.folder '\' A.name],log_name) %检查返回值
rmdir([path_nix_m '\others\' project_name],'s')
end

%% Optional file fixtures  
function setupOnce(testCase)  % do not change function name
% set a new path, for example
path_list=get_path_init('test000');
% 测试前先删除test000这一测试用临时文件夹
path_dir=[path_list.root '\others\test000'];
if 7==exist(path_dir,'dir')
    rmdir(path_dir,'s') %删除非空文件夹
end
end

function teardownOnce(testCase)  % do not change function name
% change back to original path, for example
end

%% Optional fresh fixtures  
function setup(testCase)  % do not change function name
% open a figure, for example
end

function teardown(testCase)  % do not change function name
% close figure, for example
end
