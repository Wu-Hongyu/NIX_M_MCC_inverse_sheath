%% Main function to generate tests
function tests = test_get_log_init
% test get_log_init
tests = functiontests(localfunctions);
end

%% Test Functions
function test_different_current_folder(testCase)
% test different_current_folder
% get_log_init.m可能的工作目录主要是根目录或./packages
% 以下主要测试  这两种情况下 是否能在./other中生成log文件
simulation=get_simulation('default'); %仿真参数 结构体
simulation.temp='测试用';
get_log_init( simulation.temp, 'g' );
current_folder = pwd;
if strcmp(current_folder(end-4:end),'Nix_M')
    % 若当前在根目录，则已在./other中生成第一个log
    cd ./packages
    % 进入./packages，准备在./other中生成第二个log
    pause(65)
    % //log命名的时间分辨率为1min，因此延时>1min
elseif strcmp(current_folder(end-7:end),'packages')
    % 若当前在./packages，则已在./other中生成第一个log
    cd ..
    % 进入根目录，准备在./other中生成第二个log
    pause(65)
end
get_log_init( simulation, 'g' );

answer = questdlg('Nix_M/others/下有两个g前缀输出日志？','人工判断','Y','N','Y');
verifyEqual(testCase,answer,'Y')
end

%% Optional file fixtures  
function setupOnce(testCase)  % do not change function name
% set a new path, for example
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
