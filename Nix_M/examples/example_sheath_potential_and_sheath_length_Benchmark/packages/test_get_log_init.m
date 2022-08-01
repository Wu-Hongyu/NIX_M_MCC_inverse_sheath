%% Main function to generate tests
function tests = test_get_log_init
% test get_log_init
tests = functiontests(localfunctions);
end

%% Test Functions
function test_different_current_folder(testCase)
% test different_current_folder
% get_log_init.m���ܵĹ���Ŀ¼��Ҫ�Ǹ�Ŀ¼��./packages
% ������Ҫ����  ����������� �Ƿ�����./other������log�ļ�
simulation=get_simulation('default'); %������� �ṹ��
simulation.temp='������';
get_log_init( simulation.temp, 'g' );
current_folder = pwd;
if strcmp(current_folder(end-4:end),'Nix_M')
    % ����ǰ�ڸ�Ŀ¼��������./other�����ɵ�һ��log
    cd ./packages
    % ����./packages��׼����./other�����ɵڶ���log
    pause(65)
    % //log������ʱ��ֱ���Ϊ1min�������ʱ>1min
elseif strcmp(current_folder(end-7:end),'packages')
    % ����ǰ��./packages��������./other�����ɵ�һ��log
    cd ..
    % �����Ŀ¼��׼����./other�����ɵڶ���log
    pause(65)
end
get_log_init( simulation, 'g' );

answer = questdlg('Nix_M/others/��������gǰ׺�����־��','�˹��ж�','Y','N','Y');
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
