%% Main function to generate tests
function tests = test_exception
% test �쳣����package���﷨
tests = functiontests(localfunctions);
end

%% Test Functions
function test_get_exception(testCase)
% test get_exception
% type: empty
actual_e=get_exception('empty','');
excepted_e=MException('Nix_M:EMPTY','[EMPTY] ');
verifyEqual(testCase,actual_e,excepted_e) %��δ�׳�����stack��Ϣ

% type: other
test_message='test get_exception';
try
throw(get_exception('error',test_message)) %��stack��Ϣ
catch actual_e
end
excepted_e=MException('Nix_M:ERROR',['[ERROR] ' test_message]); %��stack��Ϣ
verify_exception_equal(testCase,actual_e,excepted_e) %ֻ�Ƚ�type��message
end

function test_output_exception_type0(testCase)
% test output_exception output_type=0
e=get_exception('empty',''); %��չ�쳣�������÷�Χ
test_message='test output_exception';
try
throw(get_exception('error',test_message)) %��stack��Ϣ
catch e %����Ԥ�ڵز����쳣����
end
% �ض���ʽlog���ļ�
path_list=get_path_init('test000');
actual_log_name=[path_list.output '\test_output_exception_type0.log'];
diary(actual_log_name)
output_exception(0, e); 
diary off
excepted_log_name=[path_list.packages '\exception\test_output_exception_type0.log'];
actual_log=read_file( actual_log_name, 'rt2char' );
excepted_log=read_file( excepted_log_name, 'rt2char' );
verifyEqual(testCase,actual_log(1:500),excepted_log(1:500))  %����λ�ò�һ��������л᲻һ�������±���
warning('���޸ı��ļ���throw�кŸı�ʱ�����Ա������˹��Ա��ļ������Ƿ����һ��')
delete(actual_log_name)
end

%% Aid function
function verify_exception_equal(testCase,actual_e,excepted_e)
verifyEqual(testCase,actual_e.identifier,excepted_e.identifier)
verifyEqual(testCase,actual_e.message,excepted_e.message)
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
