%% Main function to generate tests
function tests = test_exception
% test 异常处理package及语法
tests = functiontests(localfunctions);
end

%% Test Functions
function test_get_exception(testCase)
% test get_exception
% type: empty
actual_e=get_exception('empty','');
excepted_e=MException('Nix_M:EMPTY','[EMPTY] ');
verifyEqual(testCase,actual_e,excepted_e) %均未抛出，无stack信息

% type: other
test_message='test get_exception';
try
throw(get_exception('error',test_message)) %有stack信息
catch actual_e
end
excepted_e=MException('Nix_M:ERROR',['[ERROR] ' test_message]); %无stack信息
verify_exception_equal(testCase,actual_e,excepted_e) %只比较type和message
end

function test_output_exception_type0(testCase)
% test output_exception output_type=0
e=get_exception('empty',''); %扩展异常变量作用范围
test_message='test output_exception';
try
throw(get_exception('error',test_message)) %有stack信息
catch e %符合预期地产生异常变量
end
% 重定向方式log到文件
path_list=get_path_init('test000');
actual_log_name=[path_list.output '\test_output_exception_type0.log'];
diary(actual_log_name)
output_exception(0, e); 
diary off
excepted_log_name=[path_list.packages '\exception\test_output_exception_type0.log'];
actual_log=read_file( actual_log_name, 'rt2char' );
excepted_log=read_file( excepted_log_name, 'rt2char' );
verifyEqual(testCase,actual_log(1:500),excepted_log(1:500))  %运行位置不一样，最后几行会不一样，导致报错
warning('当修改本文件致throw行号改变时，测试报错，需人工对比文件内容是否大致一致')
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
