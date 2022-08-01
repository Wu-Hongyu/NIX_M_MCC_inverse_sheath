%% Main function to generate tests
function tests = test_io_file
% test io file
tests = functiontests(localfunctions);
end

%% Test Functions
function test_io_char_t(testCase)
% test case 'rt2char' and 'char2wt'
path_list=get_path_init('test000');
filename=[path_list.output '\test_io_char_t.txt'];
output='test\r\nio char t';
write_file( filename, 'char2wt', output )
input=read_file( filename, 'rt2char' );
verifyEqual(testCase,input,output)
delete(filename)
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