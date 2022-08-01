%% Main function to generate tests
function tests = test_get_field_init
% test get_field_init
tests = functiontests(localfunctions);
end

%% Test Functions
function test_case00(testCase)
% test 首末第一类边界
simulation=get_simulation('default');
simulation.num_grid_point=4;
simulation.field_boundaries_type=0;
simulation.field_boundaries=[2,3];
[ u, A, b_extra ] = get_field_init( simulation );
% full(A)
verifyEqual(testCase,u(1),2)
verifyEqual(testCase,u(end),3)
constants=get_constants();% 全局常数 结构体
coeff_temp=-constants.eps0/(simulation.dx*simulation.dx);
% spdiags组装n × n 的稀疏系数矩阵
% 位于主对角线下方的对角线首先从列的顶部获取元素。
% 位于主对角线上方的对角线首先从列的底部获取元素。
verifyTrue(testCase,isequal(spdiags(A,0),coeff_temp*[1,-2,-2,1]'))
verifyTrue(testCase,isequal(spdiags(A,-1),coeff_temp*[1,1,0,0]'))
verifyTrue(testCase,isequal(spdiags(A,1),coeff_temp*[0,0,1,1]'))
verifyEqual(testCase,b_extra(1),coeff_temp*2)
verifyEqual(testCase,b_extra(end),coeff_temp*3)
end

function test_case01(testCase)
% test 首第一类，末第二类
simulation=get_simulation('default');
simulation.num_grid_point=4;
simulation.field_boundaries_type=1;
simulation.field_boundaries=[2,3];
[ u, A, b_extra ] = get_field_init( simulation );
% full(A)
verifyEqual(testCase,u(1),2)
verifyEqual(testCase,u(end),0)
constants=get_constants();% 全局常数 结构体
coeff_temp=-constants.eps0/(simulation.dx*simulation.dx);
verifyTrue(testCase,isequal(spdiags(A,0),coeff_temp*[1,-2,-2,-2]'))
verifyTrue(testCase,isequal(spdiags(A,-1),coeff_temp*[1,1,2,0]'))
verifyTrue(testCase,isequal(spdiags(A,1),coeff_temp*[0,0,1,1]'))
verifyEqual(testCase,b_extra(1),coeff_temp*2)
verifyEqual(testCase,b_extra(end),coeff_temp*(-2*3*simulation.dx))
end

function test_case10(testCase)
% test 首第二类，末第一类
simulation=get_simulation('default');
simulation.num_grid_point=4;
simulation.field_boundaries_type=2;
simulation.field_boundaries=[2,3];
[ u, A, b_extra ] = get_field_init( simulation );
% full(A)
verifyEqual(testCase,u(1),0)
verifyEqual(testCase,u(end),3)
constants=get_constants();% 全局常数 结构体
coeff_temp=-constants.eps0/(simulation.dx*simulation.dx);
verifyTrue(testCase,isequal(spdiags(A,0),coeff_temp*[-2,-2,-2,1]'))
verifyTrue(testCase,isequal(spdiags(A,-1),coeff_temp*[1,1,0,0]'))
verifyTrue(testCase,isequal(spdiags(A,1),coeff_temp*[0,2,1,1]'))
verifyEqual(testCase,b_extra(1),coeff_temp*(2*2*simulation.dx))
verifyEqual(testCase,b_extra(end),coeff_temp*3)
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
