%% Main function to generate tests
function tests = test_get_field_init
% test get_field_init
tests = functiontests(localfunctions);
end

%% Test Functions
function test_case00(testCase)
% test ��ĩ��һ��߽�
simulation=get_simulation('default');
simulation.num_grid_point=4;
simulation.field_boundaries_type=0;
simulation.field_boundaries=[2,3];
[ u, A, b_extra ] = get_field_init( simulation );
% full(A)
verifyEqual(testCase,u(1),2)
verifyEqual(testCase,u(end),3)
constants=get_constants();% ȫ�ֳ��� �ṹ��
coeff_temp=-constants.eps0/(simulation.dx*simulation.dx);
% spdiags��װn �� n ��ϡ��ϵ������
% λ�����Խ����·��ĶԽ������ȴ��еĶ�����ȡԪ�ء�
% λ�����Խ����Ϸ��ĶԽ������ȴ��еĵײ���ȡԪ�ء�
verifyTrue(testCase,isequal(spdiags(A,0),coeff_temp*[1,-2,-2,1]'))
verifyTrue(testCase,isequal(spdiags(A,-1),coeff_temp*[1,1,0,0]'))
verifyTrue(testCase,isequal(spdiags(A,1),coeff_temp*[0,0,1,1]'))
verifyEqual(testCase,b_extra(1),coeff_temp*2)
verifyEqual(testCase,b_extra(end),coeff_temp*3)
end

function test_case01(testCase)
% test �׵�һ�࣬ĩ�ڶ���
simulation=get_simulation('default');
simulation.num_grid_point=4;
simulation.field_boundaries_type=1;
simulation.field_boundaries=[2,3];
[ u, A, b_extra ] = get_field_init( simulation );
% full(A)
verifyEqual(testCase,u(1),2)
verifyEqual(testCase,u(end),0)
constants=get_constants();% ȫ�ֳ��� �ṹ��
coeff_temp=-constants.eps0/(simulation.dx*simulation.dx);
verifyTrue(testCase,isequal(spdiags(A,0),coeff_temp*[1,-2,-2,-2]'))
verifyTrue(testCase,isequal(spdiags(A,-1),coeff_temp*[1,1,2,0]'))
verifyTrue(testCase,isequal(spdiags(A,1),coeff_temp*[0,0,1,1]'))
verifyEqual(testCase,b_extra(1),coeff_temp*2)
verifyEqual(testCase,b_extra(end),coeff_temp*(-2*3*simulation.dx))
end

function test_case10(testCase)
% test �׵ڶ��࣬ĩ��һ��
simulation=get_simulation('default');
simulation.num_grid_point=4;
simulation.field_boundaries_type=2;
simulation.field_boundaries=[2,3];
[ u, A, b_extra ] = get_field_init( simulation );
% full(A)
verifyEqual(testCase,u(1),0)
verifyEqual(testCase,u(end),3)
constants=get_constants();% ȫ�ֳ��� �ṹ��
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
