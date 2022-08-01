%% Main function to generate tests
function tests = test_get_u
% test get_u
tests = functiontests(localfunctions);
end

%% Test Functions
% 测试系数矩阵方程求解器
function test_dircet_inverse(testCase)
% test case 'direct inverse'
solve_simple_equation1(testCase,'direct inverse', 1e-7)
end

function test_cg(testCase)
% test case 'CG'
solve_simple_equation1(testCase,'CG', 5e-2)
end

function test_iccg(testCase)
% test case 'ICCG'
solve_simple_equation1(testCase,'ICCG', 5e-2)
end

function test_bicgstab(testCase)
% test case 'BiCGstab'
solve_simple_equation1(testCase,'BiCGstab', 5e-2)
end

%% aid function
function solve_simple_equation1(testCase,type_solver, relative_tolerance)
rng default
n=400;
u0=zeros(n,1);
A = sprand(n,n,0.5); %创建密度为 50% 的随机稀疏矩阵
A = A'*A;
b = sum(A,2); %由 A 的行总和组成的右端向量 b
u=get_u( u0,A,b,type_solver ); %解应为全1
u=full(u);
% 误差分析
% 注意，迭代法中残差与求解结果的误差，二者概念不同，因此容差也不同
verifyEqual(testCase,u,ones(n,1),'RelTol',relative_tolerance); 
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
