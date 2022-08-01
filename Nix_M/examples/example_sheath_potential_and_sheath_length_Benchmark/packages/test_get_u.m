%% Main function to generate tests
function tests = test_get_u
% test get_u
tests = functiontests(localfunctions);
end

%% Test Functions
% 测试系数矩阵方程求解器
function test_dircet_inverse(testCase)
% test case 'direct inverse'
rng default
n=400;
u0=zeros(n,1);
A = sprand(n,n,0.5); %创建密度为 50% 的随机稀疏矩阵
A = A'*A;
b = sum(A,2); %由 A 的行总和组成的右端向量 b
% Au=b的解u为全1的列向量
type_solver='direct inverse';
u=get_u( u0,A,b,type_solver );
flag=sum(abs(u-ones(n,1))>1e-7);
verifyEqual(testCase,flag,0)
end

function test_cg(testCase)
% test case 'CG'
rng default
n=400;
u0=zeros(n,1);
A = sprand(n,n,0.5); %创建密度为 50% 的随机稀疏矩阵
A = A'*A;
b = sum(A,2); %由 A 的行总和组成的右端向量 b
type_solver='CG';
u=get_u( u0,A,b,type_solver );
flag=sum(abs(u-ones(n,1))>5e-2);
%误差容差5%，非相对残差的容差
verifyEqual(testCase,flag,0)
end

function test_iccg(testCase)
% test case 'ICCG'
rng default
n=400;
u0=zeros(n,1);
A = sprand(n,n,0.5); %创建密度为 50% 的随机稀疏矩阵
A = A'*A;
b = sum(A,2); %由 A 的行总和组成的右端向量 b
type_solver='ICCG';
u=get_u( u0,A,b,type_solver );
flag=sum(abs(u-ones(n,1))>5e-2);
verifyEqual(testCase,flag,0)
end

function test_bicgstab(testCase)
% test case 'BiCGstab'
rng default
n=400;
u0=zeros(n,1);
A = sprand(n,n,0.5); %创建密度为 50% 的随机稀疏矩阵
A = A'*A;
b = sum(A,2); %由 A 的行总和组成的右端向量 b
type_solver='BiCGstab';
u=get_u( u0,A,b,type_solver );
flag=sum(abs(u-ones(n,1))>5e-2);
verifyEqual(testCase,flag,0)
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
