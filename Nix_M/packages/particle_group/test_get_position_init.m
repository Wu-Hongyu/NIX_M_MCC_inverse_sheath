%% Main function to generate tests
function tests = test_get_position_init
% test get_position_init
tests = functiontests(localfunctions);
end

%% Test Functions
function test_e_d_u_n(testCase)
% test case 'entire domain uniform+noise'
simulation=get_simulation('default'); %仿真参数 结构体
type_init='entire domain uniform+noise';
check_uniform_distribution(testCase,simulation,type_init,1e8,0.05,simulation.Lx)
end

function test_e_d_u_r(testCase)
% test case 'entire domain uniform random'
simulation=get_simulation('default'); %仿真参数 结构体
type_init='entire domain uniform random';
check_uniform_distribution(testCase,simulation,type_init,1e8,0.05,simulation.Lx)
end

function test_s_r_u_r(testCase)
% test case 'source region uniform random' 
simulation=get_simulation('default'); %仿真参数 结构体
type_init='source region uniform random' ;
check_uniform_distribution(testCase,simulation,type_init,1e8,0.05,simulation.source_region(2)-simulation.source_region(1))
end

%% Aid function
function check_uniform_distribution(testCase,simulation,type_init,total_num,tolerance,interval)
dimension_vec=[total_num,1];
position_init=get_position_init( simulation, type_init, dimension_vec );
assert(isequal(dimension_vec,size(position_init)))

num=length(position_init);
h_plot1=histogram(position_init);
N_actual=h_plot1.Values;
N_excepted=num*h_plot1.BinWidth*1/interval*ones(1,h_plot1.NumBins);
% 误差分析
% 实际频数与概率计算的频数在小概率或小样本时相对误差容易较大
% verifyEqual(testCase,N_actual,N_excepted,'RelTol',5e-2);  %Failed
% 难以统一给定绝对误差的容差
A=max(N_excepted);
% verifyEqual(testCase,N_actual,N_excepted,'AbsTol',A*0.05); %Failed
% 分段最后一段可能频数很小
verifyLessThan(testCase,rms(N_actual-N_excepted),A*tolerance);
% % 绘图，人工判断
% hold on
% x_edges=h_plot1.BinEdges;
% x_min=x_edges(1);
% x_max=x_edges(end);
% line([x_min,x_max],[A,A],'Color','r','LineWidth',3);
% axis([x_min,x_max,-inf,inf])
% legend('实际分布','预期分布')
% answer = questdlg('实际与预期一致？',['人工判断' type_init],'Y','N','Y');
% verifyEqual(testCase,answer,'Y')
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
figure
end

function teardown(testCase)  % do not change function name
% close figure, for example
close
end
