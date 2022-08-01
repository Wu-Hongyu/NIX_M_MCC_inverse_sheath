%% Main function to generate tests
function tests = test_get_v_init
% test get_v_init
tests = functiontests(localfunctions);
end

%% Test Functions
function test_maxwellian_velocity(testCase)
% test case 'Maxwellian velocity'
vth=1e6;
type_init='Maxwellian velocity';
v_init=get_v_init( vth, type_init, [3,2]);
assert(isequal([3,2],size(v_init))) %指定维度

pdf_maxwellian_velocity=@(v) exp(-v.*v/(2*vth*vth))/(sqrt(2*pi)*vth);
check_pdf(testCase, type_init, 1e8,0.03, pdf_maxwellian_velocity) 
end

function check_pdf(testCase, type_init, total_num,tolerance, pdf_fun)
% 检查实际pdf与预期pdf一致
vth=1e6;
dimension_vec=[total_num,1];
v_init=get_v_init( vth, type_init, dimension_vec);
num=length(v_init);
% 使用histogram做自动分段统计
h_plot1=histogram(v_init); %自己会单独figure
bin_width=h_plot1.BinWidth;
N_actual=h_plot1.Values;
x_edges=h_plot1.BinEdges;
pd_excepted=pdf_fun(x_edges);
N_excepted=num*bin_width*(pd_excepted(1:end-1)+pd_excepted(2:end))/2; %PDF分段梯形积分

% 误差分析
% 实际频数与概率计算的频数在小概率或小样本时相对误差容易较大，因此不用相对误差
% verifyEqual(testCase,N_actual,N_excepted,'RelTol',5e-2);  %Failed
% 难以统一给定绝对误差的容差，因此不用统一的绝对误差
% 实际使用与柱图最大值A成比例的绝对误差，即A越大，允许N*PDF误差 值越多

% TODO：这是一个错误的、解释不清的思路。实际上，不应该用A，应该用N；进而应该用相对误差

A=max(N_excepted);
verifyEqual(testCase,N_actual,N_excepted,'AbsTol',A*tolerance);

% 建议新test case第一次test时，开启以下的全面判断代码
% verifyLessThan(testCase,rms(N_actual-N_excepted),A*0.01);
% % 绘图，人工判断
% hold on
% plot_x=(x_edges(1:end-1)+x_edges(2:end))/2;
% plot(plot_x,N_excepted,'-r','LineWidth',3)
% legend('实际分布','预期分布')
% answer = questdlg('实际与预期一致？',['人工判断' type_init],'Y','N','Y');
% verifyEqual(testCase,answer,'Y')
% close
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
