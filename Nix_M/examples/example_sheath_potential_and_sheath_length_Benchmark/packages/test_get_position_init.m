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
dimension_vec=[1e5,1];
position_init=get_position_init( simulation, type_init, dimension_vec );
assert(isequal(dimension_vec,size(position_init)))

num=length(position_init);
h_plot1=histogram(position_init);
bin_width=h_plot1.BinWidth;
hold on
n=num*bin_width*1/simulation.Lx;
line([0,simulation.Lx],[n,n],'Color','r','LineWidth',3);
axis([0,simulation.Lx,-inf,inf])
legend('实际分布','预期分布')
answer = questdlg('实际与预期一致？',['人工判断' type_init],'Y','N','Y');
verifyEqual(testCase,answer,'Y')
end

function test_e_d_u_r(testCase)
% test case 'entire domain uniform random'
simulation=get_simulation('default'); %仿真参数 结构体
type_init='entire domain uniform random';
% dimension_vec=[3,2]; % 目前仅适用于1D
% position_init=get_position_init( simulation, type_init, dimension_vec );
% assert(isequal(dimension_vec,size(position_init)))

dimension_vec=[1e5,1];
position_init=get_position_init( simulation, type_init, dimension_vec );

num=length(position_init);
h_plot1=histogram(position_init);
bin_width=h_plot1.BinWidth;
hold on
n=num*bin_width*1/simulation.Lx;
line([0,simulation.Lx],[n,n],'Color','r','LineWidth',3);
axis([0,simulation.Lx,-inf,inf])
legend('实际分布','预期分布')
answer = questdlg('实际与预期一致？',['人工判断' type_init],'Y','N','Y');
verifyEqual(testCase,answer,'Y')
end

function test_s_r_u_r(testCase)
% test case 'source region uniform random' 
simulation=get_simulation('default'); %仿真参数 结构体
type_init='source region uniform random' ;
% dimension_vec=[3,2]; % 目前仅适用于1D
% position_init=get_position_init( simulation, type_init, dimension_vec );
% assert(isequal(dimension_vec,size(position_init)))

dimension_vec=[1e5,1];
position_init=get_position_init( simulation, type_init, dimension_vec );

num=length(position_init);
h_plot1=histogram(position_init);
bin_width=h_plot1.BinWidth;
hold on
source_region_Lx=simulation.source_region(2)-simulation.source_region(1);
n=num*bin_width*1/source_region_Lx;
line([simulation.source_region(1),simulation.source_region(2)],[n,n],'Color','r','LineWidth',3);
axis([0,simulation.Lx,-inf,inf])
legend('实际分布','预期分布')
answer = questdlg('实际与预期一致？',['人工判断' type_init],'Y','N','Y');
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
figure
end

function teardown(testCase)  % do not change function name
% close figure, for example
close
end
