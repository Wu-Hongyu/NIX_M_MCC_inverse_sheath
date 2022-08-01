%% Main function to generate tests
function tests = test_get_simulation
% test get_simulation
tests = functiontests(localfunctions);
end

%% Test Functions
function test_default_basic(testCase)
% test case 'default' %A1DPIC_Ver1.0.m, 2019Montellano，不check
actual_simulation=get_simulation('default');
goal_simulation.dt=10*10^-12;%时间步长
goal_simulation.end_time=6E-6;%仿真时间总长
goal_simulation.all_timesteps=floor(goal_simulation.end_time/goal_simulation.dt);%总循环次数
goal_simulation.Lx=0.01;%仿真区域长度
goal_simulation.source_region=[0.2*goal_simulation.Lx,0.3*goal_simulation.Lx];
goal_simulation.num_grid_point=201;%网格格点数目
goal_simulation.dx=goal_simulation.Lx/(goal_simulation.num_grid_point-1);%空间步长
goal_simulation.ne0=1E16;%等离子体密度
goal_simulation.Te0=1;%单位eV
goal_simulation.field_boundaries_type=0;%电位边界条件类型
goal_simulation.field_boundaries=[0,0];%电位边界条件值
goal_simulation.B0=0; %代表性磁感应强度 T
goal_simulation.particle_boundaries_type=0;%粒子边界条件类型
verifyEqual(testCase,actual_simulation,goal_simulation)
end

function test_auto1_basic(testCase)
% test case 'auto1' %输入ne, Te, Lx, end_time，不check
simulation1=get_simulation('auto1');
simulation2=get_simulation('auto1','ne0',1E16,'Te0',1,'end_time',6E-6,'Lx',0.01);
verifyEqual(testCase,simulation1,simulation2)
end

function test_check_default_parameters(testCase)
% test check default simulation parameters
simulation=get_simulation('default');
result=check_simulation_parameters( simulation, 3 );
verifyTrue(testCase,result.flag_check)
end

function test_check_BATMAN_old_MF(testCase)
% test check default simulation parameters
simulation=get_simulation('BATMAN_old_MF');
result=check_simulation_parameters( simulation, 2 );
verifyTrue(testCase,result.flag_check)
end

function test_check_auto0_parameters(testCase)
% test check simulation parameters from auto1 type
simulation=get_simulation('auto0','ne0',1E17,'Te0',10,'end_time',6E-6,'Lx',0.01);
result=check_simulation_parameters( simulation, 3 );
verifyTrue(testCase,result.flag_check)
end

function test_check_auto1_parameters(testCase)
% test check simulation parameters from auto1 type
simulation=get_simulation('auto1','ne0',1E17,'Te0',10,'end_time',6E-6,'Lx',0.01);
result=check_simulation_parameters( simulation, 2 );
verifyTrue(testCase,result.flag_check)
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
