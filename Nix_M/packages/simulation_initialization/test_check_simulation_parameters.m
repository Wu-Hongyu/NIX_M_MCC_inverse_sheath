%% Main function to generate tests
function tests = test_check_simulation_parameters
% test check_simulation_parameters
tests = functiontests(localfunctions);
end

%% Test Functions
function test_check_simulation_default(testCase)
% test check_discrete_parameters ok
simulation=get_simulation('default');
result_default=check_simulation_parameters( simulation, 3 );
verifyTrue(testCase,result_default.flag_check)
end
function test_check_simulation_with_B_warn(testCase)
% test check_discrete_parameters for simulation_with_B
simulation=get_simulation('BATMAN_old_MF');
result=check_simulation_parameters( simulation, 2 );
verifyTrue(testCase,result.flag_check)
end
function test_unbound_Lx_s_debye_ratio(testCase)
% test unbounded plasma and Lx < λ_De
simulation=get_simulation('default');
simulation.field_boundaries_type=3;
simulation.Lx=0.9*get_debye_length(simulation.ne0,simulation.Te0);
ecpected_e=get_exception('warn','Lx < λ_De，无法识别德拜球');
check_exception_during_check_simulation(testCase,simulation,3,ecpected_e)
end
function test_bound_Lx_s_debye_ratio(testCase)
% test bounded plasma and Lx < λ_De
simulation=get_simulation('default');
simulation.field_boundaries_type=0;
simulation.Lx=0.9*get_debye_length(simulation.ne0,simulation.Te0);
ecpected_e=get_exception('error','Lx < λ_De，等离子体判据要求Lx>>λ_De');
check_exception_during_check_simulation(testCase,simulation,3,ecpected_e)
end
function test_bound_Lx_s_10debye_ratio(testCase)
% test bounded plasma and Lx < 10*λ_De
simulation=get_simulation('default');
simulation.field_boundaries_type=0;
simulation.Lx=9*get_debye_length(simulation.ne0,simulation.Te0);%导致异常
ecpected_e=get_exception('warn','λ_De <= Lx < 10*λ_De，等离子体判据要求Lx>>λ_De');
check_exception_during_check_simulation(testCase,simulation,3,ecpected_e)
end
function test_N_D_s1(testCase)
% test if result.N_D<1
constants=get_constants();
simulation=get_simulation('default');
N_D=0.9; 
simulation.Te0=(3*N_D/4/pi)^(2/3)*simulation.ne0^(1/3)*constants.e/constants.eps0;%导致异常
simulation.Lx=10*get_debye_length(simulation.ne0,simulation.Te0); %规避已测试过的异常
ecpected_e=get_exception('error','N_D < 1，等离子体判据要求N_D>>>1');
check_exception_during_check_simulation(testCase,simulation,3,ecpected_e)
end
function test_N_D_s100(testCase)
% test elseif result.N_D<100
constants=get_constants();
simulation=get_simulation('default');
N_D=99; 
simulation.Te0=(3*N_D/4/pi)^(2/3)*simulation.ne0^(1/3)*constants.e/constants.eps0;%导致异常
simulation.Lx=10*get_debye_length(simulation.ne0,simulation.Te0); %规避已测试过的异常
ecpected_e=get_exception('warn','1 <= N_D < 100，等离子体判据要求N_D>>>1');
check_exception_during_check_simulation(testCase,simulation,3,ecpected_e)
end
% 暂未实现其他check的test

function test_type2(testCase)
% test ignore warn in bound_Lx_s_10debye_ratio
simulation=get_simulation('default');
simulation.field_boundaries_type=0;
simulation.Lx=9*get_debye_length(simulation.ne0,simulation.Te0);%导致异常
ecpected_e=get_exception('empty','');
check_exception_during_check_simulation(testCase,simulation,2,ecpected_e)
end

function test_type1(testCase)
% test ignore error in bound_Lx_s_debye_ratio
simulation=get_simulation('default');
simulation.field_boundaries_type=0;
simulation.Lx=0.9*get_debye_length(simulation.ne0,simulation.Te0);
ecpected_e=get_exception('empty','');
check_exception_during_check_simulation(testCase,simulation,1,ecpected_e)
end
%% Aid function
function check_exception_during_check_simulation(testCase,simulation,type,ecpected_e)
actual_e=get_exception('empty','');
try
check_simulation_parameters( simulation, type );
catch actual_e
end
verify_exception_equal(testCase,actual_e,ecpected_e)
end

function verify_exception_equal(testCase,actual_e,excepted_e)
verifyEqual(testCase,actual_e.identifier,excepted_e.identifier)
% verifyEqual(testCase,actual_e.message,excepted_e.message)
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
