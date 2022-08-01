%% Main function to generate tests
function tests = test_plasma_physic_tools
% test get_plasma_derived_parameters
tests = functiontests(localfunctions);
end

%% Test Functions
function test_basic(testCase)
% test basic
constants=get_constants();
tolerance=1/100;
verify_equal=@(actual, expected) verifyEqual(testCase,actual,expected,'RelTol',tolerance);

n0=4e17;
T0=2;
debye_length=get_debye_length( n0, T0 );
omega_pe= get_omega_pe( n0 );
omega_pHp= get_omega_pi( n0, 1 ,1 );
omega_pHn= get_omega_pi( n0,-1,1 );
vth_e=sqrt(-constants.q_m_ration_e*T0);
N_D=get_N_D(n0, T0);

verify_equal(debye_length, 1.66e-5)
verify_equal(omega_pe/2/pi, 5.679e9)
verify_equal(omega_pHp/2/pi, 1.325e8)
verify_equal(omega_pHn/2/pi, 1.325e8)
verify_equal(vth_e, 5.93e5)
verify_equal(N_D, n0*debye_length^3*4*pi/3)
end

function test_interp_linear_loglog(testCase)
% test interp_linear_loglog
exp_list=1:9;
known_x0=10.^exp_list;
known_y0=2.^flip(exp_list);
% 则该曲线关系为log2(y)=10-log10(x)
known_y1=known_y0;
known_y1([3,4,9])=0;
known_y1(1)=known_y1(2);
known_y1(7)=known_y1(6);
new_exp_list=sort([exp_list,0.5,2.5,3.5,6.5,8.5,9.5]);
new_x=10.^new_exp_list;
% 直角系中线性插值
known_y20=interp1(known_x0,known_y0,new_x,'linear','extrap');
known_y21=interp1(known_x0,known_y1,new_x,'linear','extrap');
% loglog系中线性插值
known_y30=interp_linear_loglog(known_x0,known_y0,new_x);
known_y31=interp_linear_loglog(known_x0,known_y1,new_x);

% 基础插值能力
verifyEqual(testCase,log2(known_y30),10-new_exp_list, 'RelTol', 1e-4)
% 对零值、水平线的插值能力
verifyEqual(testCase,known_y31([5:7,end-1:end]),zeros(1,5))
verifyEqual(testCase,known_y31(end-3),known_y30(end-3))
a=flip(exp_list);
b=a(6)*ones(1,3);
verifyEqual(testCase,log2(known_y31(9:11)),b, 'RelTol', 1e-4)
b=a(2)*ones(1,3);
verifyEqual(testCase,log2(known_y31(1:3)),b, 'RelTol', 1e-4)

% 人工对比查看
% 基础插值能力
figure
subplot(1,2,1)
plot(known_x0,known_y0,'-r','LineWidth',5)
% semilogx
hold on
plot(new_x,known_y20,'--k','LineWidth',4)
plot(new_x,known_y30,'--y')
grid on
legend('origin', '线性坐标中插值', 'log-log中插值')
subplot(1,2,2)
loglog(known_x0,known_y0,'-r','LineWidth',5)
hold on
loglog(new_x,known_y20,'--k','LineWidth',4)
loglog(new_x,known_y30,'--y')
grid on
legend('origin', '线性坐标中插值', 'log-log中插值')
title('基础插值能力')
% 对零值、水平线的插值能力
figure
subplot(1,2,1)
semilogx(known_x0,known_y1,'--b','LineWidth',5)
hold on
semilogx(new_x,known_y21,'--c','LineWidth',4)
semilogx(new_x,known_y31,'--m')
grid on
legend('origin', '线性坐标中插值', 'log-log中插值')
subplot(1,2,2)
loglog(known_x0,known_y1,'-b','LineWidth',5)
hold on
loglog(new_x,known_y21,'--c','LineWidth',4)
loglog(new_x,known_y31,'--m')
legend('origin', '线性坐标中插值', 'log-log中插值')
grid on
title('对零值、水平线的插值能力')

end
% 人工查看：
% 对于loglog坐标系中呈线性的数据，如截面数据
%%%%%% 基础插值能力
% interp_linear_loglog ok
% 应在loglog中插值：
% 1. loglog中插值，能维持loglog坐标系中线性关系
% 2. loglog中插值最小值为0+，而线性坐标插值会出现负数
%%%%%% 对零值、水平线的插值能力
% interp_linear_loglog ok
% 在0的邻域内插值均得到0.
% 所以原始数据中间不应该有0

function test_isequal_with_error(testCase)
% test isequal_with_error
verifyTrue(testCase,isequal_with_error(0,0,'AbsTol',1e-7));
verifyFalse(testCase,isequal_with_error(0,1,'AbsTol',1e-7));
verifyTrue(testCase,isequal_with_error(0,1,'AbsTol',1.1));
verifyFalse(testCase,isequal_with_error(0,1,'AbsTol',1));
verifyTrue(testCase,isequal_with_error(1,1,'RelTol',0.2));
verifyFalse(testCase,isequal_with_error(1.1,1,'RelTol',0.01));
verifyTrue(testCase,isequal_with_error(1.1,1,'RelTol',0.2));
verifyFalse(testCase,isequal_with_error(1.1,1,'RelTol',0.1));
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
