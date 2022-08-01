%% Main function to generate tests
function tests = test_get_v_and_position
% test get_v_andposition
tests = functiontests(localfunctions);
end

%% Test Functions
% TODO: 长时间运行后的误差分析

function test_leap_frog_1D(testCase)
% test leap_frog_1D
% 电场下简谐运动
simulation=get_simulation('default');
num_period=2; 
% 可改为1e3，观察累积误差
% 模型参数
L=simulation.Lx;
U0=30;

v0=0; % H从x'=A处无初速释放
position0=L;
constants=get_constants();
q_m_ratio=constants.q_m_ration_H;

% x: 位置
x_sh=@(x) x-L/2; %位移
E= @(x) -8*U0*x_sh(x)/L/L;
A=L/2;
w0=2*sqrt(2*U0*q_m_ratio)/L;
x_sh_analytic=@(t) A*cos(w0*t); 
T_analytic=2*pi/w0;

% 蛙跳法
% 离散限制
simulation.dx=simulation.Lx/20;
simulation.num_grid_point=simulation.Lx/simulation.dx+1;
simulation.dt=pi*L/sqrt(2*U0*q_m_ratio)/100;
% 循环初始化
simulation.all_timesteps=num_period*ceil(T_analytic/simulation.dt);
v_record=zeros(simulation.all_timesteps+1,1);
position_record=zeros(simulation.all_timesteps+1,1);

constants=get_constants();
simulation.weight=1;
Ek=@(v,m) simulation.weight*sum(m*v.*v/(2*constants.e)); %动能，单位eV
Eke0=sum(Ek(v_record,constants.me)); %初始1D总动能
if 0==Eke0
    Eke0=1;
end

v_sh_analytic=@(t) -A*w0*sin(w0*t); % 速度解析解
v_n_half_timestep_analytic=v_sh_analytic(-simulation.dt/2); %半时间节点速度解析解
% 半时间节点 速度初始化
% v_record(1)=v0; % 直接使用整时间节点 速度值
% 蛙跳法 推进半时间节点
ti=1;
E_particle =E(position_record(ti));
v_n_half_timestep_leapfrog=v0+q_m_ratio*(-simulation.dt/2)*E_particle;
fprintf('v(-0.5Δt)= 解析%d, 数值%d\n',v_n_half_timestep_analytic,v_n_half_timestep_leapfrog)
% v_record(1)=v_n_half_timestep_leapfrog;
v_record(1)=v_n_half_timestep_analytic;
% v_record(1)=v0;
Eke(1)=sum(Ek(v_record,constants.me));

position_record(1)=position0;
for ti=1:simulation.all_timesteps
E_particle =E(position_record(ti));
[ v_record(ti+1) , position_record(ti+1) ] = pusher1D_leap_frog( v_record(ti) , position_record(ti) , q_m_ratio, E_particle ,simulation );
% [ v, position ] = pusher1D_leap_frog( v, position, q_m_ratio, E_particle ,simulation );

Eke(ti+1)=sum(Ek(v_record(ti+1),constants.me))/Eke0; %归一化1D总动能
end
plot_timestep=(0:simulation.all_timesteps)';
plot_x_sh_leapfrog=position_record-L/2;
plot_x_sh_analytic=x_sh_analytic(simulation.dt*plot_timestep);

% 误差分析
% 如果使用相对误差，由于在零点附近真值绝对值小，相对误差容易被放大
% verifyEqual(testCase,plot_x_sh_leapfrog,plot_x_sh_analytic,'RelTol',1/100); %Failed
% 最大绝对误差小于A/100，大于A/1000；RMS绝对误差小于A/1000
verifyEqual(testCase,plot_x_sh_leapfrog,plot_x_sh_analytic,'AbsTol',A/100);
verifyLessThan(testCase,rms(plot_x_sh_leapfrog-plot_x_sh_analytic),A/1000);

% % 轨迹 画图，人工判断
% figure
% plot(plot_timestep, plot_x_sh_leapfrog,'-b','linewidth',3)
% hold on
% plot(plot_timestep, plot_x_sh_analytic,'--r','linewidth',3)
% xlabel('时间 [时步]');
% ylabel('x [m]');
% % title('简谐运动');
% L1=legend('蛙跳法','解析解');
% set(L1,'location','northeast');
% set(L1,'AutoUpdate','off');
% line([T_analytic/simulation.dt,T_analytic/simulation.dt],[-A,A],'linestyle',':','linewidth',3,'color','k');
% answer = questdlg('数值解与解析解一致？','人工判断','Y','N','Y');
% verifyEqual(testCase,answer,'Y')
% close

% % 动能分析
% % 20201127 陈鹏 搁置
% figure
% yyaxis left;
% h_fig1=plot(plot_timestep, plot_x_sh_analytic,'-b','linewidth',3);
% xlabel('时间 [时步]');
% ylabel('x [m]');
% yyaxis right;
% h_fig2=plot(plot_timestep, Eke,'--r','linewidth',3);
% 
% % 解析
% % v_sh_analytic
% 
% ylabel('Ek [eV]');
% % title('简谐运动');
% L1=legend([h_fig1,h_fig2],'位移','动能');
% set(L1,'location','northeast');
% set(L1,'AutoUpdate','off');
% % line([T_analytic/simulation.dt,T_analytic/simulation.dt],[-A,A],'linestyle',':','linewidth',3,'color','k');
% % close

cumulative_phase_error=simulation.all_timesteps*(w0*simulation.dt)^3/24; %累积相差
cumulative_timestep_error=(T_analytic/simulation.dt)*cumulative_phase_error/2/pi;
fprintf('%d 步后，累积误差= %d 周期= %d 时步 \n',simulation.all_timesteps,cumulative_phase_error/2/pi,cumulative_timestep_error)
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