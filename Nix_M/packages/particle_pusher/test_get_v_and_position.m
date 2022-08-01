%% Main function to generate tests
function tests = test_get_v_and_position
% test get_v_andposition
tests = functiontests(localfunctions);
end

%% Test Functions
% TODO: ��ʱ�����к��������

function test_leap_frog_1D(testCase)
% test leap_frog_1D
% �糡�¼�г�˶�
simulation=get_simulation('default');
num_period=2; 
% �ɸ�Ϊ1e3���۲��ۻ����
% ģ�Ͳ���
L=simulation.Lx;
U0=30;

v0=0; % H��x'=A���޳����ͷ�
position0=L;
constants=get_constants();
q_m_ratio=constants.q_m_ration_H;

% x: λ��
x_sh=@(x) x-L/2; %λ��
E= @(x) -8*U0*x_sh(x)/L/L;
A=L/2;
w0=2*sqrt(2*U0*q_m_ratio)/L;
x_sh_analytic=@(t) A*cos(w0*t); 
T_analytic=2*pi/w0;

% ������
% ��ɢ����
simulation.dx=simulation.Lx/20;
simulation.num_grid_point=simulation.Lx/simulation.dx+1;
simulation.dt=pi*L/sqrt(2*U0*q_m_ratio)/100;
% ѭ����ʼ��
simulation.all_timesteps=num_period*ceil(T_analytic/simulation.dt);
v_record=zeros(simulation.all_timesteps+1,1);
position_record=zeros(simulation.all_timesteps+1,1);

constants=get_constants();
simulation.weight=1;
Ek=@(v,m) simulation.weight*sum(m*v.*v/(2*constants.e)); %���ܣ���λeV
Eke0=sum(Ek(v_record,constants.me)); %��ʼ1D�ܶ���
if 0==Eke0
    Eke0=1;
end

v_sh_analytic=@(t) -A*w0*sin(w0*t); % �ٶȽ�����
v_n_half_timestep_analytic=v_sh_analytic(-simulation.dt/2); %��ʱ��ڵ��ٶȽ�����
% ��ʱ��ڵ� �ٶȳ�ʼ��
% v_record(1)=v0; % ֱ��ʹ����ʱ��ڵ� �ٶ�ֵ
% ������ �ƽ���ʱ��ڵ�
ti=1;
E_particle =E(position_record(ti));
v_n_half_timestep_leapfrog=v0+q_m_ratio*(-simulation.dt/2)*E_particle;
fprintf('v(-0.5��t)= ����%d, ��ֵ%d\n',v_n_half_timestep_analytic,v_n_half_timestep_leapfrog)
% v_record(1)=v_n_half_timestep_leapfrog;
v_record(1)=v_n_half_timestep_analytic;
% v_record(1)=v0;
Eke(1)=sum(Ek(v_record,constants.me));

position_record(1)=position0;
for ti=1:simulation.all_timesteps
E_particle =E(position_record(ti));
[ v_record(ti+1) , position_record(ti+1) ] = pusher1D_leap_frog( v_record(ti) , position_record(ti) , q_m_ratio, E_particle ,simulation );
% [ v, position ] = pusher1D_leap_frog( v, position, q_m_ratio, E_particle ,simulation );

Eke(ti+1)=sum(Ek(v_record(ti+1),constants.me))/Eke0; %��һ��1D�ܶ���
end
plot_timestep=(0:simulation.all_timesteps)';
plot_x_sh_leapfrog=position_record-L/2;
plot_x_sh_analytic=x_sh_analytic(simulation.dt*plot_timestep);

% ������
% ���ʹ���������������㸽����ֵ����ֵС�����������ױ��Ŵ�
% verifyEqual(testCase,plot_x_sh_leapfrog,plot_x_sh_analytic,'RelTol',1/100); %Failed
% ���������С��A/100������A/1000��RMS�������С��A/1000
verifyEqual(testCase,plot_x_sh_leapfrog,plot_x_sh_analytic,'AbsTol',A/100);
verifyLessThan(testCase,rms(plot_x_sh_leapfrog-plot_x_sh_analytic),A/1000);

% % �켣 ��ͼ���˹��ж�
% figure
% plot(plot_timestep, plot_x_sh_leapfrog,'-b','linewidth',3)
% hold on
% plot(plot_timestep, plot_x_sh_analytic,'--r','linewidth',3)
% xlabel('ʱ�� [ʱ��]');
% ylabel('x [m]');
% % title('��г�˶�');
% L1=legend('������','������');
% set(L1,'location','northeast');
% set(L1,'AutoUpdate','off');
% line([T_analytic/simulation.dt,T_analytic/simulation.dt],[-A,A],'linestyle',':','linewidth',3,'color','k');
% answer = questdlg('��ֵ���������һ�£�','�˹��ж�','Y','N','Y');
% verifyEqual(testCase,answer,'Y')
% close

% % ���ܷ���
% % 20201127 ���� ����
% figure
% yyaxis left;
% h_fig1=plot(plot_timestep, plot_x_sh_analytic,'-b','linewidth',3);
% xlabel('ʱ�� [ʱ��]');
% ylabel('x [m]');
% yyaxis right;
% h_fig2=plot(plot_timestep, Eke,'--r','linewidth',3);
% 
% % ����
% % v_sh_analytic
% 
% ylabel('Ek [eV]');
% % title('��г�˶�');
% L1=legend([h_fig1,h_fig2],'λ��','����');
% set(L1,'location','northeast');
% set(L1,'AutoUpdate','off');
% % line([T_analytic/simulation.dt,T_analytic/simulation.dt],[-A,A],'linestyle',':','linewidth',3,'color','k');
% % close

cumulative_phase_error=simulation.all_timesteps*(w0*simulation.dt)^3/24; %�ۻ����
cumulative_timestep_error=(T_analytic/simulation.dt)*cumulative_phase_error/2/pi;
fprintf('%d �����ۻ����= %d ����= %d ʱ�� \n',simulation.all_timesteps,cumulative_phase_error/2/pi,cumulative_timestep_error)
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