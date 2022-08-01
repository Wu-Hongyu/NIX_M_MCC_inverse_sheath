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
assert(isequal([3,2],size(v_init))) %ָ��ά��

pdf_maxwellian_velocity=@(v) exp(-v.*v/(2*vth*vth))/(sqrt(2*pi)*vth);
check_pdf(testCase, type_init, 1e8,0.03, pdf_maxwellian_velocity) 
end

function check_pdf(testCase, type_init, total_num,tolerance, pdf_fun)
% ���ʵ��pdf��Ԥ��pdfһ��
vth=1e6;
dimension_vec=[total_num,1];
v_init=get_v_init( vth, type_init, dimension_vec);
num=length(v_init);
% ʹ��histogram���Զ��ֶ�ͳ��
h_plot1=histogram(v_init); %�Լ��ᵥ��figure
bin_width=h_plot1.BinWidth;
N_actual=h_plot1.Values;
x_edges=h_plot1.BinEdges;
pd_excepted=pdf_fun(x_edges);
N_excepted=num*bin_width*(pd_excepted(1:end-1)+pd_excepted(2:end))/2; %PDF�ֶ����λ���

% ������
% ʵ��Ƶ������ʼ����Ƶ����С���ʻ�С����ʱ���������׽ϴ���˲���������
% verifyEqual(testCase,N_actual,N_excepted,'RelTol',5e-2);  %Failed
% ����ͳһ�������������ݲ��˲���ͳһ�ľ������
% ʵ��ʹ������ͼ���ֵA�ɱ����ľ�������AԽ������N*PDF��� ֵԽ��

% TODO������һ������ġ����Ͳ����˼·��ʵ���ϣ���Ӧ����A��Ӧ����N������Ӧ����������

A=max(N_excepted);
verifyEqual(testCase,N_actual,N_excepted,'AbsTol',A*tolerance);

% ������test case��һ��testʱ���������µ�ȫ���жϴ���
% verifyLessThan(testCase,rms(N_actual-N_excepted),A*0.01);
% % ��ͼ���˹��ж�
% hold on
% plot_x=(x_edges(1:end-1)+x_edges(2:end))/2;
% plot(plot_x,N_excepted,'-r','LineWidth',3)
% legend('ʵ�ʷֲ�','Ԥ�ڷֲ�')
% answer = questdlg('ʵ����Ԥ��һ�£�',['�˹��ж�' type_init],'Y','N','Y');
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
