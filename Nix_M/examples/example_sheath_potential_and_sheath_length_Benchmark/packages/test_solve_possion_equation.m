%% Main function to generate tests
function tests = test_solve_possion_equation
% test solve_possion_equation
tests = functiontests(localfunctions);
end

%% Test Functions
function test_uniform_rho_boundary00(testCase)
% test uniform charge density and boundary type 00 ���
constants=get_constants();% ȫ�ֳ��� �ṹ��
simulation=get_simulation('default'); %������� �ṹ��
simulation.field_boundaries_type=0;%��λ�߽���������
U1=0;
Un=0;
simulation.field_boundaries=[U1,Un];%��λ�߽�����ֵ
[ u, A, b_extra ] = get_field_init( simulation );
grid_point_x=(0:simulation.dx:simulation.Lx)';
rho0=-simulation.num0_macro_e*simulation.weight*constants.e/simulation.Lx;
rho=rho0*ones(simulation.num_grid_point,1);
b=get_b( rho,b_extra,simulation );
u=get_u( u,A,b,'direct inverse');
% rcond(full(A))=5.0000e-05��A��̬��ichol(A) ��Ԫ���ָ���������ʹ��ICCG

coeff_a2=-rho0/(2*constants.eps0);
coeff_a1=rho0*simulation.Lx/2/constants.eps0+(Un-U1)/simulation.Lx;
coeff_a0=U1;
u_analytic_fun=@(x) coeff_a2*x.*x+coeff_a1*x+coeff_a0;
u_analytic=u_analytic_fun(grid_point_x);
flag=sum(abs(u-u_analytic)>5e-2);
verifyEqual(testCase,flag,0)
% ��ͼ���˹��ж�
figure
plot(grid_point_x,u,'-r','LineWidth',3)
hold on
plot(grid_point_x,u_analytic,'-.b','LineWidth',3)
% line([0,simulation.Lx],[n,n],'Color','r','LineWidth',3);
% axis([0,simulation.Lx,-inf,inf])
legend('ʵ�ʷֲ�','Ԥ�ڷֲ�')
xlabel('{\itx} [m]')
ylabel('{\it\phi} [V]')
answer = questdlg('ʵ����Ԥ��һ�£�','�˹��ж�','Y','N','Y');
verifyEqual(testCase,answer,'Y')
close
end

function test_parabolic_rho_boundary00(testCase)
% test parabolic charge density and boundary type 00 �����
constants=get_constants();% ȫ�ֳ��� �ṹ��
simulation=get_simulation('default'); %������� �ṹ��
simulation.field_boundaries_type=0;%��λ�߽���������
U1=2;
Un=3;
simulation.field_boundaries=[U1,Un];%��λ�߽�����ֵ
[ u, A, b_extra ] = get_field_init( simulation );
grid_point_x=(0:simulation.dx:simulation.Lx)';
rho0=2*simulation.num0_macro_H*simulation.weight*constants.e/simulation.Lx;
rho=rho0*grid_point_x.*(grid_point_x-simulation.Lx);
b=get_b( rho,b_extra,simulation );
u=get_u( u,A,b,'direct inverse');

coeff_a4=-rho0/(12*constants.eps0);
coeff_a3=rho0*simulation.Lx/(6*constants.eps0);
coeff_a1=-rho0*simulation.Lx^3/12/constants.eps0+(Un-U1)/simulation.Lx;
coeff_a0=U1;
u_analytic_fun=@(x) coeff_a4*x.^4+coeff_a3*x.^3+coeff_a1*x+coeff_a0;
u_analytic=u_analytic_fun(grid_point_x);

flag=sum(abs(u-u_analytic)>5e-2);
verifyEqual(testCase,flag,0)
% ��ͼ���˹��ж�
figure
plot(grid_point_x,u,'-r','LineWidth',3)
hold on
plot(grid_point_x,u_analytic,'-.b','LineWidth',3)
% line([0,simulation.Lx],[n,n],'Color','r','LineWidth',3);
% axis([0,simulation.Lx,-inf,inf])
legend('ʵ�ʷֲ�','Ԥ�ڷֲ�')
xlabel('{\itx} [m]')
ylabel('{\it\phi} [V]')
answer = questdlg('ʵ����Ԥ��һ�£�','�˹��ж�','Y','N','Y');
verifyEqual(testCase,answer,'Y')
close
end

function test_uniform_rho_boundary01(testCase)
% test uniform charge density and boundary type 01 �����
constants=get_constants();% ȫ�ֳ��� �ṹ��
simulation=get_simulation('default'); %������� �ṹ��
simulation.field_boundaries_type=1;%��λ�߽���������
U1=2;
Un=3;
simulation.field_boundaries=[U1,Un];%��λ�߽�����ֵ
[ u, A, b_extra ] = get_field_init( simulation );
grid_point_x=(0:simulation.dx:simulation.Lx)';
rho0=simulation.num0_macro_e*simulation.weight*constants.e/simulation.Lx;
rho=rho0*ones(simulation.num_grid_point,1);
b=get_b( rho,b_extra,simulation );
u=get_u( u,A,b,'direct inverse');

coeff_a2=-rho0/(2*constants.eps0);
coeff_a1=rho0*simulation.Lx/constants.eps0+Un;
coeff_a0=U1;
u_analytic_fun=@(x) coeff_a2*x.*x+coeff_a1*x+coeff_a0;
u_analytic=u_analytic_fun(grid_point_x);

flag=sum(abs(u-u_analytic)>5e-2);
verifyEqual(testCase,flag,0)
% ��ͼ���˹��ж�
figure
plot(grid_point_x,u,'-r','LineWidth',3)
hold on
plot(grid_point_x,u_analytic,'-.b','LineWidth',3)
% line([0,simulation.Lx],[n,n],'Color','r','LineWidth',3);
% axis([0,simulation.Lx,-inf,inf])
legend('ʵ�ʷֲ�','Ԥ�ڷֲ�')
xlabel('{\itx} [m]')
ylabel('{\it\phi} [V]')
answer = questdlg('ʵ����Ԥ��һ�£�','�˹��ж�','Y','N','Y');
verifyEqual(testCase,answer,'Y')
close
end

function test_uniform_rho_boundary11(testCase)
% test uniform charge density and boundary type 11 ����
constants=get_constants();% ȫ�ֳ��� �ṹ��
simulation=get_simulation('default'); %������� �ṹ��
simulation.field_boundaries_type=3;%��λ�߽���������
simulation.field_boundaries=[2,3];%������ֵ
[ u, A, b_extra ] = get_field_init( simulation );
grid_point_x=(0:simulation.dx:simulation.Lx)';
rho0=-simulation.num0_macro_e*simulation.weight*constants.e/simulation.Lx;
rho=rho0*ones(simulation.num_grid_point,1);
b=get_b( rho,b_extra,simulation );
u=get_u( u,A,b,'direct inverse');
% rcond(full(A))=5.0000e-05��A��̬��ichol(A) ��Ԫ���ָ���������ʹ��ICCG

U1=0;
Un=0;
coeff_a2=-rho0/(2*constants.eps0);
coeff_a1=rho0*simulation.Lx/2/constants.eps0+(Un-U1)/simulation.Lx;
coeff_a0=U1;
u_analytic_fun=@(x) coeff_a2*x.*x+coeff_a1*x+coeff_a0;
u_analytic=u_analytic_fun(grid_point_x);

flag=sum(abs(u-u_analytic)>5e-2);
verifyEqual(testCase,flag,0)
% ��ͼ���˹��ж�
figure
plot(grid_point_x,u,'-r','LineWidth',3)
hold on
plot(grid_point_x,u_analytic,'-.b','LineWidth',3)
% line([0,simulation.Lx],[n,n],'Color','r','LineWidth',3);
% axis([0,simulation.Lx,-inf,inf])
legend('ʵ�ʷֲ�','Ԥ�ڷֲ�')
xlabel('{\itx} [m]')
ylabel('{\it\phi} [V]')
answer = questdlg('ʵ����Ԥ��һ�£�','�˹��ж�','Y','N','Y');
verifyEqual(testCase,answer,'Y')
close
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
