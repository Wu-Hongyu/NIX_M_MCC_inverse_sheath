%% Main function to generate tests
function tests = test_solve_possion_equation
% test solve_possion_equation
tests = functiontests(localfunctions);
end

%% Test Functions
function test_uniform_rho_boundary00(testCase)
% test uniform charge density and boundary type 00 齐次
constants=get_constants();% 全局常数 结构体
simulation=get_simulation('default'); %仿真参数 结构体
simulation.field_boundaries_type=0;%电位边界条件类型
U1=0;
Un=0;
simulation.field_boundaries=[U1,Un];%电位边界条件值
[ u, A, b_extra ] = get_field_init( simulation );
grid_point_x=(0:simulation.dx:simulation.Lx)';
rho0=-simulation.num0_macro_e*simulation.weight*constants.e/simulation.Lx;
rho=rho0*ones(simulation.num_grid_point,1);
b=get_b( rho,b_extra,simulation );
u=get_u( u,A,b,'direct inverse');
% rcond(full(A))=5.0000e-05，A病态，ichol(A) 主元出现负数，不能使用ICCG

coeff_a2=-rho0/(2*constants.eps0);
coeff_a1=rho0*simulation.Lx/2/constants.eps0+(Un-U1)/simulation.Lx;
coeff_a0=U1;
u_analytic_fun=@(x) coeff_a2*x.*x+coeff_a1*x+coeff_a0;
u_analytic=u_analytic_fun(grid_point_x);
flag=sum(abs(u-u_analytic)>5e-2);
verifyEqual(testCase,flag,0)
% 画图，人工判断
figure
plot(grid_point_x,u,'-r','LineWidth',3)
hold on
plot(grid_point_x,u_analytic,'-.b','LineWidth',3)
% line([0,simulation.Lx],[n,n],'Color','r','LineWidth',3);
% axis([0,simulation.Lx,-inf,inf])
legend('实际分布','预期分布')
xlabel('{\itx} [m]')
ylabel('{\it\phi} [V]')
answer = questdlg('实际与预期一致？','人工判断','Y','N','Y');
verifyEqual(testCase,answer,'Y')
close
end

function test_parabolic_rho_boundary00(testCase)
% test parabolic charge density and boundary type 00 非齐次
constants=get_constants();% 全局常数 结构体
simulation=get_simulation('default'); %仿真参数 结构体
simulation.field_boundaries_type=0;%电位边界条件类型
U1=2;
Un=3;
simulation.field_boundaries=[U1,Un];%电位边界条件值
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
% 画图，人工判断
figure
plot(grid_point_x,u,'-r','LineWidth',3)
hold on
plot(grid_point_x,u_analytic,'-.b','LineWidth',3)
% line([0,simulation.Lx],[n,n],'Color','r','LineWidth',3);
% axis([0,simulation.Lx,-inf,inf])
legend('实际分布','预期分布')
xlabel('{\itx} [m]')
ylabel('{\it\phi} [V]')
answer = questdlg('实际与预期一致？','人工判断','Y','N','Y');
verifyEqual(testCase,answer,'Y')
close
end

function test_uniform_rho_boundary01(testCase)
% test uniform charge density and boundary type 01 非齐次
constants=get_constants();% 全局常数 结构体
simulation=get_simulation('default'); %仿真参数 结构体
simulation.field_boundaries_type=1;%电位边界条件类型
U1=2;
Un=3;
simulation.field_boundaries=[U1,Un];%电位边界条件值
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
% 画图，人工判断
figure
plot(grid_point_x,u,'-r','LineWidth',3)
hold on
plot(grid_point_x,u_analytic,'-.b','LineWidth',3)
% line([0,simulation.Lx],[n,n],'Color','r','LineWidth',3);
% axis([0,simulation.Lx,-inf,inf])
legend('实际分布','预期分布')
xlabel('{\itx} [m]')
ylabel('{\it\phi} [V]')
answer = questdlg('实际与预期一致？','人工判断','Y','N','Y');
verifyEqual(testCase,answer,'Y')
close
end

function test_uniform_rho_boundary11(testCase)
% test uniform charge density and boundary type 11 周期
constants=get_constants();% 全局常数 结构体
simulation=get_simulation('default'); %仿真参数 结构体
simulation.field_boundaries_type=3;%电位边界条件类型
simulation.field_boundaries=[2,3];%无意义值
[ u, A, b_extra ] = get_field_init( simulation );
grid_point_x=(0:simulation.dx:simulation.Lx)';
rho0=-simulation.num0_macro_e*simulation.weight*constants.e/simulation.Lx;
rho=rho0*ones(simulation.num_grid_point,1);
b=get_b( rho,b_extra,simulation );
u=get_u( u,A,b,'direct inverse');
% rcond(full(A))=5.0000e-05，A病态，ichol(A) 主元出现负数，不能使用ICCG

U1=0;
Un=0;
coeff_a2=-rho0/(2*constants.eps0);
coeff_a1=rho0*simulation.Lx/2/constants.eps0+(Un-U1)/simulation.Lx;
coeff_a0=U1;
u_analytic_fun=@(x) coeff_a2*x.*x+coeff_a1*x+coeff_a0;
u_analytic=u_analytic_fun(grid_point_x);

flag=sum(abs(u-u_analytic)>5e-2);
verifyEqual(testCase,flag,0)
% 画图，人工判断
figure
plot(grid_point_x,u,'-r','LineWidth',3)
hold on
plot(grid_point_x,u_analytic,'-.b','LineWidth',3)
% line([0,simulation.Lx],[n,n],'Color','r','LineWidth',3);
% axis([0,simulation.Lx,-inf,inf])
legend('实际分布','预期分布')
xlabel('{\itx} [m]')
ylabel('{\it\phi} [V]')
answer = questdlg('实际与预期一致？','人工判断','Y','N','Y');
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
