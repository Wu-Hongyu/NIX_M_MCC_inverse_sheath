clear
close all;

%文件路径
addpath('./packages')
% 全局变量
constants=get_constants();% 全局常数 结构体

%--------仿真参数-----------------------------------------------------------------------------------
simulation=get_simulation('default'); %仿真参数 结构体
Te=2;
n0=1E15;
potential=get_sheath_potential(Te);
x_final=(0:1:200)*simulation.dx;
figure
plot(x_final,potential,'-b','LineWidth',3)%多个步长电压平均值
axis([0,simulation.Lx,-inf,inf])
hold on
ub=-Te*log(sqrt(constants.mH/(2*pi*constants.me)));
uW=ones(size(potential))*ub;
plot(x_final,uW,'-r','LineWidth',1)%多个步长电压平均值
axis([0,simulation.Lx,-inf,inf])
