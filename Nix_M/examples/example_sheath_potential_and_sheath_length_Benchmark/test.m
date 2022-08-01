clear
close all;

%�ļ�·��
addpath('./packages')
% ȫ�ֱ���
constants=get_constants();% ȫ�ֳ��� �ṹ��

%--------�������-----------------------------------------------------------------------------------
simulation=get_simulation('default'); %������� �ṹ��
Te=2;
n0=1E15;
potential=get_sheath_potential(Te);
x_final=(0:1:200)*simulation.dx;
figure
plot(x_final,potential,'-b','LineWidth',3)%���������ѹƽ��ֵ
axis([0,simulation.Lx,-inf,inf])
hold on
ub=-Te*log(sqrt(constants.mH/(2*pi*constants.me)));
uW=ones(size(potential))*ub;
plot(x_final,uW,'-r','LineWidth',1)%���������ѹƽ��ֵ
axis([0,simulation.Lx,-inf,inf])
